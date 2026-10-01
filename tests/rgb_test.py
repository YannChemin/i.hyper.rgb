"""Tests of i.hyper.rgb."""

import os
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest

import grass.script as gs
from grass.exceptions import CalledModuleError
from grass.experimental import TemporaryMapsetSession
from grass.script import array as garray
from grass.tools import Tools

# Bands every 10 nm, 400-900 nm; band k holds the value k + 1 plus a
# per-pixel offset, so each channel tells which bands it took.
WAVES = np.arange(400.0, 901.0, 10.0)
ROWS, COLS = 3, 4
OFFSET = np.arange(ROWS * COLS, dtype=float).reshape(ROWS, COLS) / 100


def write_cube(ds, name="cube", invalid=()):
    """Cube of len(WAVES) bands with i.hyper metadata; bands at the
    wavelengths in invalid are marked invalid."""
    gs.run_command(
        "g.region", n=ROWS, s=0, w=0, e=COLS, res=1, t=len(WAVES), b=0,
        res3=1, tbres=1, env=ds.env,
    )  # fmt: skip
    cube = garray.array3d(dtype=np.float32, env=ds.env)
    for k in range(len(WAVES)):
        cube[k] = k + 1 + OFFSET
    cube.write(name, overwrite=True)
    validity = [bool(w not in invalid) for w in WAVES]
    code = (
        "import sys\n"
        "from grass.script.utils import get_lib_path\n"
        "sys.path.append(get_lib_path(modname='i_hyper_lib', libname='hyper_meta'))\n"
        "from hyper_meta import HyperMetadata\n"
        f"HyperMetadata(n_bands={len(WAVES)}, wavelengths={WAVES.tolist()}, "
        f"validity={validity}).save('{name}')\n"
    )
    subprocess.run([sys.executable, "-c", code], env=ds.env, check=True)


@pytest.fixture
def ds(tmp_path):
    """Tools in a temporary mapset of an XY project."""
    gs.create_project(tmp_path / "xy")
    with (
        gs.setup.init(tmp_path / "xy", env=os.environ.copy()) as session,
        TemporaryMapsetSession(env=session.env) as mapset,
        Tools(session=mapset) as tools,
    ):
        yield SimpleNamespace(tools=tools, env=mapset.env)


def read(ds, name):
    return np.asarray(garray.array(name, env=ds.env))


def band_value(wave):
    return float(np.nonzero(WAVES == wave)[0][0] + 1) + OFFSET


def test_nearest_bands(ds):
    """Each channel is the band nearest to its wavelength; the channels make
    an image group; the input region is left unchanged."""
    write_cube(ds)
    gs.run_command("g.region", n=10, s=0, w=0, e=10, res=1, env=ds.env)
    ds.tools.i_hyper_rgb(
        input="cube", output="c", red_wavelength=652, green_wavelength=548,
        blue_wavelength=451,
    )  # fmt: skip
    region = gs.region(env=ds.env)
    assert (region["rows"], region["cols"]) == (10, 10)
    gs.run_command("g.region", raster="c_red", env=ds.env)
    np.testing.assert_allclose(read(ds, "c_red"), band_value(650), rtol=1e-6)
    np.testing.assert_allclose(read(ds, "c_green"), band_value(550), rtol=1e-6)
    np.testing.assert_allclose(read(ds, "c_blue"), band_value(450), rtol=1e-6)
    group = gs.read_command("i.group", group="c_rgb", flags="g", env=ds.env)
    assert [m.split("@")[0] for m in group.split()] == ["c_red", "c_green", "c_blue"]
    title = gs.raster_info("c_red", env=ds.env)["title"].strip('"')
    assert title == "red channel: 650.0 nm"


def test_invalid_bands_skipped(ds):
    """Bands marked invalid in the metadata are never chosen."""
    write_cube(ds, invalid=(650.0,))
    ds.tools.i_hyper_rgb(input="cube", output="c", red_wavelength=651)
    gs.run_command("g.region", raster="c_red", env=ds.env)
    np.testing.assert_allclose(read(ds, "c_red"), band_value(660), rtol=1e-6)


@pytest.mark.parametrize(
    ("statistic", "expected"),
    [("mean", 0.0), ("median", 0.0), ("min", -2.0), ("max", 2.0),
     ("sd1_pos", np.sqrt(2.0)), ("sd2_neg", -2 * np.sqrt(2.0))],
)  # fmt: skip
def test_window_statistics(ds, statistic, expected):
    """bandwidth=40 takes the five bands within 20 nm (values v - 2 .. v +
    2 around the centre band's v)."""
    write_cube(ds)
    ds.tools.i_hyper_rgb(
        input="cube", output="c", bandwidth=40, statistic=statistic,
        red_wavelength=650,
    )  # fmt: skip
    gs.run_command("g.region", raster="c_red", env=ds.env)
    np.testing.assert_allclose(
        read(ds, "c_red"), band_value(650) + expected, rtol=1e-5, atol=1e-5
    )
    assert "of 5 bands, 630.0-670.0 nm" in gs.raster_info("c_red", env=ds.env)["title"]


def test_rescale_and_colorblind(ds):
    """-n rescales each channel to 0-255; colorblind mixes the channels
    with its matrix."""
    write_cube(ds)
    ds.tools.i_hyper_rgb(input="cube", output="n", flags="n")
    gs.run_command("g.region", raster="n_red", env=ds.env)
    red, green, blue = (read(ds, f"n_{c}") for c in ("red", "green", "blue"))
    for channel in (red, green, blue):
        assert channel.min() == pytest.approx(0.0, abs=1e-9)
        assert channel.max() == pytest.approx(255.0)
    ds.tools.i_hyper_rgb(input="cube", output="m", flags="n", colorblind="protanopia")
    np.testing.assert_allclose(read(ds, "m_red"), 0.567 * red + 0.433 * green)
    np.testing.assert_allclose(read(ds, "m_blue"), 0.242 * green + 0.758 * blue)


def test_cmyk(ds):
    """colorspace=cmyk makes four channels and a _cmyk group."""
    write_cube(ds)
    ds.tools.i_hyper_rgb(input="cube", output="k", colorspace="cmyk")
    gs.run_command("g.region", raster="k_key", env=ds.env)
    np.testing.assert_allclose(read(ds, "k_key"), band_value(800), rtol=1e-6)
    group = gs.read_command("i.group", group="k_cmyk", flags="g", env=ds.env)
    assert len(group.split()) == 4


def test_errors(ds):
    """No band in the window, colorblind with CMYK and a cube without
    metadata are refused."""
    write_cube(ds)
    with pytest.raises(CalledModuleError, match="No valid band within"):
        ds.tools.i_hyper_rgb(
            input="cube", output="e", bandwidth=4, red_wavelength=1500
        )
    with pytest.raises(CalledModuleError, match="colorspace=rgb only"):
        ds.tools.i_hyper_rgb(
            input="cube", output="e", colorspace="cmyk", colorblind="tritanopia"
        )
    gs.run_command("r3.mapcalc", expression="bare = 1.0", env=ds.env)
    with pytest.raises(CalledModuleError, match="metadata"):
        ds.tools.i_hyper_rgb(input="bare", output="e")
