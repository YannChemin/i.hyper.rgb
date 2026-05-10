"""
Tests that libras3d is functional for i.hyper.rgb in standalone mode.

Run without GRASS:
    pytest testsuite/test_ras3d_rgb.py -v
"""
import os, sys, tempfile
import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from test_ras3d_common import (
    WYVERN_PATH, skip_without_ras3d, skip_without_wyvern,
    open_cube_checked, assert_band_valid, install_ras3d_shim, make_wl_sidecar,
)


@skip_without_ras3d
def test_shim_installs():
    """ras3d_grass_shim.install() populates sys.modules['grass.script']."""
    install_ras3d_shim()
    import grass.script as gs
    assert hasattr(gs, 'parser')
    assert hasattr(gs, 'fatal')
    assert hasattr(gs, 'raster3d_info')


@skip_without_ras3d
@skip_without_wyvern
def test_open_cube_geotiff():
    """open_cube() on Wyvern GeoTIFF returns correct dims."""
    import ras3d
    h, r = open_cube_checked(WYVERN_PATH)
    assert r['cols']   == 6003
    assert r['rows']   == 7825
    assert r['depths'] == 23
    ras3d.close_cube(h)


@skip_without_ras3d
@skip_without_wyvern
def test_get_band_geotiff():
    """get_band() returns a valid float32 array for each of the 23 Wyvern bands."""
    import ras3d
    h, r = open_cube_checked(WYVERN_PATH)
    for z in (0, 11, 22):
        arr = ras3d.get_band(h, z)
        assert arr.shape == (r['rows'], r['cols'])
        assert_band_valid(arr, f'Wyvern band {z}')
    ras3d.close_cube(h)


@skip_without_ras3d
@skip_without_wyvern
def test_extract_z_slice_writes_geotiff(tmp_path):
    """extract_z_slice() in ras3d mode writes a GeoTIFF and populates the cache."""
    install_ras3d_shim()
    os.environ['RAS3D_OUTDIR'] = str(tmp_path)

    sys.path.insert(0, '/home/yann/dev/i.hyper.rgb')
    import importlib, i_hyper_rgb
    importlib.invalidate_caches()

    # Directly call the module's extract_z_slice in ras3d mode
    from i_hyper_rgb import extract_z_slice as _ext   # noqa: F401 — import side-effect test
    slice_name = 'test_rgb_band1'
    _ext(WYVERN_PATH, 1, slice_name)

    from ras3d_grass_shim import get_band_cache
    assert slice_name in get_band_cache()
    arr = get_band_cache()[slice_name]
    assert_band_valid(arr, 'extract_z_slice band 1')

    out = tmp_path / (slice_name + '.tif')
    assert out.exists(), f"Expected {out} to be written"


@skip_without_ras3d
@skip_without_wyvern
def test_wavelength_sidecar(tmp_path):
    """get_band_wavelengths() reads the .wl.json sidecar in ras3d mode."""
    install_ras3d_shim()
    import ras3d
    h, r = open_cube_checked(WYVERN_PATH)
    sidecar, wl_list = make_wl_sidecar(WYVERN_PATH, r['depths'])
    ras3d.close_cube(h)

    sys.path.insert(0, '/home/yann/dev/i.hyper.rgb')
    from i_hyper_rgb import get_band_wavelengths
    wl_dict = get_band_wavelengths(WYVERN_PATH)
    assert len(wl_dict) == r['depths']
    os.unlink(sidecar)
