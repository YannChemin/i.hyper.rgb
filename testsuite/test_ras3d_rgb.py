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


def load_module():
    """i.hyper.rgb.py as a module, in libras3d mode (no GISBASE)."""
    import importlib.util
    install_ras3d_shim()
    os.environ.pop('GISBASE', None)
    path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'i.hyper.rgb.py')
    spec = importlib.util.spec_from_file_location('i_hyper_rgb', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.RAS3D
    return module


@skip_without_ras3d
@skip_without_wyvern
def test_composite_writes_geotiff(tmp_path):
    """In libras3d mode, the RGB channels are written as GeoTIFF files to
    $RAS3D_OUTDIR, each the band nearest to its wavelength."""
    import ras3d
    module = load_module()
    os.environ['RAS3D_OUTDIR'] = str(tmp_path)
    h, r = open_cube_checked(WYVERN_PATH)
    sidecar, wl_list = make_wl_sidecar(WYVERN_PATH, r['depths'])
    ras3d.close_cube(h)
    try:
        options = {
            'input': WYVERN_PATH, 'output': 'rgbtest', 'colorspace': 'rgb',
            'bandwidth': '0', 'statistic': 'mean', 'colorblind': 'none',
            'red_wavelength': '600', 'green_wavelength': '500',
            'blue_wavelength': '400',
        }
        module.main(options, {'n': False})
    finally:
        os.unlink(sidecar)
    for channel in ('red', 'green', 'blue'):
        out = tmp_path / f'rgbtest_{channel}.tif'
        assert out.exists(), f"Expected {out} to be written"
    backend = module.Ras3dBackend(WYVERN_PATH)
    assert_band_valid(backend.band(0), 'band 1')
    backend.cleanup()


@skip_without_ras3d
@skip_without_wyvern
def test_wavelength_sidecar(tmp_path):
    """The libras3d backend reads the .wl.json sidecar."""
    import ras3d
    module = load_module()
    h, r = open_cube_checked(WYVERN_PATH)
    sidecar, wl_list = make_wl_sidecar(WYVERN_PATH, r['depths'])
    ras3d.close_cube(h)
    try:
        backend = module.Ras3dBackend(WYVERN_PATH)
        waves, depths = backend.wavelengths()
        backend.cleanup()
    finally:
        os.unlink(sidecar)
    assert len(waves) == r['depths']
    assert list(waves) == wl_list
    assert depths == list(range(r['depths']))
