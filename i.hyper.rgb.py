#!/usr/bin/env python
# ── ras3d standalone detection ────────────────────────────────────────────────
import os as _os
_RAS3D = False
if not _os.environ.get('GISBASE'):
    try:
        import importlib.util as _ilu
        if _ilu.find_spec('ras3d') and _ilu.find_spec('ras3d_grass_shim'):
            from ras3d_grass_shim import install as _r3_install
            _r3_install()
            _RAS3D = True
    except Exception:
        pass
# ─────────────────────────────────────────────────────────────────────────────
##############################################################################
# MODULE:    i.hyper.rgb
# AUTHOR(S): Created for hyperspectral RGB/CMYK composite generation
# PURPOSE:   Create RGB/CMYK composites from hyperspectral imagery
# COPYRIGHT: (C) 2025 by the GRASS Development Team
# SPDX-License-Identifier: GPL-2.0-or-later
##############################################################################

# %module
# % description: Create RGB/CMYK composites from hyperspectral imagery with statistical band selection
# % keyword: imagery
# % keyword: hyperspectral
# % keyword: composite
# %end

# %option G_OPT_R3_INPUT
# % key: input
# % required: yes
# % description: Input hyperspectral 3D raster map (from i.hyper.import)
# % guisection: Input
# %end

# %option
# % key: output
# % type: string
# % required: yes
# % description: Base name for output image group (will append _rgb or _cmyk)
# % guisection: Output
# %end

# %option
# % key: colorspace
# % type: string
# % required: yes
# % options: rgb,cmyk
# % answer: rgb
# % description: Output color space
# % guisection: Output
# %end

# %option
# % key: statistic
# % type: string
# % required: no
# % options: mean,median,mode,min,max,sd1_pos,sd2_pos,sd3_pos,sd1_neg,sd2_neg,sd3_neg
# % answer: mean
# % description: Statistical method for wavelength band selection
# % guisection: Processing
# %end

# %option
# % key: colorblind
# % type: string
# % required: no
# % options: none,protanopia,deuteranopia,tritanopia
# % answer: none
# % description: Color blind safe palette adjustment
# % guisection: Accessibility
# %end

# %option
# % key: red_wavelength
# % type: double
# % required: no
# % answer: 650
# % description: Target wavelength for red channel (nm)
# % guisection: Wavelengths
# %end

# %option
# % key: green_wavelength
# % type: double
# % required: no
# % answer: 550
# % description: Target wavelength for green channel (nm)
# % guisection: Wavelengths
# %end

# %option
# % key: blue_wavelength
# % type: double
# % required: no
# % answer: 450
# % description: Target wavelength for blue channel (nm)
# % guisection: Wavelengths
# %end

# %option
# % key: cyan_wavelength
# % type: double
# % required: no
# % answer: 490
# % description: Target wavelength for cyan channel (nm) - CMYK only
# % guisection: Wavelengths
# %end

# %option
# % key: magenta_wavelength
# % type: double
# % required: no
# % answer: 580
# % description: Target wavelength for magenta channel (nm) - CMYK only
# % guisection: Wavelengths
# %end

# %option
# % key: yellow_wavelength
# % type: double
# % required: no
# % answer: 570
# % description: Target wavelength for yellow channel (nm) - CMYK only
# % guisection: Wavelengths
# %end

# %option
# % key: key_wavelength
# % type: double
# % required: no
# % answer: 800
# % description: Target wavelength for key/black channel (nm) - CMYK only
# % guisection: Wavelengths
# %end

# %flag
# % key: n
# % description: Normalize output bands to 0-255
# % guisection: Processing
# %end

import sys
import os
import re
import ctypes
import grass.script as gs
import numpy as np


# ---------------------------------------------------------------------------
# Fast Z-slice extraction via Rast3d_extract_z_slice() (ctypes)
# ---------------------------------------------------------------------------

_raster3d_lib = None


def _load_raster3d_lib():
    """Load libgrass_raster3d and its deps via ctypes (once per process)."""
    global _raster3d_lib
    if _raster3d_lib is not None:
        return _raster3d_lib

    gisbase = os.environ["GISBASE"]
    libdir = os.path.join(gisbase, "lib")

    for name in ("libgrass_gis.so", "libgrass_raster.so"):
        ctypes.CDLL(os.path.join(libdir, name), ctypes.RTLD_GLOBAL)

    lib = ctypes.CDLL(
        os.path.join(libdir, "libgrass_raster3d.so"), ctypes.RTLD_GLOBAL
    )
    lib.Rast3d_extract_z_slice.restype = ctypes.c_int
    lib.Rast3d_extract_z_slice.argtypes = [
        ctypes.c_char_p,  # name3d
        ctypes.c_char_p,  # mapset3d ("" = search)
        ctypes.c_int,     # z  (0-based)
        ctypes.c_char_p,  # name2d
    ]

    libgis = ctypes.CDLL(os.path.join(libdir, "libgrass_gis.so"))
    libgis.G_gisinit(b"i.hyper.rgb")

    _raster3d_lib = lib
    return lib


def extract_z_slice(name3d, band_num_1based, name2d):
    """Extract one band (1-based) from a 3D raster to a 2D raster.

    Uses Rast3d_extract_z_slice() which opens the map with RASTER3D_NO_CACHE
    and calls Rast3d_get_block() for tile-bulk reads — each tile is loaded
    exactly once instead of spawning r3.to.rast per band.
    """
    if _RAS3D:
        import ras3d as _r3, ras3d_write as _r3w
        _h = _r3.open_cube(name3d)
        _arr = _r3.get_band(_h, band_num_1based - 1)
        from ras3d_grass_shim import get_band_cache
        get_band_cache()[name2d] = _arr
        _r3w.write_raster2d(_r3w.outpath(name2d), _arr, _h)
        _r3.close_cube(_h)
        return
    lib = _load_raster3d_lib()
    z = band_num_1based - 1  # convert 1-based band to 0-based z index
    ret = lib.Rast3d_extract_z_slice(
        name3d.encode(), b"", ctypes.c_int(z), name2d.encode()
    )
    if ret != 0:
        gs.fatal(
            f"Rast3d_extract_z_slice failed for band {band_num_1based} of {name3d}"
        )


def get_band_wavelengths(raster3d):
    """Extract wavelength metadata from Raster3D history via r3.info -h.

    Parses lines of the form 'Band N: WL nm' written by i.hyper.atcorr
    and i.hyper.import into the map's history file.
    """
    if _RAS3D:
        import json as _json
        for _sfx in ('', '.tif', '.tiff', '.h5', '.hdf5'):
            _base = raster3d.removesuffix(_sfx) if raster3d.endswith(_sfx) else raster3d
            _wlp = _base + '.wl.json'
            if _os.path.exists(_wlp):
                with open(_wlp) as _f:
                    _wl = _json.load(_f)
                _wl_nm = [w * 1000 if w < 10 else w for w in _wl]
                return {float(wl): i + 1 for i, wl in enumerate(_wl_nm)}
        import ras3d as _r3
        _h = _r3.open_cube(raster3d); _r = _r3.get_region(_h); _r3.close_cube(_h)
        return {float(i): i for i in range(1, _r['depths'] + 1)}
    try:
        header = gs.read_command('r3.info', map=raster3d, flags='h')
    except Exception as e:
        gs.fatal(f"Cannot read r3.info for {raster3d}: {e}")

    wavelengths = {}
    for line in header.splitlines():
        m = re.match(r'\s*Band\s+(\d+):\s+([\d.]+)\s*nm', line)
        if m:
            band_num = int(m.group(1))
            wl_nm = float(m.group(2))
            wavelengths[wl_nm] = band_num

    if not wavelengths:
        gs.warning("No wavelength metadata found in r3.info -h. "
                   "Using band numbers as wavelengths.")
        info = gs.parse_command('r3.info', map=raster3d, flags='g')
        depths = int(info['depths'])
        wavelengths = {float(i): i for i in range(1, depths + 1)}

    return wavelengths


def find_closest_band(target_wavelength, wavelengths):
    """Find the band closest to target wavelength"""
    closest_wl = min(wavelengths.keys(), key=lambda x: abs(x - target_wavelength))
    return wavelengths[closest_wl], closest_wl


def calculate_statistic(raster3d, band_indices, statistic):
    """Calculate statistical composite from multiple bands"""
    temp_maps = []
    
    for idx in band_indices:
        band_name = f"{raster3d}#{idx}"
        temp_maps.append(band_name)
    
    if statistic == "mean":
        expr = f"({' + '.join(temp_maps)}) / {len(temp_maps)}"
    elif statistic == "median":
        # Use r.series for median
        return temp_maps, "median"
    elif statistic == "mode":
        return temp_maps, "mode"
    elif statistic == "min":
        return temp_maps, "minimum"
    elif statistic == "max":
        return temp_maps, "maximum"
    elif statistic.startswith("sd"):
        # Standard deviation variants
        return temp_maps, statistic
    else:
        gs.fatal(f"Unknown statistic: {statistic}")
    
    return expr, None


def apply_colorblind_adjustment(channels, colorblind_type):
    """Apply colorblind-safe adjustments to RGB channels"""
    adjustments = {
        'protanopia': {  # Red-blind
            'red': (0.567, 0.433, 0),
            'green': (0.558, 0.442, 0),
            'blue': (0, 0.242, 0.758)
        },
        'deuteranopia': {  # Green-blind
            'red': (0.625, 0.375, 0),
            'green': (0.7, 0.3, 0),
            'blue': (0, 0.3, 0.7)
        },
        'tritanopia': {  # Blue-blind
            'red': (0.95, 0.05, 0),
            'green': (0, 0.433, 0.567),
            'blue': (0, 0.475, 0.525)
        }
    }
    
    if colorblind_type not in adjustments:
        return channels
    
    adj = adjustments[colorblind_type]
    gs.message(f"Applying {colorblind_type} color adjustments...")
    
    # Apply transformation matrix
    # This would require actual raster math operations
    # For now, return original channels with warning
    gs.warning(f"Colorblind adjustment for {colorblind_type} selected but requires post-processing")
    
    return channels


def create_rgb_composite(options, flags):
    """Create RGB composite from hyperspectral data"""
    input_raster = options['input']
    output_base = options['output']
    statistic = options['statistic']
    colorblind = options['colorblind']
    normalize = flags['n']
    
    # Get wavelengths
    wavelengths = get_band_wavelengths(input_raster)
    gs.message(f"Found {len(wavelengths)} wavelength bands")
    
    # Find bands for RGB
    red_wl = float(options['red_wavelength'])
    green_wl = float(options['green_wavelength'])
    blue_wl = float(options['blue_wavelength'])
    
    red_band, actual_red = find_closest_band(red_wl, wavelengths)
    green_band, actual_green = find_closest_band(green_wl, wavelengths)
    blue_band, actual_blue = find_closest_band(blue_wl, wavelengths)
    
    gs.message(f"Selected bands - R: {actual_red}nm (band {red_band}), "
               f"G: {actual_green}nm (band {green_band}), "
               f"B: {actual_blue}nm (band {blue_band})")
    
    # Create output rasters for each channel
    channels = {'red': red_band, 'green': green_band, 'blue': blue_band}
    output_maps = []
    
    for color, band in channels.items():
        output_map = f"{output_base}_{color}"

        gs.message(f"Creating {color} channel: {output_map}")

        extract_z_slice(input_raster, band, output_map)

        # Normalize if requested
        if normalize:
            gs.run_command('r.rescale', input=output_map, output=output_map,
                          to='0,255', overwrite=True, quiet=True)

        output_maps.append(output_map)

    # Apply colorblind adjustments if requested
    if colorblind != 'none':
        apply_colorblind_adjustment(output_maps, colorblind)

    # Create image group
    group_name = f"{output_base}_rgb"
    gs.run_command('i.group', group=group_name, subgroup=group_name,
                  input=','.join(output_maps))
    
    gs.message(f"Created RGB image group: {group_name}")
    gs.message(f"Individual bands: {', '.join(output_maps)}")


def create_cmyk_composite(options, flags):
    """Create CMYK composite from hyperspectral data"""
    input_raster = options['input']
    output_base = options['output']
    statistic = options['statistic']
    normalize = flags['n']
    
    # Get wavelengths
    wavelengths = get_band_wavelengths(input_raster)
    gs.message(f"Found {len(wavelengths)} wavelength bands")
    
    # Find bands for CMYK
    cyan_wl = float(options['cyan_wavelength'])
    magenta_wl = float(options['magenta_wavelength'])
    yellow_wl = float(options['yellow_wavelength'])
    key_wl = float(options['key_wavelength'])
    
    cyan_band, actual_cyan = find_closest_band(cyan_wl, wavelengths)
    magenta_band, actual_magenta = find_closest_band(magenta_wl, wavelengths)
    yellow_band, actual_yellow = find_closest_band(yellow_wl, wavelengths)
    key_band, actual_key = find_closest_band(key_wl, wavelengths)
    
    gs.message(f"Selected bands - C: {actual_cyan}nm (band {cyan_band}), "
               f"M: {actual_magenta}nm (band {magenta_band}), "
               f"Y: {actual_yellow}nm (band {yellow_band}), "
               f"K: {actual_key}nm (band {key_band})")
    
    # Create output rasters for each channel
    channels = {
        'cyan': cyan_band,
        'magenta': magenta_band,
        'yellow': yellow_band,
        'key': key_band
    }
    output_maps = []
    
    for color, band in channels.items():
        output_map = f"{output_base}_{color}"

        gs.message(f"Creating {color} channel: {output_map}")

        extract_z_slice(input_raster, band, output_map)

        # Normalize if requested
        if normalize:
            gs.run_command('r.rescale', input=output_map, output=output_map,
                          to='0,255', overwrite=True, quiet=True)

        output_maps.append(output_map)

    # Create image group
    group_name = f"{output_base}_cmyk"
    gs.run_command('i.group', group=group_name, subgroup=group_name,
                  input=','.join(output_maps))
    
    gs.message(f"Created CMYK image group: {group_name}")
    gs.message(f"Individual bands: {', '.join(output_maps)}")


def main(options, flags):
    """Main function"""
    colorspace = options['colorspace']
    
    gs.message(f"Creating {colorspace.upper()} composite...")
    
    if colorspace == 'rgb':
        create_rgb_composite(options, flags)
    elif colorspace == 'cmyk':
        create_cmyk_composite(options, flags)
    else:
        gs.fatal(f"Unknown colorspace: {colorspace}")
    
    gs.message("Composite creation complete!")


if __name__ == "__main__":
    options, flags = gs.parser()
    sys.exit(main(options, flags))
