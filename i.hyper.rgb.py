#!/usr/bin/env python3

##############################################################################
# MODULE:    i.hyper.rgb
# AUTHOR(S): Yann Chemin, after i.hyper.composite by Alen Mangafic and
#            Tomaz Zagar (Geodetic Institute of Slovenia)
# PURPOSE:   Create RGB or CMYK composites from a hyperspectral 3D raster,
#            each channel the nearest band to its wavelength or a statistic
#            of the bands around it.
# COPYRIGHT: (C) 2025 by the GRASS Development Team
# SPDX-License-Identifier: GPL-2.0-or-later
##############################################################################

# %module
# % description: Creates RGB or CMYK composites from a hyperspectral 3D raster map.
# % keyword: imagery
# % keyword: hyperspectral
# % keyword: composite
# %end

# %option G_OPT_R3_INPUT
# % key: input
# % description: Input hyperspectral 3D raster map (with i.hyper metadata)
# % guisection: Input
# %end

# %option
# % key: output
# % type: string
# % required: yes
# % description: Base name of the channel maps (<output>_red, ...) and of the image group (<output>_rgb or <output>_cmyk)
# % guisection: Output
# %end

# %option
# % key: colorspace
# % type: string
# % options: rgb,cmyk
# % answer: rgb
# % description: Output colour space
# % guisection: Output
# %end

# %option
# % key: bandwidth
# % type: double
# % answer: 0
# % description: Width (nm) of the window of bands around each channel's wavelength; 0 for the nearest band only
# % guisection: Processing
# %end

# %option
# % key: statistic
# % type: string
# % options: mean,median,mode,min,max,sd1_pos,sd2_pos,sd3_pos,sd1_neg,sd2_neg,sd3_neg
# % answer: mean
# % description: Statistic of the bands within bandwidth combined into a channel
# % descriptions: mean;mean;median;median;mode;most frequent value;min;minimum;max;maximum;sd1_pos;mean plus one standard deviation;sd2_pos;mean plus two standard deviations;sd3_pos;mean plus three standard deviations;sd1_neg;mean minus one standard deviation;sd2_neg;mean minus two standard deviations;sd3_neg;mean minus three standard deviations
# % guisection: Processing
# %end

# %option
# % key: colorblind
# % type: string
# % options: none,protanopia,deuteranopia,tritanopia
# % answer: none
# % description: Mix the RGB channels with a dichromacy simulation matrix (the composite as seen with that colour vision deficiency)
# % guisection: Accessibility
# %end

# %option
# % key: red_wavelength
# % type: double
# % answer: 650
# % description: Wavelength of the red channel (nm)
# % guisection: Wavelengths
# %end

# %option
# % key: green_wavelength
# % type: double
# % answer: 550
# % description: Wavelength of the green channel (nm)
# % guisection: Wavelengths
# %end

# %option
# % key: blue_wavelength
# % type: double
# % answer: 450
# % description: Wavelength of the blue channel (nm)
# % guisection: Wavelengths
# %end

# %option
# % key: cyan_wavelength
# % type: double
# % answer: 490
# % description: Wavelength of the cyan channel (nm), CMYK only
# % guisection: Wavelengths
# %end

# %option
# % key: magenta_wavelength
# % type: double
# % answer: 580
# % description: Wavelength of the magenta channel (nm), CMYK only
# % guisection: Wavelengths
# %end

# %option
# % key: yellow_wavelength
# % type: double
# % answer: 570
# % description: Wavelength of the yellow channel (nm), CMYK only
# % guisection: Wavelengths
# %end

# %option
# % key: key_wavelength
# % type: double
# % answer: 800
# % description: Wavelength of the key (black) channel (nm), CMYK only
# % guisection: Wavelengths
# %end

# %flag
# % key: n
# % description: Rescale each channel linearly to 0-255
# % guisection: Processing
# %end

import builtins
import importlib.util
import json
import os
import sys
import uuid

import numpy as np

# Without GRASS, run on files through libras3d when it is installed: its
# shim stands in for grass.script.
RAS3D = False
if not os.environ.get("GISBASE"):
    try:
        if importlib.util.find_spec("ras3d") and importlib.util.find_spec(
            "ras3d_grass_shim"
        ):
            from ras3d_grass_shim import install as _ras3d_install

            _ras3d_install()
            RAS3D = True
    except ImportError:
        pass

import grass.script as gs  # noqa: E402

if not hasattr(builtins, "_"):  # no GRASS translation catalogue (libras3d)
    builtins._ = str

CHANNELS = {
    "rgb": ("red", "green", "blue"),
    "cmyk": ("cyan", "magenta", "yellow", "key"),
}

# Dichromacy simulation matrices (rows: output red, green, blue).
COLORBLIND = {
    "protanopia": ((0.567, 0.433, 0), (0.558, 0.442, 0), (0, 0.242, 0.758)),
    "deuteranopia": ((0.625, 0.375, 0), (0.7, 0.3, 0), (0, 0.3, 0.7)),
    "tritanopia": ((0.95, 0.05, 0), (0, 0.433, 0.567), (0, 0.475, 0.525)),
}

R_SERIES = {
    "mean": "average",
    "median": "median",
    "mode": "mode",
    "min": "minimum",
    "max": "maximum",
}


def select_bands(target, bandwidth, wavelengths):
    """Indices of the valid bands of a channel: within bandwidth / 2 of the
    target, or the nearest one when bandwidth is 0."""
    wavelengths = np.asarray(wavelengths, dtype=float)
    if bandwidth > 0:
        inside = np.flatnonzero(np.abs(wavelengths - target) <= bandwidth / 2)
        if inside.size == 0:
            gs.fatal(
                _(
                    "No valid band within {width} nm of {target} nm (bands "
                    "{low:.1f}-{high:.1f} nm)"
                ).format(
                    width=bandwidth / 2,
                    target=target,
                    low=wavelengths.min(),
                    high=wavelengths.max(),
                )
            )
        return inside.tolist()
    return [int(np.argmin(np.abs(wavelengths - target)))]


def sd_factor(statistic):
    """Signed number of standard deviations of the sd* statistics."""
    return int(statistic[2]) * (1 if statistic.endswith("_pos") else -1)


class GrassBackend:
    """Channels as 2D raster maps of the current mapset. The current region
    must be the cube's (main() sets a temporary one)."""

    def __init__(self, cube):
        found = gs.find_file(cube, element="raster_3d")
        if not found["name"]:
            gs.fatal(_("3D raster map <{}> not found").format(cube))
        self.cube = cube
        self.name, self.mapset = found["name"], found["mapset"]
        self.depths = int(gs.raster3d_info(cube)["depths"])
        self.tmp = f"tmp_ihrgb_{uuid.uuid4().hex[:8]}"
        self.temporary = []
        self.slicer = self._load_slicer()

    def wavelengths(self):
        """Wavelengths and depth indices of the valid bands, from the i.hyper
        metadata (i_hyper_lib.hyper_meta)."""
        from grass.script.utils import get_lib_path

        path = get_lib_path(modname="i_hyper_lib", libname="hyper_meta")
        if not path:
            gs.fatal(_("i_hyper_lib (i.hyper suite) is not installed"))
        if path not in sys.path:
            sys.path.append(path)
        import hyper_meta

        try:
            meta = hyper_meta.HyperMetadata.load(self.name, self.mapset)
            axis = meta.resolve_band_axis(self.depths)
        except (OSError, ValueError) as error:
            gs.fatal(
                _("Cannot read the metadata of <{map}>: {error}").format(
                    map=self.cube, error=error
                )
            )
        wavelengths = axis.get("wavelengths")
        if wavelengths is None:
            gs.fatal(_("<{}> has no band wavelengths").format(self.cube))
        valid = np.flatnonzero(axis["validity"])
        depth = [int(axis["source_to_depth"][i]) for i in valid]
        keep = [k for k, d in enumerate(depth) if d >= 0]
        wavelengths = np.asarray(wavelengths, dtype=float)[valid][keep]
        if np.isnan(wavelengths).any():
            gs.fatal(_("<{}> has valid bands without wavelength").format(self.cube))
        return wavelengths, [depth[k] for k in keep]

    @staticmethod
    def _load_slicer():
        """Rast3d_extract_z_slice() of the GRASS raster3d library, which
        reads the tiles of one depth once, when this GRASS has it (through
        the grass.lib ctypes bindings); None otherwise."""
        try:
            import grass.lib.gis as libgis
            import grass.lib.raster3d as libraster3d

            slicer = libraster3d.Rast3d_extract_z_slice
        except (ImportError, OSError, AttributeError):
            gs.verbose(_("No Rast3d_extract_z_slice(): bands read with r3.to.rast"))
            return None
        libgis.G_gisinit("i.hyper.rgb")
        gs.verbose(_("Bands read with Rast3d_extract_z_slice()"))
        return slicer

    def band(self, depth):
        """Band at a depth (0-based) as a temporary 2D map."""
        name = f"{self.tmp}_z{depth}"
        if name in self.temporary:
            return name
        if self.slicer is not None:
            ret = self.slicer(
                self.name.encode(), self.mapset.encode(), depth, name.encode()
            )
            if ret != 0:
                gs.fatal(
                    _("Cannot extract band {band} of <{map}>").format(
                        band=depth + 1, map=self.cube
                    )
                )
        else:
            # One depth of the cube's region: r3.to.rast reads only it.
            gs.run_command("g.region", b=depth, t=depth + 1, tbres=1)
            gs.run_command("r3.to.rast", input=self.cube, output=name, quiet=True)
            gs.run_command("g.region", raster_3d=self.cube)
            gs.run_command("g.rename", raster=(f"{name}_00001", name), quiet=True)
        self.temporary.append(name)
        return name

    def combine(self, name, depths, statistic):
        """Channel map name from the bands at depths."""
        maps = [self.band(d) for d in depths]
        if len(maps) == 1:
            gs.mapcalc(f"{name} = {maps[0]}", quiet=True)
        elif statistic in R_SERIES:
            gs.run_command(
                "r.series",
                input=maps,
                output=name,
                method=R_SERIES[statistic],
                quiet=True,
            )
        else:
            mean, sd = f"{self.tmp}_mean", f"{self.tmp}_sd"
            gs.run_command(
                "r.series",
                input=maps,
                output=(mean, sd),
                method=("average", "stddev"),
                overwrite=True,
                quiet=True,
            )
            self.temporary += [mean, sd]
            gs.mapcalc(f"{name} = {mean} + ({sd_factor(statistic)}) * {sd}", quiet=True)

    def rescale(self, name):
        """Linear rescaling of a channel to 0-255."""
        info = gs.raster_info(name)
        low, high = info["min"], info["max"]
        if low is None or high is None or not high > low:
            gs.fatal(_("Channel <{}> is constant: cannot rescale it").format(name))
        tmp = f"{self.tmp}_rescale"
        gs.mapcalc(f"{tmp} = 255.0 * ({name} - {low}) / ({high} - {low})", quiet=True)
        gs.run_command("g.rename", raster=(tmp, name), overwrite=True, quiet=True)

    def mix(self, names, matrix):
        """Replace the channels by their mix with a 3 x 3 matrix."""
        tmp = [f"{self.tmp}_mix{k}" for k in range(3)]
        for out, row in zip(tmp, matrix):
            expr = " + ".join(f"{w} * {n}" for w, n in zip(row, names) if w)
            gs.mapcalc(f"{out} = {expr}", quiet=True)
        for out, name in zip(tmp, names):
            gs.run_command("g.rename", raster=(out, name), overwrite=True, quiet=True)

    def finish(self, names, group, description):
        for name in names:
            gs.run_command("r.support", map=name, title=description[name])
            gs.raster_history(name, overwrite=True)
        gs.run_command("i.group", group=group, subgroup=group, input=names)
        gs.message(_("Image group <{}> created").format(group))

    def cleanup(self):
        if self.temporary:
            gs.run_command(
                "g.remove", type="raster", name=self.temporary, flags="f", quiet=True
            )


class Ras3dBackend:
    """Without GRASS: channels as arrays, written as GeoTIFF files to
    $RAS3D_OUTDIR by libras3d."""

    def __init__(self, cube):
        import ras3d

        self.ras3d = ras3d
        self.cube = cube
        self.handle = ras3d.open_cube(cube)
        self.depths = ras3d.get_region(self.handle)["depths"]
        self.arrays = {}
        self.bands = {}

    def wavelengths(self):
        """Wavelengths from the <cube>.wl.json sidecar (µm or nm)."""
        path = os.path.splitext(self.cube)[0] + ".wl.json"
        if not os.path.exists(path):
            gs.fatal(_("No wavelength sidecar {} for the cube").format(path))
        with open(path, encoding="utf-8") as stream:
            values = [float(w) for w in json.load(stream)]
        if len(values) != self.depths:
            gs.fatal(
                _("{path} has {n} wavelengths for {depths} bands").format(
                    path=path, n=len(values), depths=self.depths
                )
            )
        nm = np.array([w * 1000 if w < 10 else w for w in values])
        return nm, list(range(self.depths))

    def band(self, depth):
        if depth not in self.bands:
            self.bands[depth] = self.ras3d.get_band(self.handle, depth).astype(
                np.float64
            )
        return self.bands[depth]

    def combine(self, name, depths, statistic):
        stack = np.stack([self.band(d) for d in depths])
        if len(depths) == 1:
            out = stack[0]
        elif statistic == "mean":
            out = np.nanmean(stack, axis=0)
        elif statistic == "median":
            out = np.nanmedian(stack, axis=0)
        elif statistic == "min":
            out = np.nanmin(stack, axis=0)
        elif statistic == "max":
            out = np.nanmax(stack, axis=0)
        elif statistic == "mode":
            out = np.apply_along_axis(
                lambda v: np.unique(v, return_counts=True)[0][
                    np.argmax(np.unique(v, return_counts=True)[1])
                ],
                0,
                stack,
            )
        else:
            out = np.nanmean(stack, axis=0) + sd_factor(statistic) * np.nanstd(
                stack, axis=0
            )
        self.arrays[name] = out

    def rescale(self, name):
        a = self.arrays[name]
        low, high = np.nanmin(a), np.nanmax(a)
        if not high > low:
            gs.fatal(_("Channel {} is constant: cannot rescale it").format(name))
        self.arrays[name] = 255.0 * (a - low) / (high - low)

    def mix(self, names, matrix):
        old = [self.arrays[n] for n in names]
        for name, row in zip(names, matrix):
            self.arrays[name] = sum(w * a for w, a in zip(row, old))

    def finish(self, names, group, description):
        import ras3d_write

        for name in names:
            path = ras3d_write.outpath(name)
            ras3d_write.write_raster2d(
                path, self.arrays[name].astype(np.float32), self.handle
            )
            gs.message(_("{name}: {path}").format(name=description[name], path=path))

    def cleanup(self):
        self.ras3d.close_cube(self.handle)


def main(options, flags):
    space = options["colorspace"]
    bandwidth = float(options["bandwidth"])
    statistic = options["statistic"]
    colorblind = options["colorblind"]
    if bandwidth < 0:
        gs.fatal(_("bandwidth must be 0 or more"))
    if colorblind != "none" and space != "rgb":
        gs.fatal(_("colorblind= applies to colorspace=rgb only"))

    if RAS3D:
        backend = Ras3dBackend(options["input"])
    else:
        # The cube's region, for its bands and the channel maps.
        gs.use_temp_region()
        gs.run_command("g.region", raster_3d=options["input"])
        backend = GrassBackend(options["input"])
    try:
        wavelengths, depths = backend.wavelengths()
        gs.verbose(_("{} valid bands").format(len(wavelengths)))
        names, description = [], {}
        for channel in CHANNELS[space]:
            target = float(options[f"{channel}_wavelength"])
            chosen = select_bands(target, bandwidth, wavelengths)
            name = f"{options['output']}_{channel}"
            backend.combine(name, [depths[k] for k in chosen], statistic)
            if len(chosen) == 1:
                description[name] = _("{channel} channel: {wave:.1f} nm").format(
                    channel=channel, wave=wavelengths[chosen[0]]
                )
            else:
                description[name] = _(
                    "{channel} channel: {statistic} of {n} bands, "
                    "{low:.1f}-{high:.1f} nm"
                ).format(
                    channel=channel,
                    statistic=statistic,
                    n=len(chosen),
                    low=wavelengths[chosen[0]],
                    high=wavelengths[chosen[-1]],
                )
            gs.message(description[name])
            names.append(name)
        if flags["n"]:
            for name in names:
                backend.rescale(name)
        if colorblind != "none":
            backend.mix(names, COLORBLIND[colorblind])
            gs.message(_("Channels mixed for {}").format(colorblind))
        backend.finish(names, f"{options['output']}_{space}", description)
    finally:
        backend.cleanup()


if __name__ == "__main__":
    sys.exit(main(*gs.parser()))
