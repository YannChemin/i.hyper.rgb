## DESCRIPTION

*i.hyper.rgb* creates RGB or CMYK composites from a hyperspectral 3D
raster map (`raster_3d`) of the *i.hyper* module family. Each colour
channel is made from the band nearest to the channel's wavelength, or
from a statistic of the bands within a window around it.

The band wavelengths and their validity come from the map's *i.hyper*
metadata (`hyper.json`, read with the *i.hyper* library `i_hyper_lib`),
as written by [i.hyper.import](i.hyper.import.md),
[i.hyper.preproc](i.hyper.preproc.md) or
[i.hyper.metadata](i.hyper.metadata.md). Bands marked invalid (water
absorption bands, for example) are never used. A map without this
metadata is refused.

With **bandwidth** 0 (the default), each channel is the valid band
nearest to its wavelength. With a **bandwidth** in nm, each channel
combines all the valid bands within half the bandwidth on either side of
its wavelength, with the **statistic**:

- **mean**, **median**, **mode** (most frequent value), **min**,
  **max** of the bands;
- **sd1_pos**, **sd2_pos**, **sd3_pos**: the mean plus 1, 2 or 3
  standard deviations of the bands;
- **sd1_neg**, **sd2_neg**, **sd3_neg**: the mean minus 1, 2 or 3
  standard deviations.

A channel without any valid band in its window is an error.

The **-n** flag rescales each channel linearly to 0-255 (its minimum to
0, its maximum to 255). Without it, the channels keep the cube's values
(radiance or reflectance).

With **colorblind** (RGB only), the three channels are mixed with the
dichromacy simulation matrix of protanopia, deuteranopia or tritanopia
(Viénot et al. 1999): the composite then shows how it appears with that
colour vision deficiency, so that a map can be checked for readability.
Use it with **-n**, so that the channels have the same scale.

The output is one raster map per channel (`<output>_red`, `_green`,
`_blue`, or `_cyan`, `_magenta`, `_yellow`, `_key`), titled with the
wavelengths used, and an image group of them (`<output>_rgb` or
`<output>_cmyk`).

## NOTES

The bands are read one depth at a time, never the whole cube. When the
GRASS raster3d library provides `Rast3d_extract_z_slice()` (reached
through the `grass.lib.raster3d` bindings), each band is read with it:
the tiles of that depth are read once, in bulk. Otherwise, each band is
extracted by *r3.to.rast* from a temporary 3D region limited to that
depth. Both give the same maps. The computational region of the session
is left unchanged.

Without GRASS (no `GISBASE`), when the
[libras3d](https://github.com/yannchemin/libras3d) Python package is
installed, the module runs on a GeoTIFF or HDF5 cube file instead
(`input=<file>`). The wavelengths then come from the `<cube>.wl.json`
sidecar (in µm or nm), every band is valid, and each channel is written
as a GeoTIFF file to `$RAS3D_OUTDIR`; no image group is made.

Default wavelengths:

- **Red** 650 nm, **green** 550 nm, **blue** 450 nm.
- **Cyan** 490 nm, **magenta** 580 nm, **yellow** 570 nm, **key
  (black)** 800 nm (near-infrared, for contrast).

## EXAMPLES

A true-colour composite of the nearest bands, rescaled for display:

```sh
i.hyper.rgb input=enmap output=enmap_display -n
d.rgb red=enmap_display_red green=enmap_display_green \
    blue=enmap_display_blue
```

Each channel the mean of the bands within 20 nm of its wavelength, with
custom wavelengths:

```sh
i.hyper.rgb input=prisma output=prisma_mean bandwidth=40 statistic=mean \
    red_wavelength=660 green_wavelength=560 blue_wavelength=470
```

The same composite as seen with deuteranopia:

```sh
i.hyper.rgb input=tanager output=tanager_deuteranopia -n \
    colorblind=deuteranopia
```

A CMYK composite with custom wavelengths:

```sh
i.hyper.rgb input=prisma output=prisma_cmyk colorspace=cmyk \
    cyan_wavelength=480 magenta_wavelength=590 yellow_wavelength=560 \
    key_wavelength=850
```

Without GRASS, on a GeoTIFF cube with libras3d:

```sh
RAS3D_OUTDIR=/tmp/out i.hyper.rgb input=/data/wyvern.tif output=wyvern
```

## SEE ALSO

*[i.hyper.composite](i.hyper.composite.md),
[i.hyper.explore](i.hyper.explore.md),
[i.hyper.export](i.hyper.export.md),
[i.hyper.import](i.hyper.import.md),
[i.hyper.metadata](i.hyper.metadata.md),
[i.hyper.preproc](i.hyper.preproc.md),
[d.rgb](https://grass.osgeo.org/grass-stable/manuals/d.rgb.html),
[i.group](https://grass.osgeo.org/grass-stable/manuals/i.group.html),
[r3.to.rast](https://grass.osgeo.org/grass-stable/manuals/r3.to.rast.html)*

## REFERENCES

- Viénot, F., Brettel, H., & Mollon, J. D. (1999). Digital video
  colourmaps for checking the legibility of displays by dichromats.
  *Color Research & Application*, 24(4), 243-252.

## AUTHORS

Yann Chemin, after *i.hyper.composite* by Alen Mangafić and Tomaž Žagar,
Geodetic Institute of Slovenia
