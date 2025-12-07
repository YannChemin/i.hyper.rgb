## DESCRIPTION

*i.hyper.rgb* creates RGB or CMYK composites from hyperspectral imagery
imported as 3D raster maps (`raster_3d`) by
[i.hyper.import](i.hyper.import.html).

The module reads wavelength metadata from hyperspectral 3D raster bands
and automatically selects the bands closest to specified target
wavelengths for each color channel. Users can choose between RGB
(red-green-blue) and CMYK (cyan-magenta-yellow-key) color spaces.

*i.hyper.rgb* is part of the **i.hyper** module family designed for
hyperspectral data import, processing, and analysis in GRASS. It
complements [i.hyper.composite](i.hyper.composite.html) by providing
color space conversions with advanced statistical band selection and
accessibility features.

The module offers multiple statistical methods for band selection when
multiple bands fall within the target wavelength range:

- **mean** -- Average of candidate bands
- **median** -- Median value of candidate bands
- **mode** -- Most frequent value of candidate bands
- **min** -- Minimum value of candidate bands
- **max** -- Maximum value of candidate bands
- **sd1_pos, sd2_pos, sd3_pos** -- Mean plus 1, 2, or 3 standard
  deviations
- **sd1_neg, sd2_neg, sd3_neg** -- Mean minus 1, 2, or 3 standard
  deviations

Color blind safe palette adjustments are available for protanopia
(red-blind), deuteranopia (green-blind), and tritanopia (blue-blind)
vision deficiencies, improving accessibility of hyperspectral
visualizations.

The resulting output consists of individual raster maps for each color
channel (e.g., `output_red`, `output_green`, `output_blue`) organized
into an image group (`output_rgb` or `output_cmyk`) for convenient
display and further processing.

## NOTES

The module expects input data to be a 3D raster map created by
*i.hyper.import* or any 3D raster with wavelength metadata stored in
band-level metadata following the *i.hyper* standard format:
**wavelength**, **FWHM**, **valid**, and **unit**.

When wavelength metadata is not found, the module falls back to using
band indices as wavelength values, which may produce unexpected results.
It is recommended to use properly imported hyperspectral data.

The **-n** flag normalizes output bands to the 0-255 range, suitable for
8-bit display formats. Without normalization, output bands retain their
original reflectance or radiance values.

Default wavelengths are chosen to match typical RGB color perception:

- **Red**: 650 nm
- **Green**: 550 nm
- **Blue**: 450 nm

For CMYK composites:

- **Cyan**: 490 nm
- **Magenta**: 580 nm
- **Yellow**: 570 nm
- **Key (black)**: 800 nm (near-infrared for contrast)

These defaults can be customized using the wavelength options to match
specific sensor characteristics or visualization needs.

Color blind adjustments apply transformation matrices to the output
channels, simulating how the composite would appear to individuals with
different types of color vision deficiency. This feature helps ensure
that visualizations convey information effectively to all viewers.

## EXAMPLES

::: code

    # Create a standard RGB composite from PRISMA data
    i.hyper.rgb input=prisma \
                output=prisma_composite \
                colorspace=rgb \
                statistic=mean

    # Result: prisma_composite_red, prisma_composite_green, prisma_composite_blue
    # Image group: prisma_composite_rgb
:::

::: code

    # Create a normalized RGB composite suitable for display
    i.hyper.rgb input=enmap \
                output=enmap_display \
                colorspace=rgb \
                statistic=mean \
                -n

    # Display the composite
    d.rgb red=enmap_display_red \
          green=enmap_display_green \
          blue=enmap_display_blue
:::

::: code

    # Create a deuteranopia-adjusted RGB composite
    # Useful for ensuring visualization accessibility
    i.hyper.rgb input=tanager \
                output=tanager_colorblind \
                colorspace=rgb \
                statistic=median \
                colorblind=deuteranopia \
                -n
:::

::: code

    # Create CMYK composite with custom wavelengths
    i.hyper.rgb input=prisma \
                output=prisma_cmyk \
                colorspace=cmyk \
                cyan_wavelength=480 \
                magenta_wavelength=590 \
                yellow_wavelength=560 \
                key_wavelength=850 \
                statistic=mean
:::

::: code

    # Use standard deviation for enhanced contrast
    i.hyper.rgb input=enmap \
                output=enmap_enhanced \
                colorspace=rgb \
                statistic=sd2_pos \
                red_wavelength=660 \
                green_wavelength=560 \
                blue_wavelength=470 \
                -n
:::

::: code

    # Create multiple versions for comparison
    # Standard RGB
    i.hyper.rgb input=hyperspectral output=standard colorspace=rgb -n

    # Protanopia-adjusted
    i.hyper.rgb input=hyperspectral output=protanopia colorspace=rgb \
                colorblind=protanopia -n

    # Tritanopia-adjusted
    i.hyper.rgb input=hyperspectral output=tritanopia colorspace=rgb \
                colorblind=tritanopia -n
:::

## SEE ALSO

[i.hyper.import](i.hyper.import.html),
[i.hyper.composite](i.hyper.composite.html),
[i.hyper.preproc](i.hyper.preproc.html),
[i.hyper.explore](i.hyper.explore.html),
[i.hyper.export](i.hyper.export.html),
[d.rgb](https://grass.osgeo.org/grass-stable/manuals/d.rgb.html),
[i.group](https://grass.osgeo.org/grass-stable/manuals/i.group.html),
[r3.support](https://grass.osgeo.org/grass-stable/manuals/r3.support.html),
[r.rescale](https://grass.osgeo.org/grass-stable/manuals/r.rescale.html)

## REFERENCES

- Viénot, F., Brettel, H., & Mollon, J. D. (1999). Digital video
  colourmaps for checking the legibility of displays by dichromats.
  *Color Research & Application*, 24(4), 243-252.
- Machado, G. M., Oliveira, M. M., & Fernandes, L. A. (2009). A
  physiologically-based model for simulation of color vision deficiency.
  *IEEE Transactions on Visualization and Computer Graphics*, 15(6),
  1291-1298.

## AUTHORS

Created for the i.hyper module family

Based on work by Alen Mangafić and Tomaž Žagar, Geodetic Institute of
Slovenia
