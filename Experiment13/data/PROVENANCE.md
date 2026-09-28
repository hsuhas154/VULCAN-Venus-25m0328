# Provenance of the rhomboclase and acid ferric sulfate extinction coefficients

## Source

Jiang, C.Z., Rimmer, P.B., Lozano, G.G., et al. (2024), "Iron-sulfur chemistry can
explain the ultraviolet absorber in the clouds of Venus", *Science Advances* 10,
eadg8826. DOI 10.1126/sciadv.adg8826.

The two files here hold the mass extinction coefficients plotted in Figures S9
and S10 of that paper's supplementary material:

- `jiang_acid-ferric-sulfate_extracted.dat` : Figure S9, (H3O)Fe(SO4)2
- `jiang_rhomboclase_extracted.dat` : Figure S10, (H5O2)Fe(SO4)2.3H2O

## Why they had to be extracted

The paper's data availability statement says that all data needed to evaluate the
conclusions are present in the paper or the supplementary materials. For these two
quantities the supplement contains only the plots. There is no table, no data file
and no repository deposit. The figures are therefore the data release.

Unlike the ClSSCl and SSCl2 cross sections used in Experiment 11, these were not
requested from the authors. They did not need to be: the figures are vector
graphics, so the curve is stored in the PDF as an exact ordered list of
coordinates, and reading that list recovers the numbers the plotting program was
given. Nothing is read off a picture.

## How they were produced

`../extract_jiang_spectra_25m0328.py`, run against the supplementary PDF. The
script converts each figure page to SVG with `pdftocairo`, which preserves path
geometry exactly, identifies the smoothed fit curve by its stroke colour, applies
each path's affine transform to reach page coordinates, and maps through the axis
box.

## Validation

Three checks, all fatal, all performed by the script on every run.

1. **Frame identification.** The axis box is taken as the largest inner
   rectangle on the page.

2. **Wavelength scale.** The extracted curve's endpoints must land on the stated
   axis limits of 200 and 800 nm. The endpoints are not used in the calibration,
   so recovering them is independent. Achieved: 199.68 and 800.11 nm for acid
   ferric sulfate, 199.68 and 800.12 nm for rhomboclase, both within 0.4 nm out
   of a 600 nm span.

3. **Extinction scale.** The minor tick marks on the vertical axis fix the scale
   a second time, independently of the frame. Five minor intervals make one major
   interval of 100 cm2/g, so the frame height divided by that must equal the whole
   number of major divisions the axis shows. Achieved: 8.000 against 8 for acid
   ferric sulfate, 5.000 against 5 for rhomboclase.

Check 3 exists because an earlier version of this extraction selected an inner
clipping rectangle rather than the axis box and inflated the acid ferric sulfate
peak by 6.4 per cent, from 746.8 to 794.5 cm2/g. Under check 3 that rectangle
yields 7.52 major divisions instead of 8 and the script stops. This was verified
by deliberately forcing the wrong rectangle and confirming the failure.

## File format

Two columns, whitespace separated, three header lines beginning with `#`:

    wavelength_nm  extinction_cm2_per_g

Acid ferric sulfate: 187 points, 199.68 to 800.11 nm, peak 746.8 cm2/g at 328.2 nm.
Rhomboclase: 331 points, 199.68 to 800.12 nm, peak 453.9 cm2/g at 358.9 nm.

## What the numbers are, and one caveat

They are the authors' smoothed fits, not their raw measurements. Both figures plot
scatter alongside the fitted curve and the scatter is visibly a few per cent wide.
The curve extracted here is the fit.

They are also **extinction** coefficients, derived from optical depth measured
through a cuvette holding a suspension, so they contain scattering as well as
absorption. Any use of them as absorption coefficients overstates the absorption,
which makes every shortfall computed from them a lower bound.
