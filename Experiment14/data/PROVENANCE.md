# Provenance of the ferric chloride absorptivity curves

## Source

Egan, J. V., James, A. D., and Plane, J. M. C. (2025), "Laboratory Measurements of
Ferric Chloride (FeCl3) under Venusian Conditions", *ACS Earth and Space Chemistry*,
doi:10.1021/acsearthspacechem.5c00132. Published open access under CC-BY 4.0.

The five files here hold the molar absorptivity curves plotted in Figure 5a of that
paper, the ferric chloride ion spectra at five HCl concentrations:

| file | curve | HCl |
|---|---|---|
| `fecl3_ClA_eps.dat` | Cl A | 12.2 M |
| `fecl3_ClB_eps.dat` | Cl B | 8.96 M |
| `fecl3_ClC_eps.dat` | Cl C | 6.09 M |
| `fecl3_ClD_eps.dat` | Cl D | 4.15 M |
| `fecl3_ClE_eps.dat` | Cl E | 1.54 M |

## Why they had to be extracted

The values appear in the paper only as a figure. There is no table and no data
deposit. The figure is therefore the data release.

They are needed because the published rejection of ferric chloride as the Venus
absorber (Spacek et al. 2026) rests on an **aqueous** absorbance: 1 g/L of aqueous
FeCl3 giving a(375) of about 1 cm-1, linearly extrapolated to 1.3 kg/L. Venus cloud
droplets are roughly 75 wt% sulfuric acid. Egan et al. state that the spectrum in
concentrated sulfuric acid "does not appear to have been measured previously" and
that "the choice of solvent has a significant effect on the measured spectrum".
Figure 5a is that measurement.

## There is no single ferric chloride spectrum

In sulfuric acid the iron speciates into chloro complexes whose populations are set
by chloride activity. The paper describes the progression with increasing chloride:
a single peak above 300 nm (FeCl2+), broader and stronger (FeCl2+), two separate
peaks (FeCl3), then a pronounced minimum and higher overall absorbance (FeCl4-).
All five curves are therefore extracted and every downstream result is reported as a
range rather than a single number.

## How they were produced

`../extract_egan_fig5a_25m0328.py`, run against the published PDF. Page 6 is
rendered at 600 dpi with `pdftoppm` and each curve is followed by its saturated
stroke colour. The shaded uncertainty bands are lighter blends of the same hue and
are excluded by requiring saturation. Nothing is read off by eye.

Unlike the Jiang figures of Experiment 13, the vector paths here are fragmented into
more than a hundred pieces per curve by the uncertainty bands, so colour tracing of
the SVG is unreliable and raster extraction at high resolution is used instead.

## Validation

Three checks, all fatal, all performed on every run.

1. **Panel identification.** The page must contain four panel verticals, and the
   horizontal edges must span the panel width. A page rule elsewhere on the page can
   otherwise be mistaken for the panel top.

2. **Tick uniformity.** Spacing on both axes must be uniform to better than 1 per
   cent. Ticks are anchored on the axis end and accepted only if they sit on the
   regular grid, which rejects a spurious tick dragged in by a misidentified edge.
   Achieved: 4.6829 px/nm horizontally and 260.57 px per 1000 units vertically.

3. **Baseline return.** Every curve must come back to zero at the long-wavelength
   end, within a tenth of a tick. The calibration knows nothing about this, so
   recovering it is independent evidence that the vertical scale is right.
   Achieved: tails of 14 to 36 against peaks of 4931 to 7413.

### Sabotage testing

The checks were verified by deliberately breaking the extraction three ways:

- pointing it at the wrong page: caught by check 1 (two verticals, not four);
- disabling the grid filter on ticks: caught by check 2 (17.42 per cent non-uniform);
- displacing the vertical zero by 40 pixels: **initially passed**, because check 3
  allowed 3 per cent of peak while the real tails are near 20 units. The tolerance
  was tightened to a tenth of a tick and the check now fires on that shift. The
  sabotage test found a real gap in the guard rather than confirming one.

## File format

Two columns, whitespace separated, three header lines beginning with `#`:

    wavelength_nm   eps_Fe_L_per_mol_per_cm

The absorptivity is per mole of iron, as the figure's axis states.

## Caveat

These are the paper's derived component spectra for ferric chloride ions, not raw
sample measurements. Which curve applies to a Venus droplet depends on the chloride
activity there, which this extraction does not determine. That is why the result is
a range.
