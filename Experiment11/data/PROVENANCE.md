# Provenance of the ClSSCl and SSCl2 photoabsorption cross sections

## Source

Trabelsi, T., and Francisco, J. S. (2026), "Chlorine-sulfur isomers as parents of
ClS2 and SCl2 on Venus: Spectroscopy and photochemistry of ClSSCl, SSCl2, and
(ClS)2", Journal of Chemical Physics 164, 084302. DOI 10.1063/5.0317311.

These are the numerical data underlying Figure 2 of that paper: state-resolved
photoabsorption cross sections computed at the EOM-CCSD/aug-cc-pV(T+d)Z level
using the nuclear ensemble approach, with 1000 Wigner-sampled geometries at
T = 0 K and 10 electronic states.

## How they were obtained

The paper's Supplementary Material section states that the numerical data
underlying Figure 2 are provided as separate files named `csec-clsscl.docx` and
`csec-sscl2.docx`. Those files were not present in the figshare deposit for the
paper: collection 8293312 contains a single article (31266793) holding a single
file, `Supplementary_information.docx`, which contains Figures S1 and S2 as
images only.

The files here were supplied directly by the corresponding authors on request,
in response to an email sent on 26 September 2026 under the paper's stated data
availability policy. They should be cited as the paper's supplementary data,
with acknowledgement to the authors for supplying them.

## File format

Two columns, whitespace separated, one header line:

    # Wavelength (nm) Cross-section (cm^2 molecule^-1)

1001 data rows each, ordered from long to short wavelength. The grid is uniform
in photon energy from 0.1 to 10 eV at a spacing of 0.00990 eV, which corresponds
to 12398.42 nm down to 123.98 nm. Within 300 to 620 nm there are 216 points,
denser than the model's 2 nm spacing in that range, so mapping onto the VULCAN
wavelength grid is a downsampling.

## Validation against the published figure

| Quantity | Published Figure 2 | These files |
| --- | --- | --- |
| SSCl2 global peak | about 3.37e-17 near 262 nm | 3.3738e-17 at 264.2 nm |
| SSCl2 near-UV band | about 9.3e-18 near 345 nm | 9.3255e-18 at 344.0 nm |
| ClSSCl peak | about 2.6e-17 near 240 nm | 2.6279e-17 at 240.3 nm |

## Note relevant to this project

ClSSCl is the stronger absorber overall, but its strength lies at 240 nm, below
the 300 to 620 nm range in which the Venus unknown UV absorber operates. Inside
that range ClSSCl peaks at only 5.9654e-19 at 300.3 nm, while SSCl2 peaks at
9.3255e-18 at 344.0 nm, a factor of 15.6 stronger.

The two isomers are therefore not interchangeable for the absorber question, and
the VULCAN-Venus network carries them as a single lumped S2Cl2 species.

## Checksums

    csec-sscl2.dat   sha256 e17d40e876d5c0f864f7...
    csec-clsscl.dat  sha256 bb3013dbf03238e96430...
