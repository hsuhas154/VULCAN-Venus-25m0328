# VULCAN-Venus: Setup and Reproduction Guide

## Project

Project: Venus photochemistry / Unknown UV Absorber investigation

Repository:
`https://github.com/hsuhas154/VULCAN-Venus-25m0328.git`

Upstream VULCAN-Venus repository:
`https://github.com/dai-lk/VULCAN-Venus.git`

The local repository keeps the upstream repository as the `origin` remote and the personal research repository as the `research` remote.

## Working directory

The intended project root is:

`~/projects-25m0328/venus-photochemistry/VULCAN-Venus`

The parent project directory is:

`~/projects-25m0328/venus-photochemistry`

The embedded FastChem component is located at:

`fastchem_vulcan/`

## Operating system and build environment

The laboratory setup used:

- WSL2
- Ubuntu 24.04.4 LTS
- GCC/G++/GFortran 13.3.0
- GNU Make 4.3
- Miniforge
- Conda environment named `venus`

CMake was not required for this repository because FastChem uses the supplied Makefile system.

## Python environment

The primary environment is `venus`.

Python version: `3.14.7`

The repository contains two environment records:

- `environment/venus-environment-25m0328.txt`: the package snapshot generated with `python -m pip freeze`.
- `environment/environment.yml`: the Conda environment specification without a machine-specific prefix.

The environment should be recreated from `environment.yml` on a new machine, followed by the compatibility steps below.

## Creating the Conda environment

From the repository root:

```bash
conda env create -f environment/environment.yml
conda activate venus
```

If an environment named `venus` already exists, update it from the YAML rather than creating a duplicate environment:

```bash
conda env update -n venus -f environment/environment.yml --prune
conda activate venus
```

## Important compatibility correction: `make_chem_funs.py`

The supplied VULCAN-Venus code predates the current Python/NumPy stack.

In `make_chem_funs.py`, a value was being unconditionally processed with:

```python
sp.decode("utf-8")
```

Under the current environment, some values were already NumPy string objects rather than byte strings.

The compatibility correction used in this project was:

```python
compo_row = [sp.decode("utf-8") if isinstance(sp, bytes) else sp for sp in compo_row]
```

The original file was preserved as:

`make_chem_funs-25m0328.py.bak`

After the correction, running:

```bash
python make_chem_funs.py
```

completed successfully and reported:

```text
Including condensation reactions.
Elements conserved in the network.
No duplicates in the network.
```

This was a compatibility correction only. No reaction or scientific parameter was intentionally changed.

## Important compatibility correction: PyMieScatt

The installed PyMieScatt version is `1.8.1.1`.

The supplied package attempted to import:

```python
from scipy.integrate import trapz
```

The installed SciPy version provides the corresponding operation as `trapezoid`.

The affected PyMieScatt source files were backed up before modification:

- `Mie.py.25m0328.bak`
- `Inverse.py.25m0328.bak`

The imports were changed to use:

```python
from scipy.integrate import trapezoid as trapz
```

After this correction:

```bash
python -c "import PyMieScatt"
```

succeeded, and a basic Mie calculation also returned numerical output.

This was a compatibility correction only and did not change the scientific inputs.

## FastChem build

FastChem is embedded in `fastchem_vulcan/`. The repository uses the supplied Makefile system.

Build with:

```bash
cd fastchem_vulcan
make
```

A successful build ends with the repository's completion message:

```text
everything is done and fine. enjoy your day!
```

The generated executable is:

`fastchem_vulcan/fastchem`

Compiler-generated object files and dependency files are intentionally ignored by Git and can be regenerated with `make`.

## Chemistry generation

After the compatibility correction to `make_chem_funs.py`, run:

```bash
python make_chem_funs.py
```

This generates the repository's active `chem_funs.py`.

The generated chemistry file is part of the research checkpoint because it is the chemistry code actually used by the successful baseline calculation.

## Nominal VULCAN-Venus configuration

The nominal calculation uses:

- `thermo/Networks/Nominal.txt`
- `atm/atm_Venus_Bkzz.txt`

Important configuration properties include:

- 57 atmospheric layers
- photochemistry enabled
- eddy diffusion enabled
- molecular diffusion disabled
- condensation enabled
- supplied Venus top actinic flux
- Rayleigh scattering species CO2, N2, H2, and O2
- temperature-dependent absorption species H2O2, CO2, and H2O
- photochemical wavelength grid spanning approximately 96 to 698 nm
- 0.1 nm bins below the transition and 2 nm bins above it

The UUV and particle input files are:

- `atm/UV_absorber.txt`
- `atm/mode1+2.txt`

## Preserving the supplied reference output

The supplied output originally located at:

`output/Nominal_Bkzz_SO2.vul`

was preserved before any fresh integration. The preserved copy is:

`output/Nominal_Bkzz_SO2-25m0328-reference.vul`

The original and preserved reference were verified byte-for-byte identical.

## Baseline reproduction

The nominal model was run with:

```bash
python3 vulcan.py
```

The fresh baseline output was saved as:

`output/Nominal_Bkzz_SO2-25m0328-baseline.vul`

The configuration used for that calculation was archived as:

`output/cfg_Nominal_Bkzz_SO2-25m0328-baseline.txt`

The baseline reported successful convergence after 296 integration steps and a model time of approximately 5.33e8 s.

## Baseline diagnostic sequence

The reference and baseline were compared through the following stages:

1. Saved-state and dictionary inventory
2. Photolysis-rate structure and numerical comparison
3. Volume mixing ratios
4. Optical depth
5. Actinic flux
6. Photolysis consistency
7. Sulfur and absorber tracing:
   - 7A. Absorber decomposition
   - 7B. SO2 atmospheric-field trace
   - 7C. SO2 and H2SO4 chemical tracing

The corresponding diagnostic reports and scripts are retained in the repository.

The key reproducibility finding was that the atmospheric grid, pressure, temperature, altitude, layer thickness, Kzz, wavelength grid, and photochemical cross sections are reproduced consistently.

The principal localized discrepancy is instead an atmospheric sulfur-state difference around levels 27 to 31. Gas-phase H2SO4 and SO3 differ strongly between the reference and baseline, H2SO4 liquid is identical, and the resulting SO2 difference produces the localized optical-depth discrepancy around 200 to 220 nm.

## UUV and aerosol implementation audit

The supplied `op.py` reads:

- `atm/UV_absorber.txt`
- `atm/mode1+2.txt`

and defines the corresponding UUV and aerosol quantities.

However, the Longkang Dai aerosol/UUV optical-depth addition inside `compute_tau()` is currently enclosed in a triple-quoted block and is therefore inactive in the supplied source state.

The active saved optical depth was independently reconstructed from gas absorption and Rayleigh scattering and matched the saved field to very small numerical residuals.

The explicit UUV optical-depth contribution was reconstructed independently from the supplied UUV table and mode-1 particle profile. The reconstructed UUV field was exactly identical between the preserved reference and reproduced baseline.

The mode-1 and mode-2 Mie contributions were also reconstructed independently. They are much smaller than the localized reference-baseline optical-depth discrepancy and do not explain it.

These findings establish that the reproducibility discrepancy should not be interpreted as evidence for the Unknown UV Absorber.

## Pre-UUV experiment checkpoint

Before modifying `op.py` for the first controlled UUV experiment, the source file was copied to:

`op-25m0328-pre-UUV-experiment.py`

The two files were verified to have identical SHA-256 hashes before any experimental modification.

The first planned UUV experiment uses a dimensionless amplitude parameter `A_UUV`, with conceptual definition:

```text
tau_UUV* = A_UUV * tau_UUV
```

The intended first comparison is:

- `A_UUV = 0` as the implementation control
- `A_UUV = 1` as the supplied nominal UUV representation

Only the explicit UUV optical-depth term is to be activated for the first experiment. The mode-1 and mode-2 Mie aerosol terms are not to be activated simultaneously.

## Current scientific status

The project has completed the reproducibility and baseline-audit phase.

The baseline has been preserved, reproduced, and diagnosed.

The UUV experiment has not yet been executed.

The next scientific task is to activate only the explicit UUV optical-depth pathway, introduce the controlled `A_UUV` parameter, verify that `A_UUV = 0` reproduces the existing baseline, and then perform the first `A_UUV = 1` experiment.

## Portability principle

Do not copy the laboratory Conda environment directory itself to another machine.

Instead:

1. Clone the private research repository.
2. Recreate the `venus` environment from `environment/environment.yml`.
3. Apply or verify the documented compatibility corrections.
4. Build FastChem with `make`.
5. Generate `chem_funs.py`.
6. Verify the nominal model execution.
7. Continue from the latest Git checkpoint.

Machine-specific credentials such as SSH private keys must never be committed to the repository.
