# VULCAN-Venus

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3](https://img.shields.io/badge/python-3.x-blue.svg)](https://www.python.org/)

**An open-access chemistry-transport model for the middle and lower atmosphere of Venus, written in Python 3.**

VULCAN-Venus solves the coupled photochemistry and vertical transport of a Venus-like
atmosphere (nitrogen, carbon, hydrogen, oxygen, sulfur, and chlorine chemistry) and
integrates it to a steady state. It is built on the exoplanet chemical kinetics code
[VULCAN](https://github.com/exoclime/VULCAN) by
[Tsai et al. (2017)](https://arxiv.org/abs/1607.00409).

> ### About this fork
> This repository is a personal research fork maintained by **Suhas H. (Roll No. 25m0328)**
> for an M.Tech. Project (MTP) investigating the **Unknown UV Absorber (UUV)** in the
> Venusian cloud region. It preserves the upstream model unchanged in behaviour and adds a
> controlled experiment framework on top of it (see [This fork: MTP experiment
> framework](#this-fork-mtp-experiment-framework)). All scientific credit for the model
> itself belongs to the original author. If you use the model, **cite Dai et al. (2024)** as
> described below.

---

## Citation and credit

The model and its theory are described in:

> Dai, L. et al. (2024). *VULCAN-Venus.* Astronomy & Astrophysics.
> [doi:10.1051/0004-6361/202450552](https://doi.org/10.1051/0004-6361/202450552)

The underlying VULCAN framework is described in
[Tsai et al. (2017)](https://arxiv.org/abs/1607.00409).

- **Original author of VULCAN-Venus:** Longkang Dai (`dailongkang2022@126.com`)
- **This MTP fork:** Suhas H. (25m0328)

---

## Table of contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Repository structure](#repository-structure)
- [Configuration file](#configuration-file)
- [Input files](#input-files)
- [Editing or using a different chemical network](#editing-or-using-a-different-chemical-network)
- [Boundary conditions](#boundary-conditions)
- [Reading output files](#reading-output-files)
- [This fork: MTP experiment framework](#this-fork-mtp-experiment-framework)
- [Modelling notes](#modelling-notes)
- [Remarks and funding](#remarks-and-funding)
- [License](#license)

---

## Requirements

VULCAN-Venus runs on **Python 3**. Two convenient ways to manage the environment:

- [pip](https://pip.pypa.io/en/stable/) : the Python package installer
- [Anaconda / conda](https://docs.continuum.io/) : virtual-environment manager

Python packages:

| Package | Required? | Notes |
| --- | --- | --- |
| numpy | yes | |
| scipy | yes | |
| sympy | yes | Used by `make_chem_funs.py` to build the chemistry |
| matplotlib | yes | |
| PIL / Pillow | optional | Only for interactive (real-time) plotting via the OS image viewer |
| [PyMieScatt](https://github.com/bsumlin/PyMieScatt) | optional | Only needed if cloud formation is enabled. If unused, remove `import PyMieScatt as mie` from `op.py` |

A standard **C++ compiler** (`g++` or `clang`) is required to build the embedded
[FastChem](https://github.com/exoclime/FastChem) equilibrium-chemistry code.

## Installation

If any packages are missing, install the SciPy stack with pip (use `pip` instead of `pip3`
on systems where Python 3 is the default):

```bash
pip3 install --upgrade pip
pip3 install --user numpy scipy matplotlib ipython jupyter pandas sympy nose
pip3 install PyMieScatt        # optional, only for cloud formation
```

Further help: <https://scipy.org/install/> and, for Pillow,
<https://github.com/python-pillow/Pillow>.

### Reproducible environment (this fork)

This fork pins its environment so runs can be reproduced. Recreate it with conda:

```bash
conda env create -f environment/environment.yml
conda activate venus
```

`environment/venus-environment-25m0328.txt` is a frozen snapshot of the exact package
versions used for the MTP runs. A full setup-and-reproduction walkthrough lives in
[`documentation/SETUP_AND_REPRODUCTION.md`](documentation/SETUP_AND_REPRODUCTION.md).

## Quick start

1. **Build FastChem** (equilibrium-chemistry initialiser). From the `fastchem_vulcan/`
   folder:

   ```bash
   cd fastchem_vulcan
   make
   cd ..
   ```

2. **Run the model** from the project root:

   ```bash
   python3 vulcan.py
   ```

   On startup `make_chem_funs.py` regenerates `chem_funs.py` from the active chemical
   network (a few seconds), then the main integration begins.

3. **Read the output.** Results are written to `output/` as a pickled `.vul` file by
   default. `output/reading_data.py` shows how to unpack and inspect it.

4. **Change the setup.** Almost every run is configured by editing `vulcan_cfg.py` alone.
   For example, to switch the vertical mixing (K<sub>zz</sub>) profile, edit:

   ```python
   atm_file = 'atm/atm_Venus_B+Dkzz.txt'
   ```

   Set `use_live_plot = True` / `False` to toggle real-time plotting. Have fun!

---

## Repository structure

The tree below highlights the files you are most likely to touch. Large data
directories (per-species thermodynamic and cross-section files) are summarised rather
than listed in full.

```
VULCAN-Venus-25m0328/
│
│  ── Core model ────────────────────────────────────────────────
├── vulcan.py                     # Top-level driver script
├── vulcan_cfg.py                 # ← main configuration file (edit this per run)
├── vulcan_cfg_README.txt         # Cheatsheet: what every cfg parameter does
├── op.py                         # Numerical core: reaction rates, radiative transfer,
│                                 #   ODE solver, and the UUV opacity term
├── store.py                      # Declaration of all classes and variables
├── build_atm.py                  # Builds atmospheric structure + initial composition
├── chem_funs.py                  # AUTO-GENERATED chemistry: sources/sinks, Jacobian,
│                                 #   equilibrium constants (do not edit by hand)
├── make_chem_funs.py             # Generates chem_funs.py from the chemical network
├── modify_chem.py                # Manual adjustments to complex rate coefficients
├── phy_const.py                  # Physical constants
│
├── atm/                          # Input atmospheres and radiative inputs
│   ├── atm_Venus_Bkzz.txt        #   T-P-Kzz profiles
│   ├── atm_Venus_B+Dkzz.txt
│   ├── ini_atm_Venus.txt
│   ├── alt_venus.txt
│   ├── UV_absorber.txt           #   Unknown UV Absorber input
│   ├── mode1+2.txt               #   Cloud mode-1 / mode-2 particle profile
│   ├── BC_top_Venus.txt          #   Top-boundary flux
│   └── stellar_flux/             #   Incident / solar spectra (+ helper scripts)
│
├── thermo/                       # Chemical kinetics + thermodynamic data
│   ├── Networks/Nominal.txt      #   Default N-C-H-O-S-Cl reaction network
│   ├── NASA9/                    #   ~230 NASA-9 Gibbs-energy files (one per species)
│   ├── photo_cross/              #   Photolysis cross sections per species + rayleigh/
│   ├── all_compose.txt           #   Atom counts and molecular weights
│   ├── gibbs_text.txt            #   Template used by make_chem_funs.py
│   └── make_compose.py
│
├── fastchem_vulcan/              # Embedded FastChem (C++) equilibrium-chemistry code
│   ├── fastchem_src/  input/  model_main/  obj/
│   └── fastchem                  #   compiled executable (after `make`)
│
├── output/                       # Model outputs (.vul) + reading_data.py
│
│  ── Fork additions (MTP / 25m0328) ─────────────────────────────
├── documentation/
│   ├── SETUP_AND_REPRODUCTION.md                 # End-to-end setup + reproduction guide
│   └── VULCAN_Venus_MTP_Progress_..._25m0328.tex # Progress checkpoint write-up
├── environment/
│   ├── environment.yml                           # conda environment spec
│   └── venus-environment-25m0328.txt             # frozen package list
│
├── op-25m0328-pre-UUV-experiment.py              # Pristine op.py snapshot before UUV work
├── vulcan_cfg-25m0328-pre-UUV-experiment.py      # Pristine cfg snapshot before UUV work
├── make_chem_funs-backup.py                      # Backups of the chemistry generator
├── make_chem_funs-25m0328.py.bak
├── baseline-dictionary-inventory-20260910-071333.txt   # Baseline audit inventory
│
├── diagnose.py                                   # General model diagnostics
├── stage7C1_SO2_rates.py                         # SO2 reaction-rate analysis
├── stage7C2_SO2_condensation_decomposition.py    # SO2 condensation / decomposition
├── stage7C3_final.py                             # Stage-7 consolidation
├── analyze_Experiment3_sensitivity_25m0328.py    # UUV amplitude-sensitivity analysis
│
└── README.md
```

> Typically **`vulcan_cfg.py` is the only file you edit for a given run.** If you want to
> look inside or modify the code itself, `store.py` is where nearly all classes and
> variables are declared.

Purpose of the core modules at a glance:

- `build_atm.py` : constructs the atmospheric structure from the input and sets up the initial composition.
- `chem_funs.py` : chemical source/sink functions, the Jacobian matrix, and equilibrium constants (regenerated automatically).
- `op.py` : all numerical operations (reaction rates, radiative transfer, ODE solvers).
- `make_chem_funs.py` : runs first to produce `chem_funs.py` from the assigned network.
- `modify_chem.py` : extra handling for complex rate-coefficient expressions.
- `phy_const.py` : physical constants.
- `store.py` : storage of all variables and classes.
- `vulcan.py` : the top-level main script.
- `vulcan_cfg.py` : the run configuration.

---

## Configuration file

**All settings and parameters (atmospheric parameters, elemental abundances, solver
options, and so on) are prescribed in `vulcan_cfg.py`.** This is typically the only file
you edit for a specific run. A cheatsheet explaining every parameter is provided in
`vulcan_cfg_README.txt`.

## Input files

The key inputs are the chemical network, the atmospheric T-P-K<sub>zz</sub> profile, and
the incident actinic flux.

**Chemical network.** `thermo/Networks/Nominal.txt` is the default reaction network
(nitrogen, carbon, hydrogen, oxygen, sulfur, and chlorine species). Rate coefficients are
written as `A`, `B`, `C` in the Arrhenius form:

```
k = A * T^B * exp(-C / T)
```

**Temperature-pressure-K<sub>zz</sub> profile.** Required when `Kzz_prof = 'file'` in
`vulcan_cfg.py`; placed in `atm/` by default. The first line is a comment for units, and
the second line must name the columns: **`Pressure  Temp  Kzz`** (K<sub>zz</sub> optional).
So the file has two columns without K<sub>zz</sub> and three columns with it.

**Stellar / actinic flux.** Stored in `atm/stellar_flux/`, with column 1 the wavelength in
nm and column 2 the flux in erg cm<sup>-2</sup> s<sup>-1</sup> nm<sup>-1</sup>. If
`use_sflux_top = True`, the surface stellar file is not read; instead the flux in
`sflux_top_path` is used as the incident flux at the top boundary.

**Thermodynamics and cross sections.** Stored in `thermo/NASA9/` and `thermo/photo_cross/`
respectively. *Change at your own risk.* A few species that lack experimental
thermodynamic data (e.g. `ClCO3`) are filled with dummy values to keep the code stable. In
that case the reverse reactions are supplied directly, and the automatic calculation of
reverse reactions from the forward ones is switched off via
`remove_list = [i*2 for i in range(1,1000)]` in `vulcan_cfg.py`, so the thermodynamic data
do **not** influence the nominal results. The reverse-reaction machinery is retained for
future modifications.

**Fixed species fluxes.** If constant fluxes for certain species are used, those files also
live in `atm/`, in the format: species, flux (cm<sup>-2</sup> s<sup>-1</sup>), and
deposition velocity (cm s<sup>-1</sup>).

---

## Editing or using a different chemical network

VULCAN's chemical network is **not hard-coded**. `make_chem_funs.py` generates every
required function from the input network into `chem_funs.py`. You can edit the default
network (remove or add reactions, change rate constants) or supply a different one, as long
as it uses the same format: reactions are written in square brackets, e.g. `[ A + B -> C + D ]`.

By default `make_chem_funs.py` runs before the main code to (re)build `chem_funs.py`. This
step (a few seconds) can be skipped with the `-n` flag:

```bash
python3 vulcan.py -n
```

**Do not skip this step after changing the network.** For complex rate-coefficient
expressions, `modify_chem.py` is the place for further adjustment.

Changing the network is not foolproof: unrealistic values can cause numerical issues.

To switch on the calculation of reverse reactions from the forward ones, set
`remove_list = []` and remove all reverse reactions from the network file. Be careful with
species that lack thermodynamic data or contain sulfur. Ensure every species is present in
`thermo/NASA9/`; if not, add it manually by consulting `nasa9_2002_E.txt` or
`new_nasa9.txt` (both in `NASA9/`) and saving the coefficients in a text file named after
the species (e.g. `CO2.txt`). The NASA-9 polynomial format is:

```
a1 a2 a3 a4 a5
a6 a7 0. a8 a9
```

Here `a7` and `a8` are separated by `0.`. The first two rows cover low temperature
(200-1000 K) and the last two cover high temperature (1000-6000 K).

Reaction IDs are irrelevant: they are generated automatically (and written back into the
network file) by `make_chem_funs`. Three-body and dissociation reactions must be listed
separately after the comment line, as in the default network. After changing the network,
the readable reaction and species lists in `chem_funs.py` are updated whenever you run
`python3 vulcan.py` (without `-n`).

## Boundary conditions

If both `use_topflux` and `use_botflux` are `False`, the default zero-flux boundary
condition applies (nothing enters or leaves).

- `use_topflux = True` reads `top_BC_flux_file` for the incoming/outgoing flux at the top boundary.
- `use_botflux = True` reads `bot_BC_flux_file` for the surface pressure and sinks at the bottom boundary.
- `use_fix_sp_bot` fixes surface mole fractions, e.g. `use_fix_sp_bot = {'CO2': 0.96}` sets the surface CO<sub>2</sub> mixing ratio to 0.96.

## Reading output files

`output/reading_data.py` is a good example of how to access the outputs. Unpack the binary
`.vul` file with `pickle.load`; the main variables live in three classes:
`data['variable']`, `data['atm']`, and `data['parameter']`. All variable names and the
class structure are documented in `store.py`.

---

## This fork: MTP experiment framework

This section documents what the fork adds on top of the upstream model. The scientific
target is the **Unknown UV Absorber (UUV)** in the Venusian cloud region.

### Design principle

The upstream model is treated as an **immutable reference**. Before touching the UUV,
pristine snapshots of the code were preserved (`op-25m0328-pre-UUV-experiment.py`,
`vulcan_cfg-25m0328-pre-UUV-experiment.py`), a baseline was reproduced and audited, and
only then were controlled experiments layered on. Every experiment writes a uniquely named
output plus a configuration snapshot, so no result is ever overwritten.

### The UUV amplitude parameter

The experiments scale the explicit UUV optical-depth term by a dimensionless multiplier
`A_UUV`:

```
tau_UUV*(z, lambda) = A_UUV * tau_UUV(z, lambda)
```

with the mode-1 / mode-2 Mie terms held disabled so the UUV contribution is isolated.

### Experiment matrix

| Experiment | `A_UUV` | Role | Output file |
| --- | --- | --- | --- |
| Experiment 1 | 0.0 | Implementation control | `output/Experiment1_A0_UUV-control-25m0328.vul` |
| Experiment 2 | 1.0 | Nominal supplied representation | `output/Experiment2_A1_UUV-nominal-25m0328.vul` |
| Experiment 3A | 0.5 | Half nominal strength (sensitivity) | `output/Experiment3_A0p5_UUV-half-25m0328.vul` |
| Experiment 3B | 2.0 | Double nominal strength (sensitivity) | `output/Experiment3_A2_UUV-double-25m0328.vul` |

Each run has a matching `cfg_*.txt` snapshot in `output/`.

### Baselines and audit trail

- `Nominal_Bkzz_SO2.vul` : the reference steady state used as the initial composition.
- `Nominal_Bkzz_SO2-25m0328-reference.vul`, `-baseline.vul`, `-laptop-baseline.vul` :
  preserved reference / reproduced baselines for the reproducibility audit.
- `baseline-dictionary-inventory-20260910-071333.txt` : inventory of the baseline output
  dictionary, used to verify the reproduction.
- `stage4_tau_comparison-*.txt`, `stage5_aflux_comparison-*.txt`,
  `stage6_photolysis_audit-*.txt`, `stage7A_absorber_decomposition-*.txt`,
  `stage7B_SO2_atmospheric_trace-*.txt` : diagnostic logs from successive audit stages
  (optical depth, actinic flux, photolysis, absorber decomposition, SO<sub>2</sub> trace).

### Analysis and diagnostic scripts

| Script | Purpose |
| --- | --- |
| `diagnose.py` | General model diagnostics |
| `stage7C1_SO2_rates.py` | SO<sub>2</sub> reaction-rate breakdown |
| `stage7C2_SO2_condensation_decomposition.py` | SO<sub>2</sub> condensation / decomposition analysis |
| `stage7C3_final.py` | Stage-7 consolidation |
| `analyze_Experiment3_sensitivity_25m0328.py` | UUV amplitude-sensitivity analysis across the Experiment-3 runs |

### Reproducing an experiment

1. Recreate the environment (see [Reproducible environment](#reproducible-environment-this-fork)).
2. Set `A_UUV` and the output name in `vulcan_cfg.py` to match the target experiment.
3. Run `python3 vulcan.py` and let it reach steady state.
4. Compare against the control with the relevant `stage*` / `analyze_*` script.

See [`documentation/SETUP_AND_REPRODUCTION.md`](documentation/SETUP_AND_REPRODUCTION.md)
for the full procedure.

---

## Modelling notes

The dissolution and release rate coefficients of SO<sub>2</sub> were adopted from Rimmer et
al. (2021) to capture the SO<sub>2</sub> depletion in the cloud region. These coefficients
have not yet been verified experimentally. The aqueous cloud chemistry that removes liquid
SO<sub>2</sub> in Rimmer et al. (2021) is **not** included here, which may overestimate
liquid SO<sub>2</sub>. Because liquid SO<sub>2</sub> does not react with any other species,
this does not affect the results, and the issue will be resolved once cloud formation is
coupled into VULCAN-Venus.

## Remarks and funding

The original project received financial support from the National Natural Science
Foundation of China and the Natural Science Foundation of Hunan Province.

VULCAN-Venus is publicly accessible. To use it, please cite
[Dai et al. (2024, doi:10.1051/0004-6361/202450552)](https://doi.org/10.1051/0004-6361/202450552).

Questions or suggestions about the model itself: contact Longkang Dai at
`dailongkang2022@126.com`. Questions about this MTP fork specifically: Suhas H. (25m0328).

## License

Released under the **GNU General Public License v3.0** (GPLv3). See
<https://www.gnu.org/licenses/gpl-3.0> for the full text.
