# TRIP YAML Configuration File Reference

TRIP is configured through a YAML file that is read at startup. This document
describes every option recognized by the parser (`loadConfig()` in
`src/RT_utility.cpp`) together with its type, default value, and meaning.

## How to use the configuration file

The path of the configuration file is passed on the command line with the
`--yaml_config` option. If the option is omitted, the code looks for a file
named `config.yml` in the current working directory:

```bash
# Use the default ./config.yml
mpirun -n 4 ./main

# Use an explicit configuration file
mpirun -n 4 ./main --yaml_config my_config.yml
```

General rules:

- Keys marked **required** must be present, otherwise the code aborts with
  `Missing required key: <name>`.
- All other keys are optional: when a key is absent, the default listed below
  is used.
- Paths may be absolute or relative to the directory from which the executable
  is launched (see `bin/config.yml` and `tests/*.yml` for working examples).
- Unknown keys are silently ignored, so double-check spelling.

---

## Input

| Key | Type | Required / Default | Description |
|---|---|---|---|
| `input_directory` | string (path) | **required** | Directory containing the atmospheric model input files. |
| `input_file` | string (path) | **required** | 3D atmospheric model file, relative to `input_directory`. Either a PORTA PMD file (`.pmd`) or an HDF5 input file (auto-detected). Not used in FAL-C mode (see below). |
| `frequency_file` | string (path) | **required** | Frequency grid file, relative to `input_directory` (e.g. `frequency/96F/frequency.dat`). |
| `atom_file` | string (path) | none | Path to a separate YAML file with the atomic model (see [Atom file](#atom-file-atom_file)). If omitted, the built-in default two-level CaI atom is used. |
| `input_cul` | string (path) | `""` | Optional PORTA auxiliary file with collisional rates (`.cul`), relative to `input_directory`. |
| `input_qel` | string (path) | `""` | Optional PORTA auxiliary file with elastic collision rates (`.qel`). |
| `input_llp` | string (path) | `""` | Optional PORTA auxiliary file with lower-level population data (`.llp`). |
| `input_back` | string (path) | `""` | Optional PORTA auxiliary file with background/continuum data (`.back`). |

When `input_cul`, `input_qel`, and `input_llp` are all provided, the code runs
in the "PMD + CUL + QEL + LLP + BACK" input mode; otherwise only the PMD/HDF5
file is used.

### FAL-C mode (1D problems)

FAL-C is a separate input format used **only for 1D problems**: it is not the
same kind of input as the PMD/HDF5 3D atmospheric models. FAL-C mode is
selected automatically when the `input_directory` path contains the string
`FAL-C` (e.g. `../input/FAL-C/`). In this mode the 1D FAL-C atmosphere is
used to build the problem, and the horizontal grid is defined by the `N_x`,
`N_y`, and `L` options (see [Numerics](#numerics)).

## Output

| Key | Type | Default | Description |
|---|---|---|---|
| `output` | bool | `false` | Master switch: enable writing of results to disk. |
| `output_directory` | string (path) | `""` | Directory where the results directory is created. |
| `output_overwrite_prevention` | bool | `false` | If `true`, the code stops instead of overwriting an existing results directory. |
| `write_whole_3D_field_hdf5` | bool | `false` | Also write the whole 3D radiation field to an HDF5 file (large output). |
| `write_text_output` | bool | `false` | Also write results as plain-text files. |
| `reference_sol_directory` | string (path) | `""` | Directory with a reference solution, used by the tests to validate results. |
| `verbose` | bool | `true` | Print extra progress/diagnostic information. |

## Physics

| Key | Type | Default | Description |
|---|---|---|---|
| `emissivity_model` | string | **required** | Scattering/emissivity model used in the transfer equation. See [allowed values](#emissivity_model-values). |
| `preconditioner_emissivity_model` | string | `CRD_limit` | Emissivity model used inside the preconditioner. Allowed values: `NONE`, `CRD_limit`, `PRD_AA`, `PRD_AA_MAPV`, `PRD_AA_GB`, `PRD_AA_MAPV_GB`, `ZERO`. |
| `use_B` | bool | `true` | Include the magnetic field (Hanle/Zeeman effects) in the calculation. |
| `use_Vb` | bool | `true` | Include the bulk velocity field (Doppler shifts). |
| `enable_continuum` | bool | `true` | Include the continuum contribution. |
| `enable_sigma_c` | bool | `true` | Include the continuum scattering term sigma_c in the emissivity. Can be disabled for comparisons with PORTA, which does not include it. |
| `use_D2_from_input` | bool | `true` | Take the collisional depolarization rate D2 from the input file. If `false`, D2 is computed from the empirical `a` and `b` coefficients as `D2 = a * nH * (T/5000K)^b` (Derouich 2019). |
| `set_uniform_B` | bool | `false` | Override the model's magnetic field with a uniform field given by `B_field`. |
| `B_field` | sequence of 3 numbers | `[0.0, 0.0, 0.0]` | Uniform magnetic field as `[magnitude (Gauss), theta (rad), chi (rad)]`. Only used when `set_uniform_B: true`. Must contain exactly three numbers. |
| `set_uniform_Vb` | bool | `false` | Override the model's bulk velocity with a uniform velocity given by `Vb_field`. |
| `Vb_field` | sequence of 3 numbers | `[0.0, 0.0, 0.0]` | Uniform bulk velocity as `[Vx, Vy, Vz]` in cm/s. Only used when `set_uniform_Vb: true`. Must contain exactly three numbers. |

### `emissivity_model` values

| Value | Meaning |
|---|---|
| `NONE` | No emissivity model. |
| `CRD_limit` | Complete redistribution (CRD) limit. |
| `CRD_limit_VHP` | CRD limit, Voigt-profile-based variant. |
| `PRD` | Partial redistribution (full treatment). |
| `PRD_NORMAL` | PRD, standard evaluation. |
| `PRD_MEDIUM` | PRD, intermediate-accuracy/speed evaluation. |
| `PRD_FAST` | PRD, fast evaluation. |
| `PRD_AA` | PRD with angle-averaged (AA) redistribution. |
| `PRD_AA_MAPV` | PRD angle-averaged with MAPV approximation. |
| `PRD_AA_GB` | PRD angle-averaged, GB variant. |
| `PRD_TWOTERM` | PRD for a two-term atom. |
| `PRD_AA_TWOTERM` | PRD angle-averaged, two-term atom. |
| `PRD_AA_TWOTERM_MAPV` | PRD angle-averaged with MAPV, two-term atom. |
| `ZERO` | Zero emissivity (debug/testing). |

## Numerics

| Key | Type | Default | Description |
|---|---|---|---|
| `use_1_5D_approx` | bool | `false` | Use the 1.5D approximation (each column treated independently) instead of a full 3D solve. |
| `formal_solver` | string | `BESSER` | Formal solver for the transfer equation. Allowed values: `implicit_Euler`, `trapezoidal`, `CrankNicolson`, `DELO_linear`, `BESSER`. Any other value aborts the run. |
| `N_theta` | int | `8` | Number of inclination (theta) nodes of the angular quadrature. |
| `N_chi` | int | `16` | Number of azimuth (chi) nodes of the angular quadrature. |
| `N_x` | int | `1` | Number of horizontal grid points in x (used only for FAL-C input). |
| `N_y` | int | `1` | Number of horizontal grid points in y (used only for FAL-C input). |
| `L` | double | `400.0` | Horizontal domain size (used only for FAL-C input). |
| `use_prec` | bool | `true` | Use the preconditioner in the iterative solver. **Note:** automatically forced to `false` when `emissivity_model` is `CRD_limit`, `CRD_limit_VHP`, or `ZERO`. |

## Solver section (`solver`)

Nested map controlling the outer PETSc Krylov (KSP) solver:

```yaml
solver:
  ksp_solver_type: KSPFGMRES
  ksp_rtol: 1e-8
  ksp_max_it: 1000
  ksp_use_J_KQ: false
  gmres_restart: 30
```

| Key | Type | Default | Description |
|---|---|---|---|
| `ksp_solver_type` | string | `KSPFGMRES` | PETSc solver type. Allowed values: `KSPGMRES`, `KSPFGMRES`, `KSPPGMRES`, `KSPBCGS`, `KSPPREONLY`, `KSPRICHARDSON`. |
| `ksp_rtol` | double | `1e-5` | Relative convergence tolerance (use ~`1e-6` for single precision, ~`1e-8` for double). |
| `ksp_max_it` | int | `1000` | Maximum number of iterations. |
| `ksp_use_J_KQ` | bool | `false` | Solve in terms of the radiation-field tensors J^K_Q instead of the full radiation field. |
| `gmres_restart` | int | `30` | GMRES restart parameter (only relevant for GMRES-type solvers). |

## Preconditioner section (`prec`)

Nested map controlling the inner (preconditioner) solve. Only relevant when
`use_prec: true`:

```yaml
prec:
  pc_solver_type: KSPGMRES
  pc_rtol: 1e-4
  pc_max_it: 1000
  pc_use_J_KQ: true
  pc_formal_solver_approx: true
  verbose: false
```

| Key | Type | Default | Description |
|---|---|---|---|
| `pc_solver_type` | string | `KSPGMRES` | PETSc solver type for the preconditioner (same allowed values as `ksp_solver_type`). |
| `pc_rtol` | double | `1e-5` | Relative convergence tolerance of the inner solve. |
| `pc_max_it` | int | `1000` | Maximum number of inner iterations. |
| `pc_use_J_KQ` | bool | `false` | Use the J^K_Q formulation inside the preconditioner. |
| `pc_formal_solver_approx` | bool | `false` | Use an approximate (cheaper) formal solver inside the preconditioner. |
| `verbose` | bool | `false` | Print convergence information of the inner solve. |

## Arbitrary beams section (`arbitrary_beams`)

Optional list of extra ray directions along which the emergent radiation is
computed, in addition to the angular quadrature:

```yaml
arbitrary_beams:
  - mu: 1.0
    chi: 0.0
  - mu: 0.3
    chi: 1.57
```

| Key | Type | Default | Description |
|---|---|---|---|
| `mu` | double | `1.0` | Cosine of the heliocentric (inclination) angle of the beam. |
| `chi` | double | `0.0` | Azimuth of the beam in radians. |

## Atom file (`atom_file`)

The atomic model is given in a **separate YAML file** referenced by the
`atom_file` key of the main configuration. If `atom_file` is not given, the
built-in default two-level CaI model is used (the defaults below). Example:

```yaml
# atom_CaI.yml
mass: 40.078        # atomic mass (amu)
Aul: 2.18e+08       # Einstein A coefficient (s^-1)
S2: 0               # 2*S (twice the spin quantum number)
Ll2: 0              # 2*L of the lower level/term
Lu2: 2              # 2*L of the upper level/term
El: 0.0             # lower level energy(ies) (cm^-1)
Eu: 23652.304       # upper level energy(ies) (cm^-1)
gl: 0.0             # lower level Landé factor(s)
gu: 1.0             # upper level Landé factor(s)
```

| Key | Type | Required / Default | Description |
|---|---|---|---|
| `mass` | double | **required** | Atomic mass in atomic mass units. |
| `Aul` | double | **required** | Einstein coefficient for spontaneous emission, in s^-1. |
| `S2` | int | **required** | Twice the spin quantum number (2S). |
| `Ll2` | int | **required** | Twice the orbital quantum number of the lower level/term (2L_l). |
| `Lu2` | int | **required** | Twice the orbital quantum number of the upper level/term (2L_u). |
| `El` | scalar or sequence | **required** | Energy of the lower level(s) in cm^-1. A single number for a two-level atom, a list for a two-term atom. |
| `Eu` | scalar or sequence | **required** | Energy of the upper level(s) in cm^-1. Scalar or list, as above. |
| `gl` | scalar or sequence | `0.0` | Landé factor(s) of the lower level(s). |
| `gu` | scalar or sequence | `1.0` | Landé factor(s) of the upper level(s). |

`El`, `Eu`, `gl`, and `gu` accept either a single value or a YAML list; lists
are used for two-term atoms (see also the `PRD_*TWOTERM*` emissivity models).

---

## Complete example

```yaml
# I/O
input_directory: ../input/PORTA/
input_file:      cai_0Bx_0By_0Bz_1Vx_1Vy_1Vz_GT4_5x5x133_it100.pmd
frequency_file:  frequency/96F/frequency.dat
atom_file:       ../input/atom_CaI.yml

# Output
output: true
output_directory: ../../output
output_overwrite_prevention: false
write_text_output: true
verbose: true

# Physics
emissivity_model: PRD
preconditioner_emissivity_model: CRD_limit
use_B: true
enable_continuum: true
set_uniform_B: true
B_field: [12.0, 0.78, 0.78]

# Numerics
formal_solver: DELO_linear
use_1_5D_approx: false
use_prec: true

solver:
  ksp_solver_type: KSPFGMRES
  ksp_rtol: 1e-8
  ksp_max_it: 1000

prec:
  pc_solver_type: KSPGMRES
  pc_rtol: 1e-4
  pc_max_it: 1000
  pc_use_J_KQ: true
  pc_formal_solver_approx: true
```

More examples can be found in `bin/config.yml` and in the test configurations
under `tests/`.
