# TRIP
Three-dimensional Radiative transfer Including Polarization and PRD

## Citation
If you use TRIP in your research, please cite:
```bibtex
@article{benedusi2023scalable,
  title={Scalable matrix-free solver for 3D transfer of polarized radiation in stellar atmospheres},
  author={Benedusi, Pietro and Riva, Simone and Zulian, Patrick and {\v{S}}t{\v{e}}p{\'a}n, Ji{\v{r}}{\'\i} and Belluzzi, Luca and Krause, Rolf},
  journal={Journal of Computational Physics},
  volume={479},
  pages={112013},
  year={2023},
  publisher={Elsevier}
}
```

## Dependencies
* [PETSc](https://petsc.org/release/)
* [rii](https://gitlab.com/simon-r/rii) 
* MPI implementation (e.g., OpenMPI or MPICH)

Optional dependencies:
* CUDA (for accelerator support)

## Compilation and Execution
We assume the environment variable `TRIP_MAIN_DIR` points to the root directory of TRIP:
```bash
export TRIP_MAIN_DIR=/path/to/TRIP
```

### Standard build
After cloning the TRIP repository:
```bash
mkdir -p ${TRIP_MAIN_DIR}/build/
cd ${TRIP_MAIN_DIR}/build/
cmake .. -DRII_ROOT_PATH=/path_to_rii/rii/rii-c/
make
```

### Build with accelerator support
Ensure that the *rii* library is compiled with ACC support.
Then build TRIP:
```bash
mkdir -p ${TRIP_MAIN_DIR}/build/
cd ${TRIP_MAIN_DIR}/build/
cmake .. -DRII_ROOT_PATH=/path_to_rii/rii/rii-c/ \
         -DUSE_ACC=ON
make -j
```

## Troubleshooting

### PETSc not found
If PETSc is installed by the user and is not in a default system path, export the env variables before compilation:
```bash
export PETSC_DIR=${PETSC_MAIN_PATH}/petsc-install
export PKG_CONFIG_PATH=$PETSC_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
```

where `${PETSC_MAIN_PATH}` is the main installation directory of PETSc.

## Input
Input is provided via a YAML configuration file.
Example configuration files are available in `bin/` and can be changed as necessary.

The scattering module can be changed via the `emissivity_model` field, between the following options:

| `emissivity_model`  | Description                                                        |
|----------------|--------------------------------------------------------------------|
| **NONE**       | undefined                                                          |
| **CRD_limit**  | CRD limit (default set to CRD_GL)                                  |
| **CRD_limit_VHP** | CRD limit with VHP (very high precision) approximation         |
| **PRD**        | default (PRD_FAST), partial redistribution (default grid)           |
| **PRD_NORMAL** | PRD with original grid                                       |
| **PRD_MEDIUM** | PRD with the same or better accuracy w.r.t. PRD_NORMAL but faster          |
| **PRD_FAST**   | PRD with fast grid                                           |
| **PRD_AA**     | angle-averaged PRD method                                          |
| **PRD_AA_MAPV** | angle-averaged PRD storing a value map (**high memory usage**), 2–3× faster   |
| **ZERO**       | continuum                                                          |


## Running TRIP
The following command assumes a SLURM environment. Replace `srun` with `mpirun` or `mpiexec` on non-SLURM systems.
```bash
srun ./solar_3D --yaml_config ${PATH_TO_YAML_CONFIG}
```

## Output files

TRIP produces HDF5 (*.h5) files by default, as well as CSV and MATLAB data script files as options.

**Large output (full 3D radiation field):**
* Set the scratch file system as the output directory.
* In the sbatch script, source the bash scripts in `${TRIP_MAIN_DIR}/bin/pfs_env` to optimize the writing parameters for a specific parallel file system.

**Lustre file system** (Capstor at CSCS, Shaheen III at KAUST): set the stripes of the output directory **before** executing TRIP:
```bash
lfs setstripe -c 16 -S 4m -i -1 ${BASE_OUTPUT_DIR}
```
The setstripe is applied recursively to new subdirectories. Without this step, writing the full 3D radiation field can take over an hour.

**GPFS file system** (MareNostrum 5 ACC|GPP): stripes are set by default; no preparation required.

## Accelerating PRD_AA with local scratch (MAPV)

The `PRD_AA_MAPV` emissivity model saves a map of values (MAPV) to avoid recomputation, making it 2–3× faster than `PRD_AA`. To enable this, set `emissivity_model: PRD_AA_MAPV` in the YAML config file and export the following in the SLURM script:
```bash
export RII_AA_SCRATCH_DIR=$PATH_TO_THE_LOCAL_SCRATCH
export RII_AA_JOB_PATTERN="TRIP_job_$SLURM_JOB_ID"
```
On **MareNostrum 5**:
```bash
export RII_AA_SCRATCH_DIR=$TMPDIR
export RII_AA_JOB_PATTERN="TRIP_job_$SLURM_JOB_ID"
```
On **Alps CSCS**: local scratch is not available; using global scratch is not recommended.

If these variables are not set, the MAPV is stored in RAM, increasing the risk of OOM.


## Testing
```bash
cd ${TRIP_MAIN_DIR}/build/
make -j 20 build_test
../bin/run_test.sh
```

## Contributors

**TRIP** is developed by [Pietro Benedusi](https://pietrobe.github.io/), [Simone Riva](https://gitlab.com/simon-r/rii), and [Melanie Tonarelli](https://github.com/melanie-t27) in Luca Belluzzi's group at [IRSOL](https://www.irsol.usi.ch/).
