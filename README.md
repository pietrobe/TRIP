# TRIP
Three-dimensional Radiative transfer Including Polarization and PRD

## Cite 
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
* rii    (Add link)
* MPI implementation (e.g., OpenMPI or MPICH)

Optional dependencies:
* CUDA (for accelerator support)

## Compilation and Execution
We assume the environment variables `TRIP_MAIN_DIR` as the root directory of TRIP.

### Standard build
After cloning the TRIP repository:
```bash
mkdir -p ${TRIP_MAIN_DIR}/build/;
cd ${TRIP_MAIN_DIR}/build/;
cmake .. -DRII_ROOT_PATH=/path_to_rii/rii/rii-c/;
make
```

### Build with accelerator support
Ensure that the *rii* library is compiled with ACC support.
Then build TRIP:
```bash
mkdir -p ${TRIP_MAIN_DIR}/build/;
cd ${TRIP_MAIN_DIR}/build/;
cmake .. -DRII_ROOT_PATH=/path_to_rii/rii/rii-c/ \
         -DUSE_ACC=ON;
make -j
```
#### troubleshooting (Petsc)
If Petsc is installed by the user and is not in a default system path before the compilation it is necessary to export the env variables:
```bash
export PETSC_DIR=${PETSC_MAIN_PATH}/petsc-install;
export PKG_CONFIG_PATH=$PETSC_DIR/lib/pkgconfig:$PKG_CONFIG_PATH;
```

where ```${PETSC_MAIN_PATH}``` is the main installation directory of Petsc.

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
| **PRD_AA_MAPV** | angle-averaged PRD storing a value map (**high memory usage**), 2 .. 3 times faster   |
| **ZERO**       | continuum                                                          |


## Running TRIP
```bash
srun ./solar_3D --yaml_config ${PATH_TO_YAML_CONFIG}
```

## Output files

TRIP produces HDM5 (*.h5) files by default, as well as CSV and MATLAB data script files as options.

#### IMPORTANT
In order to save the entire 3D radiation field: 
* Set up the scratch file system as the output directory. 
* In the sbatch script, source the bash scripts in ```${TRIP_MAIN_DIR}/bin/pfs_env``` to optimize the writing parameters for a specific parallel file system.


**ATTENTION**: 
* For the **Lustre file system**, it is necessary to set the stripes of the output directory **before** executing TRIP with the command line ```lfs setstripe -c 16 -S 4m -i -1 ${BASE_OUTPUT_DIR}``` (the setstripe is applied recursively in new subdirectories). This is the case for the Capstor FS at CSCS, or for the Shaheen III at KAUST. 
If this preventive preparation of the output direcory is not made the writing mehtods could require more than an hour to save the entire 3D radiation field.

* In the gpfs filesystem on the in the MareNostrum 5 (ACC|GPP) stripes are set by default in the system, so that, preventive preparation is not required.

## R_II Angle Averaged MAPV

In order to accelerate the RII Angle Averaged method that is used by the emissivity model `PRD_AA` or `PRD_AA_MAPV`, it is possible to save the MAPV (map of values) on the local scratch drive.
To enable this feature you must set
`emissivity_model: PRD_AA_MAPV` in the yaml config file, and in the slurm script it is necessary to set up the environment variables:
```bash
export RII_AA_SCRATCH_DIR=$PATH_TO_THE_LOCAL_SCRATCH
export RII_AA_JOB_PATTERN="TRIP_job_$SLURM_JOB_ID"
```
On **Mare Nostrum 5** it is:
```bash
export RII_AA_SCRATCH_DIR=$TMPDIR
export RII_AA_JOB_PATTERN="TRIP_job_$SLURM_JOB_ID"
```
On **Alps CSCS** the local scratch is not available, and we do not recommend using the global scratch.
If these two environment variables are not set, the MAPV is stored in RAM, increasing the risk of OOM.


## Testing
```bash
cd ${TRIP_MAIN_DIR}/build/;
make -j 20 build_test;
../bin/run_test.sh;
```

## Contributors

**TRIP** is developed by [Pietro Benedusi](https://pietrobe.github.io/), Simone Riva, and [Melanie Tonarelli](https://github.com/melanie-t27) in Luca Belluzzi's group at [IRSOL](https://www.irsol.usi.ch/).



