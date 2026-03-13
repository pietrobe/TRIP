# TRIP
Three-dimensional Radiative transfer Including Polarization and PRD

### Cite 
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
We assume the following environment variables:
- `TRIP_MAIN_DIR` — root directory of TRIP
- `RII_MAIN_DIR` — root directory of rii

### Prerequisite: Build rii
Clone the *rii* repository, then:
```bash
mkdir -p ${RII_MAIN_DIR}/rii-c/build/;
cd ${RII_MAIN_DIR}/rii-c/build/;

cmake .. -DDYNAMIC_LIBRARY=ON \
         -DUSE_RII_CONTRIB_V4=ON \
         -DOUT_LOG_LEVEL=0 \
         -DDEBUG=OFF \
         -DUSE_KERNEL_DIRECT_VECT=ON \
         -DCMAKE_CXX_COMPILER=g++ \
         -DRII_GLOBALJUMP=8;
make -j
```

**Accelerator support (optional)**

Build the *rii* accelerated module:
```bash
cd ${RII_MAIN_DIR}/rii-c/src_acc;
make clean;
make cuda_arch=90 -j12;
```
The `cuda_arch=90` option targets the NVIDIA H100 GPUs (e.g. on the MareNostrum 5 ACC and Alps Daint systems). Adjust this flag for other GPU architectures.

For the MareNostrum 5 ACC system, use:
```bash
make cuda_arch=90 CXX=/apps/ACC/GCC/11.4.0/bin/g++ -j12
```

Then rebuild *rii* with accelerator support:
```bash
mkdir -p ${RII_MAIN_DIR}/rii-c/build/;
cd ${RII_MAIN_DIR}/rii-c/build/;

cmake .. -DDYNAMIC_LIBRARY=ON \
         -DENABLE_MPI=ON \
         -DGPU_AWARE_MPI=ON \
         -DUSE_ACCELERATOR=cuda \
         -DOUT_LOG_LEVEL=0 \
         -DRII_GLOBALJUMP=64;
make -j
```

Add additional `cmake` flags if needed.

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
| **PRD_NORMAL** | force PRD with original grid                                       |
| **PRD_MEDIUM** | force PRD with the same or better accuracy w.r.t. PRD_NORMAL but faster          |
| **PRD_FAST**   | force PRD with fast grid                                           |
| **PRD_AA**     | angle-averaged PRD method                                          |
| **PRD_AA_MAPV** | angle-averaged PRD storing a value map (**high memory usage**), 2 .. 3 times faster   |
| **ZERO**       | continuum                                                          |


## Running TRIP
```bash
srun ./solar_3D --yaml_config ${PATH_TO_YAML_CONFIG}
```

## Testing
```bash
cd ${TRIP_MAIN_DIR}/build/;
make -j 20 build_test;
../bin/run_test.sh;
```

# Contributors

**TRIP** is developed by [Pietro Benedusi](https://pietrobe.github.io/), Simone Riva, and [Melanie Tonarelli](https://github.com/melanie-t27) in Luca Belluzzi's group at [IRSOL](https://www.irsol.usi.ch/).

