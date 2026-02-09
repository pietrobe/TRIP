# TRIP
Three-dimensional Radiative transfer Including Polarization and PRD

## Dependencies
* [PETSc](https://petsc.org/release/)
* rii    (Add link)

## Input
The input is encoded in the `bin/config.yml` file, to be changed as necessary.
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


### Cite 
```
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
## Compilation
### Prerequisite
Firstly, build *rii*:
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

For **accelerator support**, build the RII accelerated module in the RII source directory:
```bash
cd ${RII_MAIN_DIR}/rii-c/src_acc;
make clean;
make cuda_arch=90 -j12;
```
The `cuda_arch=90` works for the NVIDIA H100 GPUs on the MareNostrum 5 ACC and Alps Daint systems. Modify the cuda_arch flag as needed for other GPU architectures.

For the MareNostrum 5 ACC system, use the options:
```bash
make cuda_arch=90 CXX=/apps/ACC/GCC/11.4.0/bin/g++ -j12
```

Then:
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

Note to add others flags to the `cmake` command if needed.

### Standard build
```bash
cd bin
cmake .. -DRII_ROOT_PATH=/path_to_rii/rii/rii-c/
make
```

### Build with accelerator support
For ACC support, make sure that the RII library is compiled with ACC support as well.
First build the RII accelerated module in the RII source directory:

```bash 
cd ${RII_MAIN_DIR}/rii-c/src_acc;
make clean;
make cuda_arch=90 -j12;
```
The ```cuda_arch=90``` works for the NVIDIA H100 GPUs on the MareNostrum 5 ACC and Alps Daint systems.
Modify the cuda_arch flag as needed for other GPU architectures.

Then build TRIP:
```bash
cd bin
cmake .. -DRII_ROOT_PATH=/path_to_rii/rii/rii-c/ \
         -DUSE_ACC=ON;
make
```
For the MareNostrum 5 ACC system, use the options:
```bash
make cuda_arch=90 CXX=/apps/ACC/GCC/11.4.0/bin/g++ -j12
```

Then, in the *rii* configure step of RII, use the flag:
```bash

mkdir -p ${RII_MAIN_DIR}/rii-c/build/;
cd ${RII_MAIN_DIR}/rii-c/build/;

cmake .. -DDYNAMIC_LIBRARY=ON \
         -DENABLE_MPI=ON \
         -DGPU_AWARE_MPI=ON \
         -DUSE_ACCELERATOR=cuda \
         -DOUT_LOG_LEVEL=0 \
         -DRII_GLOBALJUMP=64;
```

add the others flags as needed, and then build the RII library:
```bash
make -j12
```

### Running TRIP
```bash
srun ./solar_3D
```

### Contributors

**TRIP** is developed by [Pietro Benedusi](https://pietrobe.github.io/), Simone Riva, and [Melanie Tonarelli](https://github.com/melanie-t27) in Luca Belluzzi's group at [IRSOL](https://www.irsol.usi.ch/).


