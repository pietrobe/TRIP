# TRIP
Three-dimensional Radiative transfer Including Polarization and PRD

### Dependencies
* PETSc
* rii
* Kokkos
* sgrid

### Input
The input is encoded in the `bin/config.yml` file, to be changed as necessary.
The scattering module can be changed via the `emissivity_model` field, between the following options:

| `emissivity_model`  | Description                                                        |
|----------------|--------------------------------------------------------------------|
| **NONE**       | undefined                                                          |
| **CRD_limit**  | CRD limit (default set to CRD_GL)                                  |
| **CRD_limit_VHP** | CRD limit with VHP (very high precision) approximation         |
| **PRD**        | default (PRD_FAST), partial redistribution (default grid)           |
| **PRD_NORMAL** | force PRD with original grid                                       |
| **PRD_FAST**   | force PRD with fast grid                                           |
| **PRD_AA**     | angle-averaged PRD method                                          |
| **PRD_AA_MAPV** | angle-averaged PRD storing a value map (**high memory usage**)   |
| **ZERO**       | continuum                                                          |


### Run
```bash
srun ./solar_3D
```

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

### Compile
```bash
cd bin
cmake .. -Dsgrid_DIR=/path_to_sgrid/sgrid_installation/lib/cmake/ -DRII_ROOT_PATH=/path_to_rii/rii/rii-c/
make
```

### Compile with accelerator support
```bash
cd bin
cmake .. -Dsgrid_DIR=/path_to_sgrid/sgrid_installation/lib/cmake/ \
         -DRII_ROOT_PATH=/path_to_rii/rii/rii-c/ \
         -DUSE_ACC=ON;
make
```

For ACC support, make sure that the RII library is compiled with ACC support as well.
First build the RII accelerated module in the RII source directory:

```bash 
cd ${RII_MAIN_DIR}/rii-c/src_acc;
make clean;
make cuda_arch=90 -j12;
```
The ```cuda_arch=90``` works for the NVIDIA H100 GPUs on the MareNostrum 5 ACC and Alps Daint systems.
Modify the cuda_arch flag as needed for other GPU architectures.


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


### Contributing guidelines

**TRIP** is developed by [Pietro Benedusi](https://pietrobe.github.io/) and Simone Riva in Luca Belluzzi's group at [IRSOL](https://www.irsol.usi.ch/).


