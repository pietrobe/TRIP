# TRIP
Three-dimensional Radiative transfer Including Polarization and PRD

### Dependencies
* PETSc
* rii
* Kokkos
* sgrid

### Input
The input is encoded in the `bin/config.yml` file, to be changed as necessary.

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

### Run
```bash
srun ./solar_3D
```
