# README #

Parallel implementation for PRD radiative transfer in 3D

### Dependencies ### 
* MPI
* PETSc
* Kokkos
* rii
* sgrid

### Input Dependencies ###
* input files encoding the atmospheric model 

### Compile ###
cd bin  
cmake .. -Dsgrid_DIR=/path_to_sgrid/sgrid_installation/lib/cmake/ -DRII_ROOT_PATH=/path_to_rii/rii/rii-c/  
make 

### Run ###
srun ./solar_3D
