# Test executed on different cluster

## MN5 GPP 

```text
128 x 128 x 128
[t] open handler: 2.553 s
[t] z=0 nz=2: 5.496 s
[t] local written: 6.000 GB  approx. peak-bandwidth: 1.092 GB/s
[t] z=2 nz=2: 1.919 s
[t] local written: 6.000 GB  approx. peak-bandwidth: 3.126 GB/s
[t] write loop: 7.455 s
[t] total written: 384.0 GB  effective loop-bandwidth: 51.5 GB/s
[t] H5Fflush: 0.248 s
[t] H5Fclose: 0.018 s
[t] total: 25.538 s
[t] total size: 384.0 GB
[MPI] max local z-slab: 4  min local z-slab: 3
[MPI] decomposition: (32, 32, 40)
[MPI] mpi_size: 40960
[t] output: 384.0 GB  TTS-bandwidth: 15.0 GB/s
Output files:
-rw-r--r-- 1 usi441290 ehpc597 385G Mar 22 09:46 /gpfs/scratch/ehpc597/hdf5_output//3d_field/radiation_field_3d.h5
```

```text
export HDF_N_X=192
export HDF_N_Y=192
export HDF_N_Z=128
[t] open handler: 2.444 s
[t] z=0 nz=2: 7.883 s
[t] local written: 13.500 GB  approx. peak-bandwidth: 1.712 GB/s
[t] z=2 nz=2: 4.399 s
[t] local written: 13.500 GB  approx. peak-bandwidth: 3.069 GB/s
[t] write loop: 12.365 s
[t] total written: 864.0 GB  effective loop-bandwidth: 69.9 GB/s
[t] H5Fflush: 0.111 s
[t] H5Fclose: 0.023 s
[t] total: 42.363 s
[t] total size: 864.0 GB
[MPI] max local z-slab: 4  min local z-slab: 4
[MPI] decomposition: (32, 40, 32)
[MPI] mpi_size: 40960
[t] output: 864.0 GB  TTS-bandwidth: 20.4 GB/s
Output files:
-rw-r--r-- 1 usi441290 ehpc597 865G Mar 22 09:54 /gpfs/scratch/ehpc597/hdf5_output//3d_field/radiation_field_3d.h5
```

```text
export HDF_N_X=256
export HDF_N_Y=256
export HDF_N_Z=128
[t] open handler: 2.239 s
[t] z=0 nz=2: 16.556 s
[t] local written: 24.000 GB  approx. peak-bandwidth: 1.450 GB/s
[t] z=2 nz=2: 16.884 s
[t] local written: 24.000 GB  approx. peak-bandwidth: 1.421 GB/s
[t] z=4 nz=2: 15.091 s
[t] local written: 24.000 GB  approx. peak-bandwidth: 1.590 GB/s
[t] z=6 nz=1: 30.688 s
[t] local written: 12.000 GB  approx. peak-bandwidth: 0.391 GB/s
[t] write loop: 79.360 s
[t] total written: 1536.0 GB  effective loop-bandwidth: 19.4 GB/s
[t] H5Fflush: 0.081 s
[t] H5Fclose: 0.015 s
[t] total: 277.302 s
[t] total size: 1536.0 GB
[MPI] max local z-slab: 7  min local z-slab: 6
[MPI] decomposition: (64, 32, 20)
[MPI] mpi_size: 40960
[t] output: 1536.0 GB  TTS-bandwidth: 5.5 GB/s
Output files:
-rw-r--r-- 1 usi441290 ehpc597 1.6T Mar 22 10:02 /gpfs/scratch/ehpc597/hdf5_output//3d_field/radiation_field_3d.h5
```

## CSCS Alps

```text
[t] open handler: 0.213 s
[t] z=0 nz=4: 32.179 s
[t] local written: 12.000 GB  approx. peak-bandwidth: 0.373 GB/s
[t] z=4 nz=4: 28.316 s
[t] local written: 12.000 GB  approx. peak-bandwidth: 0.424 GB/s
[t] write loop: 63.754 s
[t] total written: 384.0 GB  effective loop-bandwidth: 6.0 GB/s
[t] H5Fflush: 2.199 s
[t] H5Fclose: 0.007 s
[t] total: 170.597 s
[sizes] N_x=128  N_y=128  N_z=128  N_incl=16  N_azimuth=8  N_frequencies=96  N_stokes=4
[t] total size: 384.0 GB
[MPI] max local z-slab: 8  min local z-slab: 8
[MPI] decomposition: (8, 8, 16) for domain (128, 128, 128)
[MPI] mpi_size: 1024
[t] output: 384.0 GB  TTS-bandwidth: 2.3 GB/s
```
