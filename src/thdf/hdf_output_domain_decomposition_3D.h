#ifndef __THDF_OUTPUT_DOMAIN_DECOMPOSITION_3D_H__
#define __THDF_OUTPUT_DOMAIN_DECOMPOSITION_3D_H__

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Factorizes 'size' into (Px, Py, Pz) such that Px*Py*Pz == size
 * and the processor grid matches the domain aspect ratio as closely as possible.
 */
void
factorize_procs(int size, int Nx, int Ny, int Nz, int *Px, int *Py, int *Pz);

/**
 * Given the processor grid (Px, Py, Pz), returns the rank of a neighbor.
 * Returns MPI_PROC_NULL for out-of-bounds neighbors (domain boundary).
 *
 * Usage example — find the +X neighbor:
 *   int nbr = get_neighbor_rank(Px, Py, Pz, ix+1, iy, iz);
 */
int
get_neighbor_rank(int Px, int Py, int Pz, int ix, int iy, int iz);

void
make_input_example_dataset_zz(double *stokes_IQUI, int N_x, int N_y, int N_z, int N_incl, int N_azimuth,
                              int N_frequencies, int N_stokes, int *stride_x, int *stride_y, int *stride_z,
                              int *stride_incl, int *stride_azimuth, int *stride_frequencies, int *stride_stokes,
                              int z_par);

/**
 * 3D Domain Decomposition —
 *
 * Maps each MPI rank to a unique (ix, iy, iz) coordinate in a
 * Px x Py x Pz processor grid using simple row-major indexing:
 *
 *   rank = ix * (Py * Pz) + iy * Pz + iz
 *
 * Inputs:
 *   N_x, N_y, N_z  — global domain extents
 *   comm            — MPI communicator
 *
 * Outputs:
 *   N_local_{x,y,z}      — size of this rank's subdomain
 *   local_start_{x,y,z}  — global start index of this rank's subdomain
 *   P{x,y,z}             — processor grid dimensions (optional, can be NULL)
 */
void
decompose_domain_3d(int mpi_rank, int mpi_size,                                  //
                    int N_x, int N_y, int N_z,                                   //
                    int *N_local_x, int *N_local_y, int *N_local_z,              //
                    int *local_start_x, int *local_start_y, int *local_start_z,  //
                    int *out_Px, int *out_Py, int *out_Pz);                      //

void
demo_decompose_domain_3d(int argc, char **argv);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // __THDF_OUTPUT_DOMAIN_DECOMPOSITION_3D_H__
