#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "hdf_output_domain_decomposition_3D.h"

int *
prime_factors(int n, int *num_factors) {
  int *factors = NULL;
  int  count   = 0;

  for (int i = 2; i <= n / i; i++) {
    while (n % i == 0) {
      factors = (int *)realloc(factors, (count + 1) * sizeof(int));
      factors[count++] = i;
      n /= i;
    }
  }

  if (n > 1) {
    factors = (int *)realloc(factors, (count + 1) * sizeof(int));
    factors[count++] = n;
  }

  *num_factors = count;
  return factors;
}

static int
read_env_int_or_default(const char *name, int default_value) {
  const char *value = getenv(name);
  if (value == NULL || *value == '\0') {
    return default_value;
  }

  char *end    = NULL;
  long  parsed = strtol(value, &end, 10);
  if (end == value || *end != '\0' || parsed <= 0 || parsed > INT_MAX) {
    return default_value;
  }

  return (int)parsed;
}

//////////////////////////////////////////////////////
// free_output_field
//////////////////////////////////////////////////////
void
factorize_procs(int size, int Nx, int Ny, int Nz, int *Px, int *Py, int *Pz) {
  int    best_px = 1, best_py = 1, best_pz = size;
  double best_score = 1e18;

  for (int px = 1; px <= size; px++) {
    if (size % px != 0)
      continue;
    for (int py = 1; py <= size / px; py++) {
      if ((size / px) % py != 0)
        continue;
      int pz = size / (px * py);

      /* Local subdomain extents for this factorization */
      double lx = (double)Nx / px;
      double ly = (double)Ny / py;
      double lz = (double)Nz / pz;

      /* Score: sum of pairwise aspect ratios (minimum = 6.0 for a cube) */
      double score = (lx / ly + ly / lx) + (ly / lz + lz / ly) + (lx / lz + lz / lx);

      // printf("  (%d x %d x %d)  local=(%.1f x %.1f x %.1f)  score=%.3f%s\n", px, py, pz, lx, ly, lz, score,
      //        score < best_score ? "  <-- new best" : "");

      if (score < best_score) {
        best_score = score;
        best_px    = px;
        best_py    = py;
        best_pz    = pz;
      }
    }
  }

  *Px = best_px;
  *Py = best_py;
  *Pz = best_pz;
}

///////////////////////////////////////////////////
// get_1d_decomposition
///////////////////////////////////////////////////
static void
get_1d_decomposition(int n_global, int p_dims, int p_coord, int *n_local, int *start) {
  int base      = n_global / p_dims;
  int remainder = n_global % p_dims;

  if (p_coord < remainder) {
    *n_local = base + 1;
    *start   = p_coord * (base + 1);
  } else {
    *n_local = base;
    *start   = remainder * (base + 1) + (p_coord - remainder) * base;
  }
}

///////////////////////////////////////////////////////
// get_neighbor_rank
///////////////////////////////////////////////////////
int
get_neighbor_rank(int Px, int Py, int Pz, int ix, int iy, int iz) {
  if (ix < 0 || ix >= Px)
    return MPI_PROC_NULL;
  if (iy < 0 || iy >= Py)
    return MPI_PROC_NULL;
  if (iz < 0 || iz >= Pz)
    return MPI_PROC_NULL;
  return ix * (Py * Pz) + iy * Pz + iz;
}

///////////////////////////////////////////////////
// decompose_domain_3d
///////////////////////////////////////////////////
void
decompose_domain_3d(int mpi_rank, int mpi_size, int N_x, int N_y, int N_z, int *N_local_x, int *N_local_y, int *N_local_z,
                    int *local_start_x, int *local_start_y, int *local_start_z, int *out_Px, int *out_Py, int *out_Pz) {

  /* --- 1. Find best (Px, Py, Pz) factorization of 'size' --- */
  int Px, Py, Pz;
  factorize_procs(mpi_size, N_x, N_y, N_z, &Px, &Py, &Pz);

  if (mpi_rank == 0)
    printf("[decompose_domain_3d] grid: %d x %d x %d procs over %d x %d x %d domain\n", Px, Py, Pz, N_x, N_y, N_z);

  /* --- 2. Map linear rank -> 3D processor coordinate (row-major) ---
   *
   *   rank = ix*(Py*Pz) + iy*Pz + iz
   */
  int ix = mpi_rank / (Py * Pz);
  int iy = (mpi_rank % (Py * Pz)) / Pz;
  int iz = mpi_rank % Pz;

  /* --- 3. Compute local subdomain size and global start offset --- */
  get_1d_decomposition(N_x, Px, ix, N_local_x, local_start_x);
  get_1d_decomposition(N_y, Py, iy, N_local_y, local_start_y);
  get_1d_decomposition(N_z, Pz, iz, N_local_z, local_start_z);

  /* --- 4. Expose grid dims if caller needs them (e.g. to find neighbors) --- */
  if (out_Px)
    *out_Px = Px;
  if (out_Py)
    *out_Py = Py;
  if (out_Pz)
    *out_Pz = Pz;
}

void
demo_decompose_domain_3d(int argc, char **argv) {

  (void)argc;
  (void)argv;

  int mpi_rank, mpi_size;
  mpi_rank = 0;
  mpi_size = 8196;

  const int N_x = read_env_int_or_default("HDF_N_X", 128);
  const int N_y = read_env_int_or_default("HDF_N_Y", 128);
  const int N_z = read_env_int_or_default("HDF_N_Z", 128);

  int N_local_x, N_local_y, N_local_z;
  int local_start_x, local_start_y, local_start_z;
  int Px, Py, Pz;

  decompose_domain_3d(mpi_rank, mpi_size,                              //
                      N_x, N_y, N_z,                                   //
                      &N_local_x, &N_local_y, &N_local_z,              //
                      &local_start_x, &local_start_y, &local_start_z,  //
                      &Px, &Py, &Pz);                                  //

  printf("Rank %d/%d: local domain X[%d+%d] Y[%d+%d] Z[%d+%d]  grid: %d x %d x %d\n",   //
         mpi_rank, mpi_size,                                                            //
         local_start_x, N_local_x, local_start_y, N_local_y, local_start_z, N_local_z,  //
         Px, Py, Pz);                                                                   //
}
