#ifndef HDF_OUTPUT_EXAMPLE_ZSLICE_3D_H
#define HDF_OUTPUT_EXAMPLE_ZSLICE_3D_H

#if defined(__cplusplus)
extern "C" {
#endif

#include "hdf_output_zslice_3D.h"

/* =========================================================================
 * Public Functions
 * ========================================================================= */

/**
 * main_3d_example_zslice
 *
 * Main entry point for z-slab HDF5 writing example.
 *
 * Orchestrates the complete workflow:
 *   1. Create slab files (collective on slab_comm)
 *   2. Write grid metadata (serial on slab_rank 0)
 *   3. Open persistent handler (collective)
 *   4. Write z-blocks in loop (collective DXPL per write)
 *   5. Close handler and file
 *   6. Build VDS master (rank 0 only)
 *
 * Args:
 *   argc: Command line argument count
 *   argv: Command line arguments
 *
 * Returns:
 *   0 on success, -1 on failure
 *
 * Environment variables:
 *   HDF_N_X              Grid dimension X (default: 8)
 *   HDF_N_Y              Grid dimension Y (default: 8)
 *   HDF_N_Z              Grid dimension Z (default: 12)
 *   HDF_N_FREQUENCIES    Number of frequencies (default: 96)
 *   HDF_N_INCL           Number of inclination angles (default: 16)
 *   HDF_N_AZIMUTH        Number of azimuthal angles (default: 8)
 *   HDF_N_STOKES         Number of Stokes parameters (default: 4)
 *   HDF_NORMALIZE_OUTPUT Enable normalization output (default: 0)
 *   HDF_STEP_Z           Z block size for writing (default: 15)
 *   HDF_OUTPUT_DIR       Output directory (default: ".")
 */
int
main_3d_example_zslice(int argc, char **argv);

int
main_3d_example_multi_zslice(int argc, char **argv);

int
main_3d_example_single_file(int argc, char **argv);

#if defined(__cplusplus)
}
#endif

#endif /* HDF_OUTPUT_EXAMPLE_ZSLICE_3D_H */
