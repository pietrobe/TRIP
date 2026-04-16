#include <stdio.h>
#include <stdlib.h>
#include "hdf_output_ar_field.h"
#include "hdf_output_example_3D.h"
#include "hdf_output_field.h"

/////////////////////////////////////////////////////////////
/// make_example_output_field_surface
/////////////////////////////////////////////////////////////
THDF_field_t *
make_example_output_field_surface(int N_frequencies, const int N_x, const int N_y, const int N_incl,
                                  const int N_azimuth) {

  const int total_size = N_x * N_y * N_incl * N_azimuth;

  THDF_field_t *output_field = (THDF_field_t *)malloc(total_size * sizeof(THDF_field_t));
  if (output_field == NULL) {
    fprintf(stderr, "Error allocating memory for output field surface\n");
    return NULL;
  }

  for (int i = 0; i < total_size; i++) {
    output_field[i].index_i       = i % N_x;
    output_field[i].index_j       = (i / N_x) % N_y;
    output_field[i].index_incl    = (i / (N_x * N_y)) % N_incl;
    output_field[i].index_azimuth = (i / (N_x * N_y * N_incl)) % N_azimuth;
    output_field[i].float_type    = HDF_OUT_FLOAT64;

    output_field[i].stokes_I  = malloc(N_frequencies * sizeof(double));
    output_field[i].stokes_QI = malloc(N_frequencies * sizeof(double));
    output_field[i].stokes_UI = malloc(N_frequencies * sizeof(double));
    output_field[i].stokes_VI = malloc(N_frequencies * sizeof(double));

    double *si  = (double *)output_field[i].stokes_I;
    double *sqi = (double *)output_field[i].stokes_QI;
    double *sui = (double *)output_field[i].stokes_UI;
    double *svi = (double *)output_field[i].stokes_VI;
    for (int f = 0; f < N_frequencies; f++) {
      si[f]  = 1.0 * (i + 1) * (f + 1);
      sqi[f] = 1.0 * (i + 1) * (f + 1);
      sui[f] = 1.0 * (i + 1) * (f + 1);
      svi[f] = 1.0 * (i + 1) * (f + 1) * 1000;
    }
  }

  return output_field;
}

/////////////////////////////////////////////////////////////
/// destroy_example_output_field_surface
/////////////////////////////////////////////////////////////
int
destroy_example_output_field_surface(THDF_field_t *output_field, const int N_x, const int N_y, const int N_incl,
                                     const int N_azimuth) {
  if (output_field == NULL) {
    return EXIT_FAILURE;
  }

  const int total_size = N_x * N_y * N_incl * N_azimuth;

  for (int i = 0; i < total_size; i++) {
    free(output_field[i].stokes_I);
    free(output_field[i].stokes_QI);
    free(output_field[i].stokes_UI);
    free(output_field[i].stokes_VI);
  }

  free(output_field);
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////
int
main_example_ar(int argc, char **argv) {

  (void)argc;
  (void)argv;

  // First, create and write the file

  int N_x = 4;
  int N_y = 4;

  int     N_frequencies = 55;
  double *frequencies   = (double *)malloc(N_frequencies * sizeof(double));

  for (int i = 0; i < N_frequencies; i++) {
    frequencies[i] = 1.0e14 + i * 1.0e12;
  }

#define N_DIRECTIONS 5
  double muv[N_DIRECTIONS]  = {0.1, 0.2, 0.3, 0.4, 0.5};
  double chiv[N_DIRECTIONS] = {0.0};

  hid_t file_id = H5Fcreate("output_field_ar.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // write frequencies grid
  THDF_ard_frequencies_t ard_freqs;
  ard_freqs.N_frequencies = N_frequencies;
  ard_freqs.frequencies   = frequencies;
  THDF_write_ard_frequencies(file_id, &ard_freqs);

  // write arbitrary directions
  THDF_ard_directions_t ard_directions;
  ard_directions.N_directions = N_DIRECTIONS;
  ard_directions.mu           = muv;
  ard_directions.chi          = chiv;
  THDF_write_ard_directions(file_id, &ard_directions);

  for (int di = 0; di < N_DIRECTIONS; di++) {

    THDF_ard_field_t field_data;

    {
      THDF_ard_field_handler *dset_handler =                 //
          THDF_create_ard_field_handler_mpi(file_id,         //
                                            muv[di],         //
                                            chiv[di],        //
                                            N_x,             //
                                            N_y,             //
                                            N_DIRECTIONS,    //
                                            N_frequencies);  //
      for (int ix = 0; ix < N_x; ix++) {
        for (int jy = 0; jy < N_y; jy++) {
          field_data.index_i       = ix;
          field_data.index_j       = jy;
          field_data.index_k       = 0;
          field_data.beam_index    = di;
          field_data.N_frequencies = N_frequencies;
          field_data.stokes_I      = (THDF_float_t *)malloc(N_frequencies * sizeof(THDF_float_t));
          field_data.stokes_QI     = (THDF_float_t *)malloc(N_frequencies * sizeof(THDF_float_t));
          field_data.stokes_UI     = (THDF_float_t *)malloc(N_frequencies * sizeof(THDF_float_t));
          field_data.stokes_VI     = (THDF_float_t *)malloc(N_frequencies * sizeof(THDF_float_t));

          for (int f = 0; f < N_frequencies; f++) {
            field_data.stokes_I[f]  = 1.0 * (ix + 1) * (jy + 1) * (f + 1);
            field_data.stokes_QI[f] = 1.0 * (ix + 1) * (jy + 1) * (f + 1);
            field_data.stokes_UI[f] = 1.0 * (ix + 1) * (jy + 1) * (f + 1);
            field_data.stokes_VI[f] = 1.0 * (ix + 1) * (jy + 1) * (f + 1) * 1000;
          }

          THDF_write_ard_field_to_hdf5(dset_handler, &field_data, ix, jy, 1, 1, N_frequencies);

          free(field_data.stokes_I);
          free(field_data.stokes_QI);
          free(field_data.stokes_UI);
          free(field_data.stokes_VI);
        }
      }
      THDF_close_ard_field_handler_mpi(dset_handler);
    }
  }

  H5Fclose(file_id);

  free(frequencies);

  return 0;
}

/////////////////////////////////////////////////////////////
/// main
/////////////////////////////////////////////////////////////
int
main(int argc, char **argv) {

  return main_example_ar(argc, argv);

  (void)argc;
  (void)argv;

  // First, create and write the file

  int     N_frequencies = 13;
  double *frequencies   = NULL;

  int     N_x                  = 4;
  int     N_y                  = 4;
  int     N_theta              = 6;
  int     N_chi                = 8;
  int     total_entries        = N_x * N_y * N_theta * N_chi;
  double *theta                = NULL;
  double *chi                  = NULL;
  int    *inclinations_indices = NULL;
  int    *azimuthal_indices    = NULL;

  theta                = (double *)malloc(N_theta * N_chi * sizeof(double));
  chi                  = (double *)malloc(N_theta * N_chi * sizeof(double));
  inclinations_indices = (int *)malloc(N_theta * N_chi * sizeof(int));
  azimuthal_indices    = (int *)malloc(N_theta * N_chi * sizeof(int));

  THDF_field_t *output_field = make_example_output_field_surface(N_frequencies, N_x, N_y, N_theta, N_chi);

  for (int i = 0; i < N_theta; i++) {
    for (int j = 0; j < N_chi; j++) {
      theta[i * N_chi + j]                = (i + 1) * 10.0;  // just example values
      chi[i * N_chi + j]                  = j * 45.0;        // just example values
      inclinations_indices[i * N_chi + j] = i;
      azimuthal_indices[i * N_chi + j]    = j;
    }
  }

  frequencies = (double *)malloc(N_frequencies * sizeof(double));

  for (int i = 0; i < N_frequencies; i++) {
    frequencies[i] = 1.0e14 + i * 1.0e12;
  }

  // Create frequencies grid struct
  THDF_frequencies_grid_t freq_grid;
  freq_grid.N_frequencies = N_frequencies;
  freq_grid.frequencies   = frequencies;

  hid_t file = H5Fcreate("output_field.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) {
    fprintf(stderr, "Error creating HDF5 file\n");
    return EXIT_FAILURE;
  }

  if (THDF_write_frequencies_grid_to_hdf5(file, &freq_grid) < 0) {
    fprintf(stderr, "Error writing frequencies grid data\n");
    H5Fclose(file);
    return EXIT_FAILURE;
  }

  // Create angular grid struct
  THDF_angular_grid_t angular_grid;
  angular_grid.N_directions         = N_theta * N_chi;
  angular_grid.N_inclination_angles = N_theta;
  angular_grid.N_azimuthal_angles   = N_chi;
  angular_grid.azimuthal_angles     = chi;
  angular_grid.inclination_angles   = theta;
  angular_grid.azimuthal_indices    = azimuthal_indices;
  angular_grid.inclinations_indices = inclinations_indices;
  if (THDF_write_angular_grid_to_hdf5(file, &angular_grid) < 0) {
    fprintf(stderr, "Error writing angular grid data\n");
    H5Fclose(file);
    return EXIT_FAILURE;
  }

  // Create output field dataset
  THDF_field_handler_t *output_dset =
      THDF_create_field_handler_mpi(file, N_x, N_y, N_theta, N_chi, N_frequencies, HDF_OUT_FLOAT64);

  for (int i = 0; i < total_entries; i++) {
    if (THDF_write_field_dataset_to_hdf5(output_dset, &output_field[i], output_field[i].index_i, output_field[i].index_j,
                                         output_field[i].index_incl, output_field[i].index_azimuth, 1, 1, 1, 1,
                                         N_frequencies) < 0) {
      fprintf(stderr, "Error writing output field surface data\n");
      H5Fclose(file);
      return EXIT_FAILURE;
    }
  }

  // Clean up
  THDF_close_field_handler_mpi(output_dset);

  H5Fclose(file);

  destroy_example_output_field_surface(output_field, N_x, N_y, N_theta, N_chi);

  free(frequencies);
  free(theta);
  free(chi);
  free(inclinations_indices);
  free(azimuthal_indices);
  return 0;
}