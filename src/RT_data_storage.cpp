#include <iostream>
#include <iterator>
#include <limits>
#include <mpi.h>

#include "RT_data_storage.hpp"

///////////////////////////////////////////////////////////////////////////////
namespace RT_data_storage {

///////////////////////////////////////////////////////////////////////////////
// normalizations_values
///////////////////////////////////////////////////////////////////////////////
std::tuple<double, double, double, double>
normalizations_values(const double *field, const int size_theta,
                      const int size_chi, const int size_u) {

  double I_max = std::numeric_limits<double>::min();
  double Q_max = std::numeric_limits<double>::min();
  double U_max = std::numeric_limits<double>::min();
  double V_max = std::numeric_limits<double>::min();

  const auto theta_off = 4 * size_chi * size_u;
  const auto chi_off = 4 * size_u;
  const auto u_off = 4;

  for (int i = 0; i < size_theta; i++) {
    for (int j = 0; j < size_chi; j++) {
      for (int k = 0; k < size_u; k++) {
        const int index = i * theta_off + j * chi_off + k * u_off;
        I_max = std::max(std::abs(I_max), field[index]);
        Q_max = std::max(std::abs(Q_max), field[index + 1]);
        U_max = std::max(std::abs(U_max), field[index + 2]);
        V_max = std::max(std::abs(V_max), field[index + 3]);
      }
    }
  }

  const long double I_norm_factor = 1.0L / (long double)(I_max);
  const long double Q_norm_factor = 1.0L / (long double)(Q_max);
  const long double U_norm_factor = 1.0L / (long double)(U_max);
  const long double V_norm_factor = 1.0L / (long double)(V_max);

  return std::make_tuple(I_norm_factor, Q_norm_factor,  //
                         U_norm_factor, V_norm_factor); //
}

///////////////////////////////////////////////////////////////////////////////
// write_fields
///////////////////////////////////////////////////////////////////////////////
bool write_fields(
    const std::filesystem::path &output_dir,
    const std::list<double *> &fields_list,
    const std::list<std::tuple<int, int, int>> &fields_cooordinates_list,
    const int size_theta, const int size_chi, const int size_u,
    const bool normalize_fields, const bool to_float32) {

  const int fileds_cnt = fields_list.size();

  if (fileds_cnt != fields_cooordinates_list.size()) {
    std::cerr << "ERROR: File: " << __FILE__ << " Line: " << __LINE__ << ": "
              << "The number of fields and the number of coordinates do not "
                 "match\n";

    return false;
  }

  // open mpi file
  MPI_File fh;
  MPI_File_open(MPI_COMM_WORLD, output_dir.c_str(),
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

  // write the number of fields
  MPI_File_write(fh, &fileds_cnt, 1, MPI_INT, MPI_STATUS_IGNORE);

  for (int field_i = 0; field_i < fileds_cnt; field_i++) {

    auto field_it = fields_list.begin();
    auto field_coord_it = fields_cooordinates_list.begin();

    std::advance(field_it, field_i);
    std::advance(field_coord_it, field_i);

    const int abs_i = std::get<0>(*field_coord_it);
    const int abs_j = std::get<1>(*field_coord_it);
    const int abs_k = std::get<2>(*field_coord_it);

    const auto norm_values =
        normalize_fields
            ? normalizations_values(*field_it, size_theta, size_chi, size_u)
            : std::make_tuple(1.0, 1.0, 1.0, 1.0);

    // write abs_i, abs_j, abs_k
    MPI_File_write(fh, &abs_i, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_write(fh, &abs_j, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_write(fh, &abs_k, 1, MPI_INT, MPI_STATUS_IGNORE);

    // write the normalization values
    MPI_File_write(fh, &std::get<0>(norm_values), 1, MPI_DOUBLE,
                   MPI_STATUS_IGNORE);
    MPI_File_write(fh, &std::get<1>(norm_values), 1, MPI_DOUBLE,
                   MPI_STATUS_IGNORE);
    MPI_File_write(fh, &std::get<2>(norm_values), 1, MPI_DOUBLE,
                   MPI_STATUS_IGNORE);
    MPI_File_write(fh, &std::get<3>(norm_values), 1, MPI_DOUBLE,
                   MPI_STATUS_IGNORE);

    // write if double or float
    const int to_float32_ = to_float32 ? 1 : 0;
    MPI_File_write(fh, &to_float32_, 1, MPI_INT, MPI_STATUS_IGNORE);

    // write the field
    const auto theta_off = 4 * size_chi * size_u;
    const auto chi_off = 4 * size_u;
    const auto u_off = 4;

    const int field_size = size_theta * size_chi * size_u * 4;

    std::vector<float> field_float32(field_size);
    std::vector<double> field_double(field_size);

    for (int i = 0; i < size_theta; i++) {
      for (int j = 0; j < size_chi; j++) {
        for (int k = 0; k < size_u; k++) {
          const int index = i * theta_off + j * chi_off + k * u_off;

          const double I = (*field_it)[index] * std::get<0>(norm_values);
          const double Q = (*field_it)[index + 1] * std::get<1>(norm_values);
          const double U = (*field_it)[index + 2] * std::get<2>(norm_values);
          const double V = (*field_it)[index + 3] * std::get<3>(norm_values);

          if (to_float32) {
            field_float32[index] = (float)I;
            field_float32[index + 1] = (float)Q;
            field_float32[index + 2] = (float)U;
            field_float32[index + 3] = (float)V;
          } else {
            field_double[index] = I;
            field_double[index + 1] = Q;
            field_double[index + 2] = U;
            field_double[index + 3] = V;
          }
        }
      }
    }

    if (to_float32) {
      MPI_File_write(fh, field_float32.data(), field_size, MPI_FLOAT,
                     MPI_STATUS_IGNORE);
    } else {
      MPI_File_write(fh, field_double.data(), field_size, MPI_DOUBLE,
                     MPI_STATUS_IGNORE);
    }
  }

  MPI_File_close(&fh);

  return true;
}

} // namespace RT_data_storage