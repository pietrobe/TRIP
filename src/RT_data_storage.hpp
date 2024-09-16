#ifndef RT_DATA_STORAGE_HPP
#define RT_DATA_STORAGE_HPP

#include <filesystem>
#include <fstream>
#include <list>
#include <tuple>

///////////////////////////////////////////////////////////////////////////////
namespace RT_data_storage {

/**
 * @brief Calucalte the normalizations values for the field for I, Q, U, and V
 * respectively
 *
 * @param field
 * @param size_theta
 * @param size_chi
 * @param size_u
 * @return std::tuple<double, double, double, double> storing the normalizations
 * values for I, Q, U, and V respectively
 */
std::tuple<double, double, double, double>
normalizations_values(const double *field, const int size_theta,
                      const int size_chi, const int size_u);

/**
 * @brief Write the fields to a file
 *
 * @param output_dir
 * @param fields_list
 * @param fields_cooordinates_list
 * @param size_theta
 * @param size_chi
 * @param size_u
 * @param normalize_fields
 * @param to_float32
 * @return true
 * @return false
 */
bool write_fields(
    const std::filesystem::path &output_dir,
    const std::list<double *> &fields_list,
    const std::list<std::tuple<int, int, int>> &fields_cooordinates_list,
    const int size_theta, const int size_chi, const int size_u,
    const bool normalize_fields = true, 
    const bool to_float32 = true);

/**
 * @brief
 *
 * @param input_dir
 * @param fields_list
 * @param fields_cooordinates_list
 * @param size_theta
 * @param size_chi
 * @param size_u
 * @param normalize_fields
 * @return true
 * @return false
 */
bool read_fields(const std::filesystem::path &input_dir,
                 std::list<double *> &fields_list,
                 std::list<std::tuple<int, int, int>> &fields_cooordinates_list,
                 int &size_theta, int &size_chi, int &size_u,
                 const bool normalize_fields = false);

} // namespace RT_data_storage

#endif // RT_DATA_STORAGE_HPP