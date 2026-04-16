/**
 * output_utils.h
 *
 * Shared utility functions for HDF5 output operations.
 */

#if defined(__cplusplus)
extern "C" {
#endif

#ifndef OUTPUT_UTILS_H
#define OUTPUT_UTILS_H

#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502984
#endif

#ifndef PATH_MAX
#define PATH_MAX 8182
#endif

/**
 * read_env_int_or_default
 *
 * Read an integer environment variable, or return a default value if not set.
 *
 * @param name The environment variable name
 * @param def  The default value to return if not set
 * @return     The environment variable value, or def if not set
 */
int
read_env_int_or_default(const char *name, int def);

/**
 * ensure_directory_exists
 *
 * Create a directory structure recursively, similar to mkdir -p.
 *
 * @param dir The directory path to create
 * @return    0 on success, -1 on failure
 */
int
ensure_directory_exists(const char *dir);

/**
 * build_output_filepath
 *
 * Build the output file path based on HDF_OUTPUT_DIR environment variable.
 * Creates directory structure as needed.
 *
 * @param buf       Output buffer for the full path
 * @param sz        Size of the output buffer
 * @param file_name The file name to append (can be NULL)
 */
void
build_output_filepath(char *buf, size_t sz, const char *file_name);

#if defined(__cplusplus)
}
#endif

#endif /* OUTPUT_UTILS_H */
