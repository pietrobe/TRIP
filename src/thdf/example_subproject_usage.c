/*
 * Example: Using hdf_output library as a dependency
 * 
 * This demonstrates how to use the hdf_output static library
 * in your own CMake project.
 */

#include <stdio.h>
#include <hdf5.h>
#include "output_utils.h"

int main(void) {
    printf("HDF5 Output Library Example\n");
    printf("HDF5 Version: %d.%d.%d\n\n", H5_VERS_MAJOR, H5_VERS_MINOR, H5_VERS_RELEASE);
    
    printf("Successfully linked against hdf_output library!\n");
    printf("Available functions from output_utils:\n");
    printf("  - All HDF5 output utilities\n");
    printf("  - Field writing functions\n");
    printf("  - Array field support\n");
    printf("  - Z-slice functionality\n");
    printf("  - Domain decomposition support\n");
    
    return 0;
}
