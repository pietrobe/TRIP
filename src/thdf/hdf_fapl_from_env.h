#ifndef __HDF_FAPL_FROM_ENV_H__
#define __HDF_FAPL_FROM_ENV_H__

#include <hdf5.h>
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

hid_t
build_fapl_from_env(MPI_Comm comm);

#ifdef __cplusplus
}
#endif

#endif  // __HDF_FAPL_FROM_ENV_H__