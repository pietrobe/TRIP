#ifndef RT_types_hpp
#define RT_types_hpp

#include <petscsys.h>

typedef PetscScalar Real; // TODO: remove and use simply PetscScalar everywhere?

template<class T>
MPI_Datatype mpi_type() {
        if (std::is_same<T, double>::value)   return MPI_DOUBLE;
        if (std::is_same<T, float>::value)    return MPI_FLOAT;
        if (std::is_same<T, int>::value)      return MPI_INT;
        if (std::is_same<T, long>::value)     return MPI_LONG;
        if (std::is_same<T, PetscInt>::value) return MPIU_INT;
        return MPI_BYTE; 
    }

#endif // RT_types_hpp