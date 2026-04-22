#ifndef RT_types_hpp
#define RT_types_hpp

#include <petscsys.h>
#include <vector>

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

// some methods used for type casting Real <--> double
template<typename T_in, typename T_out>
inline void convert_vector(const std::vector<T_in>& x, std::vector<T_out>& y)
{
    y.resize(x.size());
    for (size_t i = 0; i < x.size(); ++i)
        y[i] = static_cast<T_out>(x[i]);
}

// specialization for same type (fast path)
template<typename T>
inline void convert_vector(const std::vector<T>& x, std::vector<T>& y)
{
    y = x;  // uses optimized copy
}

// float → double 
inline std::vector<double> to_double(const float* x, int n)
{
    std::vector<double> y(n);
    for (int i = 0; i < n; ++i)
        y[i] = static_cast<double>(x[i]);
    return y;
}

// double → double (direct copy)
inline std::vector<double> to_double(const double* x, int n)
{
    return std::vector<double>(x, x + n);
}

#endif // RT_types_hpp