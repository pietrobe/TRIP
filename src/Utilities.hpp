#ifndef Utilities_hpp
#define Utilities_hpp

#include <complex>
#include <numeric>
#include "petsc.h" 
#include "Legendre_rule.hpp"
#include "Faddeeva.hpp"
#include "Rotation_matrix.hpp"

using Real = double; 

#ifndef MPI_CHECK
#define MPI_CHECK(stmt)                         \
    do                                  \
    {                                 \
    const int code = stmt;                      \
                                    \
    if (code != MPI_SUCCESS)                    \
    {                               \
      char error_string[2048];                  \
                                    \
      int length_of_error_string = sizeof(error_string);      \
                                    \
      MPI_Error_string                      \
        (code, error_string, &length_of_error_string);      \
                                    \
      fprintf(stderr,                       \
          "ERROR!\n" #stmt " mpiAssert: %s %d %s\n",      \
          __FILE__, __LINE__, error_string);          \
                                    \
      fflush(stderr);                       \
                                    \
      MPI_Abort(MPI_COMM_WORLD, code);              \
    }                               \
    }                                 \
    while(0)
#endif

// usa std::min , std::max
#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))

#define MAX3(a, b, c) MAX(MAX(a, b), c)
#define MIN3(a, b, c) MIN(MIN(a, b), c)

#define COU_IS_ODD(n) ((n)&1)

#define PI 3.1415926535897932384626

// for pmd input
#define ERR {fprintf(stderr,"ERROR reading PORTA input.\n"); exit(1);}

inline double* convert_cartesian_to_spherical(const double x, const double y, const double z)
{
    static double spherical_coordinates[3]; 

    const double r = sqrt(x*x + y*y + z*z);   
    const double theta = atan2(sqrt(x * x + y * y), z);
    const double chi   = atan2(y, x);   
    
    spherical_coordinates[0] = r;
    spherical_coordinates[1] = theta;
    spherical_coordinates[2] = chi;

    return spherical_coordinates;
}

inline int apply_periodic_bc(const int i, const size_t N)
{
    const int i_new = (i < 0 ) ? i + N : i % N;

    if (i_new < 0 ) std::cout << "WARNING: negative index in apply BC";                 
                
    return i_new;
}


void save_vec(Vec &m, const char * filename, const char * name);
void save_mat(Mat &m, const char * filename, const char * name);

void read_vec(std::string filename, std::vector<double> &vec);

void print_vec(const std::vector<double> &vec);

void print_local_sizes(const Mat &M);
void print_global_sizes(const Mat &M);

// petsc matrix
void create_identity_matrix(int size, Mat &Id);

// propagation matrix methods
std::vector<double> assemble_propagation_matrix(const std::vector<double> &etas_and_rhos);
std::vector<double> assemble_propagation_matrix(const std::vector<double> &etas, const std::vector<double> &rhos); 
std::vector<double> assemble_propagation_matrix_scaled(const std::vector<double> &etas_and_rhos); 
std::vector<double> assemble_propagation_matrix_scaled(const std::vector<double> &etas, const std::vector<double> &rhos); 
void print_propagation_matrix(const std::vector<double> &K);
void print_Stokes(const std::vector<double> &I);
std::vector<double> solve_4_by_4_system(const std::vector<double> &K, const std::vector<double> &rhs);
std::vector<double> solve_4_by_4_system_optimized(const std::vector<double> &K, const std::vector<double> &rhs);

// Wigner3j symbols, use int multiples for not int inputs 
double W3JS(int J1, int J2, int J3, int M1, int M2, int M3);

// linearly interpolate vector and double its size 
std::vector<double> refine_vector(const std::vector<double> &v);
std::vector<double> refine_vector_blocked(const std::vector<double> &v, const size_t block_size);
std::vector<double> refine_vector_blocked2(const std::vector<double> &v, const size_t block_size_fn);

double pow_gen(const double x, const double exponent);

#endif 
