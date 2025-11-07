#ifndef Utilities_hpp
#define Utilities_hpp

#include <complex>
#include <numeric>
#include <map>
#include <ostream>
#include <iomanip>
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

inline void print_PETSc_mem(const std::string& tag = "")
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    PetscLogDouble space;
    PetscMemoryGetCurrentUsage(&space);   

    // in GB
    const double mem_usage = 1e-9 * space;
    
    double total_mem, max_mem, min_mem;
    
    MPI_Reduce(&mem_usage, &min_mem,   1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mem_usage, &max_mem,   1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mem_usage, &total_mem, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            
    if (rank == 0)
    {
        std::cout << "\n===== MPI Memory Usage Summary " << (tag.empty() ? "" : "(" + tag + ")") << " =====\n";
        std::cout << "Min memory used by single processor: " << min_mem   << " GB"<< std::endl;
        std::cout << "Max memory used by single processor: " << max_mem   << " GB"<< std::endl;
        std::cout << "Total memory used: "                   << total_mem << " GB"<< std::endl;
        std::cout << "===========================================\n";
    }
}

inline void print_parallel_memory_usage(const std::string& tag = "") 
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // List of memory fields to parse
    const std::vector<std::string> fields = {"VmSize", "VmRSS", "VmData", "VmSwap"};

    std::map<std::string, long> local_values;

    // Read local /proc/self/status
    std::ifstream status("/proc/self/status");
    std::string line;
    while (std::getline(status, line)) {
        for (const auto& field : fields) {
            if (line.find(field + ":") == 0) {
                std::istringstream iss(line);
                std::string key;
                long value_kb;
                std::string unit;
                iss >> key >> value_kb >> unit;
                local_values[field] = value_kb;
            }
        }
    }

    // Gather local values into arrays
    std::vector<long> local_array(fields.size());
    for (size_t i = 0; i < fields.size(); ++i) {
        local_array[i] = local_values[fields[i]];
    }

    std::vector<long> sum(fields.size()), min(fields.size()), max(fields.size());

    MPI_Reduce(local_array.data(), sum.data(), fields.size(), MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_array.data(), min.data(), fields.size(), MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_array.data(), max.data(), fields.size(), MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "\n===== MPI Memory Usage Summary " << (tag.empty() ? "" : "(" + tag + ")") << " =====\n";
        std::cout << "Field       |   Sum (MB)   |   Min (MB)   |   Max (MB)\n";
        std::cout << "-------------------------------------------------------\n";
        for (size_t i = 0; i < fields.size(); ++i) {
            double sum_mb = sum[i] / 1024.0;
            double min_mb = min[i] / 1024.0;
            double max_mb = max[i] / 1024.0;
            std::cout << std::left << std::setw(11) << fields[i] << " | "
                      << std::setw(11) << sum_mb << " | "
                      << std::setw(11) << min_mb << " | "
                      << std::setw(11) << max_mb << "\n";
        }
        std::cout << "=======================================================\n" << std::endl;
    }
}


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

// TODO: this one is not necessary anymore with size_t -> int?
inline int apply_periodic_bc(const int i, const int N)
{    
    int i_new;

    if (i > 0)
    {
        i_new = i % N;
    }
    else
    {
        i_new = i;

        while (i_new < 0) i_new += N;
    }

    return i_new;
}


void save_vec(Vec &m, const char * filename, const char * name);
void save_mat(Mat &m, const char * filename, const char * name);

void read_vec(std::string filename, std::vector<double> &vec);

void print_vec(std::ostream& os, const std::vector<double>& v);

PetscErrorCode PrintVec(Vec &v);

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
std::vector<double> refine_vector_blocked(const std::vector<double> &v, const int block_size);
std::vector<double> refine_vector_blocked2(const std::vector<double> &v, const int block_size_fn);

double pow_gen(const double x, const double exponent);

#endif 
