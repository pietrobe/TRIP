#ifndef Utilities_hpp
#define Utilities_hpp

#include <complex>
#include <iomanip>
#include <map>
#include <numeric>
#include <ostream>
#include <algorithm>
#include <cmath>

#include "Faddeeva.hpp"
#include "Legendre_rule.hpp"
#include "RT_types.hpp"
#include "Rotation_matrix.hpp"
#include "configure.h"
#include "petsc.h"
#include "tools.h"

// Forward declarations instead of includes to break circular dependency
class RT_solver;
class RT_problem;

#ifndef MPI_CHECK
#define MPI_CHECK(stmt)                                                                                   \
	do                                                                                                    \
	{                                                                                                     \
		const int code = stmt;                                                                            \
                                                                                                          \
		if (code != MPI_SUCCESS)                                                                          \
		{                                                                                                 \
			char error_string[2048];                                                                      \
                                                                                                          \
			int length_of_error_string = sizeof(error_string);                                            \
                                                                                                          \
			MPI_Error_string(code, error_string, &length_of_error_string);                                \
                                                                                                          \
			fprintf(stderr, "ERROR!\n" #stmt " mpiAssert: %s %d %s\n", __FILE__, __LINE__, error_string); \
                                                                                                          \
			fflush(stderr);                                                                               \
                                                                                                          \
			MPI_Abort(MPI_COMM_WORLD, code);                                                              \
		}                                                                                                 \
	} while (0)
#endif

// usa std::min , std::max
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define MAX3(a, b, c) MAX(MAX(a, b), c)
#define MIN3(a, b, c) MIN(MIN(a, b), c)

#define COU_IS_ODD(n) ((n) & 1)

#define PI 3.1415926535897932384626

// for pmd input
#define ERR                                              \
	{                                                    \
		fprintf(stderr, "ERROR reading PORTA input.\n"); \
		exit(1);                                         \
	}

inline void
print_PETSc_mem(const std::string &tag = "", std::ostream &os = std::cout)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	PetscLogDouble space;
	PetscMemoryGetCurrentUsage(&space);

	// in GB
	const double mem_usage = 1e-9 * space;

	double total_mem, max_mem, min_mem;

	MPI_Reduce(&mem_usage, &min_mem, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&mem_usage, &max_mem, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&mem_usage, &total_mem, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		os << "\n===== MPI Memory Usage Summary " << (tag.empty() ? "" : "(" + tag + ")") << " =====\n";
		os << "Min memory used by single processor: " << min_mem << " GB" << std::endl;
		os << "Max memory used by single processor: " << max_mem << " GB" << std::endl;
		os << "Total memory used: " << total_mem << " GB" << std::endl;
		os << "===========================================\n";
	}
}

struct ExecutionClocks
{
	double start_time					= 0.0;
	double setup_time					= 0.0;
	double solve_end_time				= 0.0;
	double end_time						= 0.0;
	double arbitrary_beam_time			= 0.0;
	double csv_out_time					= 0.0;
	double hdf5_out_time_emergent		= 0.0;
	double hdf5_out_time_JKQ			= 0.0;
	double hdf5_out_time_whole_3D_field = 0.0;
};

inline std::string
format_execution_time_summary(const int mpi_size, const ExecutionClocks &clocks)
{
	const double total_time_hdf5 =
		clocks.hdf5_out_time_emergent + clocks.hdf5_out_time_JKQ + clocks.hdf5_out_time_whole_3D_field;

	std::stringstream ss;
	ss << std::fixed << std::setprecision(2);
	ss << std::endl;
	ss << std::string(80, '=') << std::endl;
	ss << " EXECUTION TIME SUMMARY" << std::endl;
	ss << std::string(80, '=') << std::endl;
	ss << " MPI processes:         " << mpi_size << std::endl;
	ss << " Setup time:            " << (clocks.setup_time - clocks.start_time) << " s" << std::endl;
	ss << " Solve time:            " << (clocks.solve_end_time - clocks.setup_time) << " s" << std::endl;
	ss << " Arbitrary beam time:   " << clocks.arbitrary_beam_time << " s" << std::endl;
	ss << " Post processing time:  " << (clocks.end_time - clocks.solve_end_time) << " s" << std::endl;
	ss << "     Output CSV:        " << clocks.csv_out_time << " s" << std::endl;
	ss << "     Output HDF5:       " << std::endl;
	ss << "	      Emergent field:  " << clocks.hdf5_out_time_emergent << " s" << std::endl;
	ss << "	      JKQ:             " << clocks.hdf5_out_time_JKQ << " s" << std::endl;
	ss << "	      Whole 3D field:  " << clocks.hdf5_out_time_whole_3D_field << " s" << std::endl;
	ss << "	      Total HDF5 time: " << total_time_hdf5 << " s" << std::endl;
	ss << std::string(80, '-') << std::endl;
	ss << " Total execution time:  " << (clocks.end_time - clocks.start_time) << " s" << std::endl;
	ss << std::string(80, '=') << std::endl;
	ss << std::endl;
	return ss.str();
}

inline void
print_parallel_memory_usage(const std::string &tag = "")
{
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// List of memory fields to parse
	const std::vector<std::string> fields = {"VmSize", "VmRSS", "VmData", "VmSwap"};

	std::map<std::string, long> local_values;

	// Read local /proc/self/status
	std::ifstream status("/proc/self/status");
	std::string	  line;
	while (std::getline(status, line))
	{
		for (const auto &field : fields)
		{
			if (line.find(field + ":") == 0)
			{
				std::istringstream iss(line);
				std::string		   key;
				long			   value_kb;
				std::string		   unit;
				iss >> key >> value_kb >> unit;
				local_values[field] = value_kb;
			}
		}
	}

	// Gather local values into arrays
	std::vector<long> local_array(fields.size());
	for (size_t i = 0; i < fields.size(); ++i)
	{
		local_array[i] = local_values[fields[i]];
	}

	std::vector<long> sum(fields.size()), min(fields.size()), max(fields.size());

	MPI_Reduce(local_array.data(), sum.data(), fields.size(), MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(local_array.data(), min.data(), fields.size(), MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(local_array.data(), max.data(), fields.size(), MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		std::cout << "\n===== MPI Memory Usage Summary " << (tag.empty() ? "" : "(" + tag + ")") << " =====\n";
		std::cout << "Field       |   Sum (MB)   |   Min (MB)   |   Max (MB)\n";
		std::cout << "-------------------------------------------------------\n";
		for (size_t i = 0; i < fields.size(); ++i)
		{
			double sum_mb = sum[i] / 1024.0;
			double min_mb = min[i] / 1024.0;
			double max_mb = max[i] / 1024.0;
			std::cout << std::left << std::setw(11) << fields[i] << " | " << std::setw(11) << sum_mb << " | "
					  << std::setw(11) << min_mb << " | " << std::setw(11) << max_mb << "\n";
		}
		std::cout << "=======================================================\n" << std::endl;
	}
}

// void
// compute_arbitrary_beam(RT_solver &rt_solver, std::shared_ptr<RT_problem> &rt_problem_ptr, double mu, double chi,
// 					   const std::string &output_file);

inline double *
convert_cartesian_to_spherical(const double x, const double y, const double z)
{
	static double spherical_coordinates[3];

	const double r	   = sqrt(x * x + y * y + z * z);
	const double theta = atan2(sqrt(x * x + y * y), z);
	const double chi   = atan2(y, x);

	spherical_coordinates[0] = r;
	spherical_coordinates[1] = theta;
	spherical_coordinates[2] = chi;

	return spherical_coordinates;
}

// TODO: this one is not necessary anymore with size_t -> int?
inline int
apply_periodic_bc(const int i, const int N)
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

void
save_vec(Vec &m, const char *filename, const char *name);
void
save_mat(Mat &m, const char *filename, const char *name);

void
read_vec(std::string filename, std::vector<double> &vec);

void
print_vec(std::ostream &os, const std::vector<double> &v);

PetscErrorCode
PrintVec(Vec &v);

void
print_local_sizes(const Mat &M);
void
print_global_sizes(const Mat &M);

// petsc matrix
void
create_identity_matrix(int size, Mat &Id);

// propagation matrix methods
std::vector<double>
assemble_propagation_matrix(const std::vector<double> &etas_and_rhos);
std::vector<double>
assemble_propagation_matrix(const std::vector<double> &etas, const std::vector<double> &rhos);
std::vector<double>
assemble_propagation_matrix_scaled(const std::vector<double> &etas_and_rhos);
std::vector<double>
assemble_propagation_matrix_scaled(const std::vector<double> &etas, const std::vector<double> &rhos);
void
print_propagation_matrix(const std::vector<double> &K);
void
print_Stokes(const std::vector<double> &I);
std::vector<double>
solve_4_by_4_system(const std::vector<double> &K, const std::vector<double> &rhs);
std::vector<double>
solve_4_by_4_system_optimized(const std::vector<double> &K, const std::vector<double> &rhs);

// Wigner3j symbols, use int multiples for not int inputs
double
W3JS(int J1, int J2, int J3, int M1, int M2, int M3);

// linearly interpolate vector and double its size
std::vector<double>
refine_vector(const std::vector<double> &v);
std::vector<double>
refine_vector_blocked(const std::vector<double> &v, const int block_size);
std::vector<double>
refine_vector_blocked2(const std::vector<double> &v, const int block_size_fn);

double
pow_gen(const double x, const double exponent);



namespace HeII_inelastic_collisions {

	// v values table for He II 304 from Janev (1987)
    const int N = 15;

    static const double TK[N] = {
        3000., 5000., 7000., 10000., 15000., 
		20000., 25000., 30000., 40000., 50000., 
		60000., 80000., 100000., 150000., 200000.
    };

    static const double vT[N] = {
        0.126657091063009, 0.162647119063730, 0.192237894128096, 0.229957563710787, 
		0.282637879731763, 0.327807171027031, 0.368224501041440, 0.405299561542058,
        0.472461037141359, 0.533128364313460, 0.589242272986984, 0.691950763217367,
        0.785777153064998, 0.995940776604376, 1.183297238171380
    };

	inline double b_coeff[N - 1];
    inline double c_coeff[N];
    inline double d_coeff[N - 1];
	inline bool spline_initialised = false;

    // Internal initialization function for the Spline
    void init_spline();

	/**
	 * @brief Interpolate v from table in Janev (1987)
	 * @param T temperature in K for which we want to obtain the value v
	 * @param interpolation_type type of interpolation we want to employ (0 = linear, 1 = cubic spline)
	 */
	double
	interpolate_v(const double T, const unsigned int interpolation_type = 0);

	/**
	 * @brief Compute inelastic collisions Qul (for HeII304A) in a given spatial point
	 * @param T temperature [K]
	 * @param gu degeneracy of the upper term 
	 * @note Assume that v is already the value interpolated from Janev (1987)
	 */ 
	inline double
	compute_Qul(
		const double T, 
		const int gu)
	{
		const double v = interpolate_v(T);
		return 8.63e-6 * v  / (static_cast<double>(gu) * std::sqrt(T));
	}

	/**
	 * @brief Compute inelastic collisional rate Cul=Ne*Qul [s^-1] (for HeII304A) in a given spatial point
	 * @param T temperature [K]
	 * @param gu degeneracy of the upper term 
	 * @param Ne  electron density [cm^-3]
	 * @note Assume that v is already the value interpolated from Janev (1987)
	 */ 
	inline double
	compute_Cul(
		const double T, 
		const int gu, 
		const double Ne) 
	{ 
		return Ne * compute_Qul(T, gu);
	}

}

#endif
