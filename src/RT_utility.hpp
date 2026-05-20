#ifndef __RT_UTILITY_HPP__
#define __RT_UTILITY_HPP__

#include <petsc.h>
#include <petscsys.h>
#include <rii_emission_coefficient_3D.h>
#include <yaml-cpp/yaml.h>
#include <array>
#include <filesystem>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "RT_types.hpp"
#include "Utilities.hpp"

// emissivty models
enum class emissivity_model_t
{
	NONE,		   //
	CRD_limit,	   //
	CRD_limit_VHP, //
	PRD,		   //
	PRD_NORMAL,	   //
	PRD_MEDIUM,	   //
	PRD_FAST,	   //
	PRD_AA,		   //
	PRD_AA_MAPV,   //
	ZERO		   //
}; //

inline bool is_CRD_limit(const emissivity_model_t model)
{
	return model == emissivity_model_t::CRD_limit || model == emissivity_model_t::CRD_limit_VHP;
}

enum class preconditioner_emissivity_model_t
{
	NONE,			//
	CRD_limit,		//
	PRD_AA,		//
	PRD_AA_MAPV, //
	PRD_AA_GB,		//
	PRD_AA_MAPV_GB, //
	ZERO			//
}; //

namespace TRIP_Comms
{

	/**
	 * Minimal owner/registry for MPI communicators and their MPI_Info objects.
	 *
	 * This class owns the communicators added or duplicated into it and frees them
	 * at destruction time.
	 */
	class MPI_TRIP_Communicators : public std::enable_shared_from_this<MPI_TRIP_Communicators>
	{
	   public:
		using MPI_TRIP_Communicators_shared_ptr_t = std::shared_ptr<MPI_TRIP_Communicators>;

		auto static make_shared()
		{
			return std::shared_ptr<MPI_TRIP_Communicators>(new MPI_TRIP_Communicators());
		}

		inline MPI_TRIP_Communicators() = default;

		inline ~MPI_TRIP_Communicators()
		{
			freeInfo();
			freeCommunicators();
		}

		inline MPI_TRIP_Communicators(const MPI_TRIP_Communicators &) = delete;

		inline MPI_TRIP_Communicators &
		operator=(const MPI_TRIP_Communicators &) = delete;

		inline MPI_TRIP_Communicators(MPI_TRIP_Communicators &&other) noexcept
			: communicators(std::move(other.communicators)), communicators_info(std::move(other.communicators_info))
		{
			other.communicators.clear();
			other.communicators_info.clear();
		}

		/** Duplicate a communicator and store it under the provided name. */
		inline void
		duplicateCommunicator(const std::string name, MPI_Comm comm)
		{
			auto it = communicators.find(name);
			if (it != communicators.end() && it->second != MPI_COMM_NULL)
			{
				checkMPI(MPI_Comm_free(&it->second), "MPI_Comm_free");
			}

			MPI_Comm new_comm;
			checkMPI(MPI_Comm_dup(comm, &new_comm), "MPI_Comm_dup");
			communicators[name] = new_comm;
		}

		/**
		 * Take the control of a given communicator. The communicator will be freed when the MPI_TRIP_Communicators object
		 * is destroyed.
		 */
		inline void
		addCommunicator(const std::string name, MPI_Comm comm)
		{
			if (comm == MPI_COMM_NULL)
			{
				throw std::runtime_error("MPI_COMM_NULL cannot be owned by MPI_TRIP_Communicators.");
			}

			if (isDefaultCommunicator(comm))
			{
				throw std::runtime_error(
					"Default MPI communicators (MPI_COMM_WORLD, MPI_COMM_SELF) cannot be owned by "
					"MPI_TRIP_Communicators. Use duplicateCommunicator instead."
					"MPI_COMM_NULL cannot be owned by MPI_TRIP_Communicators.");
			}

			if (communicators.find(name) != communicators.end())
			{
				throw std::runtime_error("Communicator already exists: " + name);
			}
			communicators[name] = comm;
		}

		/** Check if a communicator with the given name exists. */
		inline bool
		hasCommunicator(const std::string name) const
		{
			return communicators.find(name) != communicators.end();
		}

		/** Return true if the given MPI communicator handle is already stored. */
		inline bool
		hasCommunicator(MPI_Comm comm) const
		{
			if (comm == MPI_COMM_NULL) return false;

			for (const auto &pair : communicators)
			{
				if (pair.second == MPI_COMM_NULL) continue;

				int cmp_result = MPI_UNEQUAL;
				checkMPI(MPI_Comm_compare(pair.second, comm, &cmp_result), "MPI_Comm_compare");
				if (cmp_result == MPI_IDENT) return true;
			}

			return false;
		}

		/** Free all owned communicators and clear the communicator registry. */
		inline void
		freeCommunicators()
		{
			for (auto &pair : communicators)
			{
				if (pair.second != MPI_COMM_NULL)
				{
					checkMPI(MPI_Comm_free(&pair.second), "MPI_Comm_free");
				}
			}
			communicators.clear();
		}

		/** Free all MPI_Info objects and clear the communicator-info registry. */
		inline void
		freeInfo()
		{
			for (auto &pair : communicators_info)
			{
				if (pair.second != MPI_INFO_NULL)
				{
					checkMPI(MPI_Info_free(&pair.second), "MPI_Info_free");
				}
			}
			communicators_info.clear();
		}

		/** Set one MPI_Info key/value pair on a registered communicator. */
		inline void
		setCommunicatorInfo(const std::string name, const std::string key, const std::string value)
		{
			auto comm_it = communicators.find(name);
			if (comm_it == communicators.end())
			{
				throw std::runtime_error("Communicator not found: " + name);
			}

			auto info_it = communicators_info.find(name);
			if (info_it != communicators_info.end() && info_it->second != MPI_INFO_NULL)
			{
				checkMPI(MPI_Info_free(&info_it->second), "MPI_Info_free");
			}

			MPI_Info info;
			checkMPI(MPI_Info_create(&info), "MPI_Info_create");
			checkMPI(MPI_Info_set(info, key.c_str(), value.c_str()), "MPI_Info_set");
			checkMPI(MPI_Comm_set_info(comm_it->second, info), "MPI_Comm_set_info");

			communicators_info[name] = info;
		}

		/** Fetch a communicator by name. Throws if the name is not registered. */
		inline void
		getCommunicator(const std::string name, MPI_Comm &comm) const
		{
			auto it = communicators.find(name);
			if (it == communicators.end())
			{
				throw std::runtime_error("Communicator not found: " + name);
			}
			comm = it->second;
		}

		/** Fetch communicator info by name. Throws if no info is available. */
		inline void
		getCommunicatorInfo(const std::string name, MPI_Info &info) const
		{
			auto it = communicators_info.find(name);
			if (it == communicators_info.end())
			{
				throw std::runtime_error("Communicator info not found: " + name);
			}
			info = it->second;
		}

	   protected:
		static inline bool
		isDefaultCommunicator(MPI_Comm comm)
		{
			int cmp_result = MPI_UNEQUAL;
			checkMPI(MPI_Comm_compare(comm, MPI_COMM_WORLD, &cmp_result), "MPI_Comm_compare");
			if (cmp_result == MPI_IDENT) return true;

			checkMPI(MPI_Comm_compare(comm, MPI_COMM_NULL, &cmp_result), "MPI_Comm_compare");
			if (cmp_result == MPI_IDENT) return true;

			checkMPI(MPI_Comm_compare(comm, MPI_COMM_SELF, &cmp_result), "MPI_Comm_compare");
			return cmp_result == MPI_IDENT;
		}

		static inline void
		checkMPI(const int err, const char *where)
		{
			if (err == MPI_SUCCESS) return;

			char errstr[MPI_MAX_ERROR_STRING];
			int	 len = 0;
			MPI_Error_string(err, errstr, &len);
			throw std::runtime_error(std::string(where) + " failed: " + std::string(errstr, len));
		}

		std::map<std::string, MPI_Comm> communicators;
		std::map<std::string, MPI_Info> communicators_info;
	};

	void
	initialize_MPI_TRIP_Communicators();

	void
	finalize_MPI_TRIP_Communicators();

	MPI_TRIP_Communicators::MPI_TRIP_Communicators_shared_ptr_t
	getTRIPCommunicators();

} // namespace TRIP_Comms

namespace YAML
{
	template <>
	struct convert<emissivity_model_t>
	{
		static Node
		encode(const emissivity_model_t &rhs)
		{
			Node node;
			switch (rhs)
			{
				case emissivity_model_t::NONE:
					node = "NONE";
					break;
				case emissivity_model_t::CRD_limit:
					node = "CRD_limit";
					break;
				case emissivity_model_t::CRD_limit_VHP:
					node = "CRD_limit_VHP";
					break;
				case emissivity_model_t::PRD:
					node = "PRD";
					break;
				case emissivity_model_t::PRD_NORMAL:
					node = "PRD_NORMAL";
					break;
				case emissivity_model_t::PRD_MEDIUM:
					node = "PRD_MEDIUM";
					break;
				case emissivity_model_t::PRD_FAST:
					node = "PRD_FAST";
					break;
				case emissivity_model_t::PRD_AA:
					node = "PRD_AA";
					break;
				case emissivity_model_t::PRD_AA_MAPV:
					node = "PRD_AA_MAPV";
					break;
				case emissivity_model_t::ZERO:
					node = "ZERO";
					break;
			}
			return node;
		}

		static bool
		decode(const Node &node, emissivity_model_t &rhs)
		{
			if (!node.IsScalar()) return false;
			const std::string s = node.as<std::string>();

			if (s == "NONE")
				rhs = emissivity_model_t::NONE;
			else if (s == "CRD_limit")
				rhs = emissivity_model_t::CRD_limit;
			else if (s == "CRD_limit_VHP")
				rhs = emissivity_model_t::CRD_limit_VHP;
			else if (s == "PRD")
				rhs = emissivity_model_t::PRD;
			else if (s == "PRD_NORMAL")
				rhs = emissivity_model_t::PRD_NORMAL;
			else if (s == "PRD_MEDIUM")
				rhs = emissivity_model_t::PRD_MEDIUM;
			else if (s == "PRD_FAST")
				rhs = emissivity_model_t::PRD_FAST;
			else if (s == "PRD_AA")
				rhs = emissivity_model_t::PRD_AA;
			else if (s == "PRD_AA_MAPV")
				rhs = emissivity_model_t::PRD_AA_MAPV;
			else if (s == "ZERO")
				rhs = emissivity_model_t::ZERO;
			else
				return false;

			return true;
		}
	};

	template <>
	struct convert<preconditioner_emissivity_model_t>
	{
		static Node
		encode(const preconditioner_emissivity_model_t &rhs)
		{
			Node node;
			switch (rhs)
			{
				case preconditioner_emissivity_model_t::NONE:
					node = "NONE";
					break;
				case preconditioner_emissivity_model_t::CRD_limit:
					node = "CRD_limit";
					break;
				case preconditioner_emissivity_model_t::PRD_AA:
					node = "PRD_AA";
					break;
				case preconditioner_emissivity_model_t::PRD_AA_MAPV:
					node = "PRD_AA_MAPV";
					break;
				case preconditioner_emissivity_model_t::PRD_AA_GB:
					node = "PRD_AA_GB";
					break;
				case preconditioner_emissivity_model_t::PRD_AA_MAPV_GB:
					node = "PRD_AA_MAPV_GB";
					break;
				case preconditioner_emissivity_model_t::ZERO:
					node = "ZERO";
					break;
			}
			return node;
		}

		static bool
		decode(const Node &node, preconditioner_emissivity_model_t &rhs)
		{
			if (!node.IsScalar()) return false;
			const std::string s = node.as<std::string>();

			if (s == "NONE")
				rhs = preconditioner_emissivity_model_t::NONE;
			else if (s == "CRD_limit")
				rhs = preconditioner_emissivity_model_t::CRD_limit;
			else if (s == "PRD_AA")
				rhs = preconditioner_emissivity_model_t::PRD_AA;
			else if (s == "PRD_AA_MAPV")
				rhs = preconditioner_emissivity_model_t::PRD_AA_MAPV;
			else if (s == "PRD_AA_GB")
				rhs = preconditioner_emissivity_model_t::PRD_AA_GB;
			else if (s == "PRD_AA_MAPV_GB")
				rhs = preconditioner_emissivity_model_t::PRD_AA_MAPV_GB;
			else if (s == "ZERO")
				rhs = preconditioner_emissivity_model_t::ZERO;
			else
				return false;

			return true;
		}
	};

} // namespace YAML

template <typename T>
T
requiredField(const YAML::Node &config, const std::string &key)
{
	if (!config[key]) throw std::runtime_error("Missing required key: " + key);
	return config[key].as<T>();
}

// function to read vector or scalar
inline std::vector<double>
readDoubleVec(const YAML::Node &node)
{
	if (node.IsSequence())
	{
		return node.as<std::vector<double>>();
	}
	else if (node.IsScalar())
	{
		return {node.as<double>()}; // wrap scalar in a vector
	}
	else
	{
		throw std::runtime_error("Expected a scalar or sequence");
	}
}

// Function declarations
std::string
getCurrentDateTime();

std::map<std::string, std::string>
get_input_PORTA(const std::filesystem::path &config_file, int mpi_rank);

std::string
getOptionArgument(int argc, char *argv[], const std::string &option, const std::string &defaultValue = "");

bool
getOptionFlag(int argc, char *argv[], const std::string &option);

void
print_help();

// Structure for arbitrary beam direction
struct BeamDirection
{
	double mu  = 1.0;
	double chi = 0.0;
};

// structures for input file
struct SolverConfig
{
	KSPType ksp_solver_type = "KSPFGMRES";
	double	ksp_rtol		= 1e-5;
	int		ksp_max_it		= 1000;
	bool	ksp_use_J_KQ	= false;
	int		gmres_restart	= 30;
};

struct PrecConfig
{
	KSPType pc_solver_type = "KSPGMRES";
	double	pc_rtol		   = 1e-5;
	int		pc_max_it	   = 1000;
	bool	pc_use_J_KQ	   = false;
	bool	verbose	       = false;
};

struct AtomConfig
{
	double mass = 40.078;
	double Aul	= 2.18e+08;

	int S2 = 0;

	int Ll2 = 0;
	int Lu2 = 2;

	// vectors in case of two-terms atoms
	std::vector<double> El_vec = {0.0};
	std::vector<double> Eu_vec = {23652.304};

	std::vector<double> gl_vec = {0.0};
	std::vector<double> gu_vec = {1.0};
};

struct AppConfig
{
	// Main I/O
	std::filesystem::path input_directory;
	std::filesystem::path input_file;
	std::filesystem::path frequency_file;
	std::filesystem::path output_directory;

	// Output settings
	bool output						 = false;
	bool output_overwrite_prevention = false;
	bool write_whole_3D_field_hdf5	 = false;
	bool write_text_output			 = false;

	// testing
	std::filesystem::path reference_sol_directory;

	// emissivity
	emissivity_model_t				  emissivity_model;
	preconditioner_emissivity_model_t preconditioner_emissivity_model;

	// Physical switches
	bool use_B			   = true;
	bool use_Vb			   = true;
	bool enable_continuum  = true;
	bool use_D2_from_input = true; // if false, D2 is calculated from the a and b coefficients.

	// Set constant magnetic field
	bool				  set_uniform_B = false;
	std::array<double, 3> B_field		= {0.0, 0.0, 0.0};

	bool enable_sigma_c = true; // It could be disabled for comparison with PORTA, which does not include sigma_c in the
								// emissivity calculation.

	// Set constant bulk velocity
	bool				  set_uniform_Vb = false;
	std::array<double, 3> Vb_field		 = {0.0, 0.0, 0.0};

	// numerical inputs
	bool		use_1_5D_approx = false;
	std::string formal_solver	= "BESSER";

	// Angular grid sizes
	int N_theta = 8;
	int N_chi	= 16;

	// Horizontal grid sizes
	int	   N_x = 1;
	int	   N_y = 1;
	double L   = 400.0;

	// Optional input strings
	std::filesystem::path input_cul	 = "";
	std::filesystem::path input_qel	 = "";
	std::filesystem::path input_llp	 = "";
	std::filesystem::path input_back = "";

	// Use preconditioner
	bool use_prec = true;

	// Subsections
	SolverConfig solver;
	PrecConfig	 prec;
	AtomConfig	 atom;

	// Arbitrary beam directions
	std::vector<BeamDirection> arbitrary_beams;

	// print out
	bool verbose = true;
};

AppConfig
loadConfig(const std::string &filename);

void
writeConfigResume(const AppConfig &cfg, std::ostream &os);

// std::shared_ptr<RT_problem>
// create_rt_problem(const AppConfig &cfg, const std::filesystem::path &input_file_path,
// 				  const std::filesystem::path &frequencies_input_path, emissivity_model emissivity_model_var,
// 				  int mpi_rank);

int
acc_devices_print_info(const int mpi_rank, const int mpi_size, std::ostream &os);

int
write_emergent_field_hdf5(RT_problem &rt_problem_ptr, const std::string &output_file);

int
write_3D_whole_field_hdf5(RT_problem &rt_problem, const std::string &output_file, //
						  const int step_z_ = 2, bool normalized_output = false); //

int
write_3D_whole_field_falp_hdf5(RT_problem &rt_problem, const std::string &output_file, //
							   const int step_z_ = 2, bool normalized_output = false); //

/**
 * @brief Compute the collisional depolarization rate D2.
 *
 * Reference: Collisional Depolarization of the Solar Ca, Mg, and Ba Levels,
 * M. Derouich (2019), DOI: 10.3847/1538-4357/ab26b4,
 * From Eq. (6).
 *
 * In the paper, T_ref is set to 5000 K with a = 1.24e-8 and b = 0.37.
 * In this implementation, a and b are provided in the h5 input file.
 * D2 is calculated as D2 = a * nH * (T / T_ref)^b, where nH is the neutral hydrogen density and T is the local
 * temperature.
 *
 * @param a Empirical coefficient from input data.
 * @param b Empirical exponent from input data.
 * @param nH Neutral hydrogen density.
 * @param T Local temperature.
 * @param T_ref Reference temperature (default: 5000.0 K).
 * @return Collisional depolarization coefficient D2.
 */
double
calculate_D2(const double a, const double b, const double nH, const double T, const double T_ref = 5000.0); //

#endif
