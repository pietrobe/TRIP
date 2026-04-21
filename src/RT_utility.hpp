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
	ZERO
}; //

namespace TRIP_Comms
{

	class MPI_TRIP_Communicators : public std::enable_shared_from_this<MPI_TRIP_Communicators>
	{
	   public:
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

		inline void
		duplicateCommunicator(const std::string &name, MPI_Comm comm)
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
		addCommunicator(const std::string &name, MPI_Comm comm)
		{
			if (communicators.find(name) != communicators.end())
			{
				throw std::runtime_error("Communicator already exists: " + name);
			}
			communicators[name] = comm;
		}

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

		inline void
		setCommunicatorInfo(const std::string &name, const std::string &key, const std::string &value)
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

		inline void
		getCommunicator(const std::string &name, MPI_Comm &comm) const
		{
			auto it = communicators.find(name);
			if (it == communicators.end())
			{
				throw std::runtime_error("Communicator not found: " + name);
			}
			comm = it->second;
		}

		inline void
		getCommunicatorInfo(const std::string &name, MPI_Info &info) const
		{
			auto it = communicators_info.find(name);
			if (it == communicators_info.end())
			{
				throw std::runtime_error("Communicator info not found: " + name);
			}
			info = it->second;
		}

	   protected:
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

	inline void
	initialize_MPI_TRIP_Communicators();

	inline void
	finalize_MPI_TRIP_Communicators();

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

} // namespace YAML

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
	emissivity_model_t emissivity_model;

	// Physical switches
	bool use_B			  = true;
	bool use_Vb			  = true;
	bool enable_continuum = true;

	// Set constant magnetic field
	bool				  set_uniform_B = false;
	std::array<double, 3> B_field		= {0.0, 0.0, 0.0};

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

#endif
