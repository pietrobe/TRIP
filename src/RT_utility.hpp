#ifndef __RT_UTILITY_HPP__
#define __RT_UTILITY_HPP__

#include <petsc.h>
#include <petscsys.h>
#include <rii_emission_coefficient_3D.h>
#include <yaml-cpp/yaml.h>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "RT_types.hpp"
#include "Utilities.hpp"
#include "tools.h"

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
	int		gmres_restart	= 30;
};

struct PrecConfig
{
	KSPType pc_solver_type = "KSPGMRES";
	double	pc_rtol		   = 1e-5;
	int		pc_max_it	   = 1000;
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

	// emissivity
	emissivity_model_t emissivity_model;

	// Physical switches
	bool use_B			  = true;
	bool use_Vb 		  = true;
	bool enable_continuum = true;

	// Set constant magnetic field
	bool				  set_uniform_B = false;
	std::array<double, 3> B_field		= {0.0, 0.0, 0.0};

	// Set constant bulk velocity
	bool				  set_uniform_Vb = false;
	std::array<double, 3> Vb_field		= {0.0, 0.0, 0.0};

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

#endif
