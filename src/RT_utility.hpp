#ifndef __RT_UTILITY_HPP__
#define __RT_UTILITY_HPP__

#include <string>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <filesystem>
#include <vector>
#include "RT_problem.hpp"

// Function declarations
std::string getCurrentDateTime();

std::map<std::string, std::string> get_input_PORTA(const std::filesystem::path &config_file, int mpi_rank);

std::string getOptionArgument(int argc, char *argv[],
                              const std::string &option,
                              const std::string &defaultValue = "");

bool getOptionFlag(int argc, char *argv[], const std::string &option);

void print_help();

// Structure for arbitrary beam direction
struct BeamDirection {
    double mu  = 1.0;
    double chi = 0.0;
};

// structures for input file
struct SolverConfig {
    KSPType ksp_solver_type = "KSPFGMRES";
    double      ksp_rtol    = 1e-5;
    int         ksp_max_it  = 1000;
};

struct PrecConfig {
    KSPType pc_solver_type = "KSPGMRES";
    double      pc_rtol    = 1e-5;
    int         pc_max_it  = 1000;
};

struct AppConfig {
    // Main I/O
    std::filesystem::path input_directory;
    std::filesystem::path input_file;
    std::filesystem::path frequency_file;
    std::filesystem::path output_directory;

    // Output settings
    bool output                      = false;
    bool output_overwrite_prevention = false;

    // emissivity
    emissivity_model emissivity_model_o;
    
    // Physical switches
    bool use_B          = true;
    bool use_continuum  = false;
    std::string formal_solver = "BESSER";

    // Grid sizes
    int N_theta = 8;
    int N_chi   = 16;

    // Optional input strings
    std::filesystem::path input_cul  = "";
    std::filesystem::path input_qel  = "";
    std::filesystem::path input_llp  = "";
    std::filesystem::path input_back = "";

    // Use preconditioner
    bool use_prec = true;

    // Subsections
    SolverConfig solver;
    PrecConfig   prec;

    // Arbitrary beam directions
    std::vector<BeamDirection> arbitrary_beams;
};

AppConfig loadConfig(const std::string& filename);


std::shared_ptr<RT_problem>
create_rt_problem(const AppConfig &cfg, const std::filesystem::path &input_file_path,
				  const std::filesystem::path &frequencies_input_path, emissivity_model emissivity_model_var,
				  int mpi_rank);

#endif
