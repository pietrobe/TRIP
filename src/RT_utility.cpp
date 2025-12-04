#include "RT_utility.hpp"
#include "tools.h"
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>

inline KSPType toKSPType(const std::string& name)
{
    if (name == "KSPGMRES")  return KSPGMRES;
    if (name == "KSPFGMRES") return KSPFGMRES;  // PETSc literal
    if (name == "KSPBCGS")   return KSPBCGS;
    if (name == "KSPPREONLY") return KSPPREONLY;

    throw std::runtime_error("Unknown KSPType: " + name);
}

inline std::string validateFormalSolver(const std::string& s)
{
    static const std::vector<std::string> allowed{
        "implicit_Euler",
        "trapezoidal",
        "Crank–Nicolson",
        "DELO_linear",
        "BESSER"
    };

    for (const auto& a : allowed)
        if (s == a) return s;

    throw std::runtime_error("Invalid formal_solver: " + s);
}


AppConfig loadConfig(const std::string& filename) {
    YAML::Node config = YAML::LoadFile(filename);
    AppConfig cfg;

    // Required strings
    auto requiredString = [&](const std::string& key) {
        if (!config[key])
            throw std::runtime_error("Missing required key: " + key);
        return config[key].as<std::string>();
    };

    cfg.input_directory = std::filesystem::path(requiredString("input_directory"));
    cfg.input_file      = std::filesystem::path(requiredString("input_file"));
    cfg.frequency_file  = std::filesystem::path(requiredString("frequency_file"));

    // Optional booleans with defaults
    if (config["output"]) cfg.output = config["output"].as<bool>();

    if (config["output_overwrite_prevention"]) cfg.output_overwrite_prevention = config["output_overwrite_prevention"].as<bool>();

    // Optional string (converted to filesystem::path)
    if (config["output_directory"]) cfg.output_directory = std::filesystem::path(config["output_directory"].as<std::string>());

    // Emissivity model (required)
    cfg.emissivity_model = config["emissivity_model"].as<emissivity_model>();

    // Physical flags
    if (config["use_B"]) cfg.use_B = config["use_B"].as<bool>();    

    // Formal solver
    if (config["formal_solver"]) 
    {
        std::string fs = config["formal_solver"].as<std::string>();
        cfg.formal_solver = validateFormalSolver(fs);
    }

    // Integers
    if (config["N_theta"]) cfg.N_theta = config["N_theta"].as<int>();

    if (config["N_chi"]) cfg.N_chi = config["N_chi"].as<int>();

    // Optional strings (converted to filesystem::path)
    if (config["input_cul"])  cfg.input_cul  = std::filesystem::path(config["input_cul"].as<std::string>());
    if (config["input_qel"])  cfg.input_qel  = std::filesystem::path(config["input_qel"].as<std::string>());
    if (config["input_llp"])  cfg.input_llp  = std::filesystem::path(config["input_llp"].as<std::string>());
    if (config["input_back"]) cfg.input_back = std::filesystem::path(config["input_back"].as<std::string>());

    // use_prec (default is true)
    if (config["use_prec"]) cfg.use_prec = config["use_prec"].as<bool>();

    // override logic: do not use preconditioner for CRD and ZERO
    switch (cfg.emissivity_model) {
        case emissivity_model::CRD_limit:
        case emissivity_model::CRD_limit_VHP:
        case emissivity_model::ZERO:
            cfg.use_prec = false;
            break;
        default:
            break;
    }

    // Solver section
    if (config["solver"]) 
    {
        auto s = config["solver"];   

        if (s["ksp_solver_type"]) cfg.solver.ksp_solver_type = toKSPType(s["ksp_solver_type"].as<std::string>());        
        if (s["ksp_rtol"])        cfg.solver.ksp_rtol        = s["ksp_rtol"].as<double>();
        if (s["ksp_max_it"])      cfg.solver.ksp_max_it      = s["ksp_max_it"].as<int>();
    }

    // Preconditioner section
    if (config["prec"]) 
    {
        auto p = config["prec"];

        if (p["pc_solver_type"]) cfg.prec.pc_solver_type = toKSPType(p["pc_solver_type"].as<std::string>());
        if (p["pc_rtol"])        cfg.prec.pc_rtol        = p["pc_rtol"].as<double>();
        if (p["pc_max_it"])      cfg.prec.pc_max_it      = p["pc_max_it"].as<int>();
    }

    return cfg;
}


// std::string get_arg(const std::string &input, const std::string &word) {
//   std::regex pattern("^\\s+" + word + ":\\s+(.*)$", std::regex::multiline);
//   std::smatch match;

//   if (std::regex_search(input, match, pattern)) {
//     return match.str(1);
//   }

//   return std::string();
// }

std::string get_arg(const std::string &input, const std::string &word) {
  std::regex pattern("^\\s*" + word + ":\\s*([^\\s]*)\\s*$");
  std::smatch match;
  std::istringstream stream(input);
  std::string line;

  while (std::getline(stream, line)) {
    if (std::regex_search(line, match, pattern)) {
      return match.str(1);
    }
  }

  return std::string();
}

std::map<std::string, std::string>
get_input_PORTA(const std::filesystem::path &config_file, int mpi_rank) {

  std::map<std::string, std::string> input_PORTA_map;

  std::ifstream input_file(config_file.string());

  if (input_file) {

    // read the input file
    std::stringstream buffer;
    buffer << input_file.rdbuf();
    std::string input_string = buffer.str();

    std::string pmd_file = get_arg(input_string, "pmd");
    std::string cul_file = get_arg(input_string, "cul");
    std::string qel_file = get_arg(input_string, "qel");
    std::string llp_file = get_arg(input_string, "llp");
    std::string back_file = get_arg(input_string, "back");

    if (pmd_file.empty()) {
      if (mpi_rank == 0)
        std::cerr << "Error reading PORTA input file: " << config_file
                  << std::endl;
      exit(1);
    }

    input_PORTA_map["pmd"] = pmd_file;
    input_PORTA_map["cul"] = cul_file;
    input_PORTA_map["qel"] = qel_file;
    input_PORTA_map["llp"] = llp_file;
    input_PORTA_map["back"] = back_file;

  } else {

    if (mpi_rank == 0)
      std::cerr << "Error in opening PORTA input file: " << config_file
                << std::endl;
    exit(1);
  }

  input_file.close();
  return input_PORTA_map;
}

std::string getCurrentDateTime() {
  auto now = std::chrono::system_clock::now();
  auto in_time_t = std::chrono::system_clock::to_time_t(now);

  std::stringstream ss;
  ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %H:%M:%S");
  return ss.str();
}

// compile and run with:
// make -j 32 && srun -n 4 ./main

std::string getOptionArgument(int argc, char *argv[],
                              const std::string &option) {
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == option && i + 1 < argc) {
      return argv[i + 1];
    }
  }
  return std::string(); // Option not found or argument missing
}

bool getOptionFlag(int argc, char *argv[], const std::string &option) {
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == option) {
      return true;
    }
  }
  return false;
}

// TODO
void print_help() {

  std::string help_string = R"(
----------------------------------------------------------------

Usage: mpirun ./solar_3D [PETSC options]

In the output directory, the code creates a results directory with the name of the problem input 
file and the extension ".CRD" or ".PRD", depending on the --CRD or the --epsilon_line options.
If the results output directory already exists, the code will stop.
Default solver is the PRD.
)";

  std::cout << help_string << std::endl;
}
