#include "tools.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <regex>
#include <sstream>
#include <string>

// std::string get_arg(const std::string &input, const std::string &word) {
//   std::regex pattern("^\\s+" + word + ":\\s+(.*)$", std::regex::multiline);
//   std::smatch match;

//   if (std::regex_search(input, match, pattern)) {
//     return match.str(1);
//   }

//   return std::string();
// }

std::string get_arg(const std::string &input, const std::string &word) {
    std::regex pattern("^\\s*" + word + ":\\s*(.*)$");
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

void print_help() {
  std::cout
      << "----------------------------------------------------------------"
      << std::endl
      << std::endl;
  std::cout << "Usage: mpirun ./solar_3D [options] [PETSC options]"
            << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  --CRD: use CRD" << std::endl;
  std::cout << "  --input_dir <input_dir>: input directory" << std::endl;
  std::cout << "  --output_dir <output_dir>: output directory" << std::endl;
  std::cout << "  --problem_pmd_file <problem_pmd_file>: problem pmd input file"
            << std::endl;
  std::cout << "  --problem_input_config <problem_input_config>: problem input "
               "config file";
  std::cout << "  --help: print help and exit" << std::endl << std::endl;

  std::cout
      << "----------------------------------------------------------------"
      << std::endl
      << std::endl;

  std::cout << "Example: mpirun ./solar_3D --CRD --input_dir /path/to/input "
               "--output_dir /path/to/output --problem_input_file "
               "problem_input_file.pmd  -ksp_type gmres -ksp_max_it 100 "
               "-ksp_rtol 1e-10"
            << std::endl;
  std::cout << std::endl;
  std::cout << "In the output directory, the code creates a results directory "
               "with the name of the problem input file and the extension "
               "\".CRD\" or \".PRD\", depending on the --CRD option."
            << std::endl;
  std::cout
      << "If the results output directory already exists, the code will stop."
      << std::endl;
  std::cout << "Default solver is the PRD." << std::endl;
  std::cout << std::endl;
  std::cout << "The config file is a text file in the input directory with the "
               "following format:"
            << std::endl
            << std::endl;
  std::cout << "  pmd: problem.pmd" << std::endl;
  std::cout << "  cul: cul_file.cul" << std::endl;
  std::cout << "  qel: qel_file.qel" << std::endl;
  std::cout << "  llp: llp_file.llp" << std::endl;
  std::cout << std::endl;
}
