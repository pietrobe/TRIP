#include "RT_utility.hpp"
#include <iomanip>
#include "RT_problem.hpp"

inline KSPType
toKSPType(const std::string &name)
{
	if (name == "KSPGMRES") return KSPGMRES;
	if (name == "KSPFGMRES") return KSPFGMRES; // PETSc literal
	if (name == "KSPPGMRES") return KSPPGMRES;
	if (name == "KSPBCGS") return KSPBCGS;
	if (name == "KSPPREONLY") return KSPPREONLY;

	throw std::runtime_error("Unknown KSPType: " + name);
}

inline std::string
KSPTypeToString(KSPType type)
{
	if (type == KSPGMRES) return "KSPGMRES";
	if (type == KSPFGMRES) return "KSPFGMRES";
	if (type == KSPPGMRES) return "KSPPGMRES";
	if (type == KSPBCGS) return "KSPBCGS";
	if (type == KSPPREONLY) return "KSPPREONLY";

	return "UNKNOWN";
}

inline std::string
validateFormalSolver(const std::string &s)
{
	static const std::vector<std::string> allowed{"implicit_Euler", "trapezoidal", "Crank–Nicolson", "DELO_linear",
												  "BESSER"};

	for (const auto &a : allowed)
		if (s == a) return s;

	throw std::runtime_error("Invalid formal_solver: " + s);
}

// std::shared_ptr<RT_problem>
// create_rt_problem(const AppConfig &cfg, const std::filesystem::path &input_file_path,
// 				  const std::filesystem::path &frequencies_input_path, emissivity_model emissivity_model_var,
// 				  int mpi_rank)
// {
// 	if (cfg.input_directory.string().find("FAL-C") != std::string::npos)
// 	{
// 		if (mpi_rank == 0) std::cout << "Using FAL-C input file:  " << input_file_path << std::endl;

// 		return std::make_shared<RT_problem>(input_file_path.string(), cfg.N_theta, cfg.N_chi, emissivity_model_var,
// 											cfg.use_B);
// 	}
// 	else if (cfg.input_cul.string().empty() || cfg.input_qel.empty() || cfg.input_llp.empty())
// 	{
// 		if (mpi_rank == 0) std::cout << "Using PORTA PMD input file ONLY:  " << input_file_path << std::endl;

// 		return std::make_shared<RT_problem>(input_file_path.string().c_str(), frequencies_input_path.string().c_str(),
// 											emissivity_model_var, cfg.use_B);
// 	}
// 	else
// 	{
// 		if (mpi_rank == 0) std::cout << "Using PORTA PMD + CUL + QEL + LLP + BACK input files" << std::endl;

// 		auto input_cul_path	 = cfg.input_directory / cfg.input_cul;
// 		auto input_qel_path	 = cfg.input_directory / cfg.input_qel;
// 		auto input_llp_path	 = cfg.input_directory / cfg.input_llp;
// 		auto input_back_path = cfg.input_directory / cfg.input_back;

// 		return std::make_shared<RT_problem>(input_file_path.string().c_str(), input_cul_path.string().c_str(),
// 											input_qel_path.string().c_str(), input_llp_path.string().c_str(),
// 											input_back_path.string().c_str(), frequencies_input_path.string(),
// 											emissivity_model_var, cfg.use_B);
// 	}
// }

AppConfig
loadConfig(const std::string &filename)
{
	YAML::Node config = YAML::LoadFile(filename);
	AppConfig  cfg;

	// Required strings
	auto requiredString = [&](const std::string &key)
	{
		if (!config[key]) throw std::runtime_error("Missing required key: " + key);
		return config[key].as<std::string>();
	};

	cfg.input_directory = std::filesystem::path(requiredString("input_directory"));
	cfg.input_file		= std::filesystem::path(requiredString("input_file"));
	cfg.frequency_file	= std::filesystem::path(requiredString("frequency_file"));

	// Optional booleans with defaults
	if (config["output"]) cfg.output = config["output"].as<bool>();

	if (config["output_overwrite_prevention"])
		cfg.output_overwrite_prevention = config["output_overwrite_prevention"].as<bool>();

	// Optional string (converted to filesystem::path)
	if (config["output_directory"])
		cfg.output_directory = std::filesystem::path(config["output_directory"].as<std::string>());

	// Emissivity model (required)
	cfg.emissivity_model = config["emissivity_model"].as<emissivity_model_t>();

	// Flags
	if (config["use_B"]) cfg.use_B = config["use_B"].as<bool>();
	if (config["use_1_5D_approx"]) cfg.use_1_5D_approx = config["use_1_5D_approx"].as<bool>();
	if (config["enable_continuum"]) cfg.enable_continuum = config["enable_continuum"].as<bool>();

	// Formal solver
	if (config["formal_solver"])
	{
		std::string fs	  = config["formal_solver"].as<std::string>();
		cfg.formal_solver = validateFormalSolver(fs);
	}

	// Integers
	if (config["N_theta"]) cfg.N_theta = config["N_theta"].as<int>();
	if (config["N_chi"]) cfg.N_chi = config["N_chi"].as<int>();
	if (config["N_x"]) cfg.N_x = config["N_x"].as<int>();
	if (config["N_y"]) cfg.N_y = config["N_y"].as<int>();
	if (config["L"]) cfg.L = config["L"].as<double>();

	// Optional strings (converted to filesystem::path)
	if (config["input_cul"]) cfg.input_cul = std::filesystem::path(config["input_cul"].as<std::string>());
	if (config["input_qel"]) cfg.input_qel = std::filesystem::path(config["input_qel"].as<std::string>());
	if (config["input_llp"]) cfg.input_llp = std::filesystem::path(config["input_llp"].as<std::string>());
	if (config["input_back"]) cfg.input_back = std::filesystem::path(config["input_back"].as<std::string>());

	// use_prec (default is true)
	if (config["use_prec"]) cfg.use_prec = config["use_prec"].as<bool>();

	// override logic: do not use preconditioner for CRD and ZERO
	switch (cfg.emissivity_model)
	{
		case emissivity_model_t::CRD_limit:
		case emissivity_model_t::CRD_limit_VHP:
		case emissivity_model_t::ZERO:
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
		if (s["ksp_rtol"]) cfg.solver.ksp_rtol = s["ksp_rtol"].as<double>();
		if (s["ksp_max_it"]) cfg.solver.ksp_max_it = s["ksp_max_it"].as<int>();
		if (s["gmres_restart"]) cfg.solver.gmres_restart = s["gmres_restart"].as<int>();
	}

	// Preconditioner section
	if (config["prec"])
	{
		auto p = config["prec"];

		if (p["pc_solver_type"]) cfg.prec.pc_solver_type = toKSPType(p["pc_solver_type"].as<std::string>());
		if (p["pc_rtol"]) cfg.prec.pc_rtol = p["pc_rtol"].as<double>();
		if (p["pc_max_it"]) cfg.prec.pc_max_it = p["pc_max_it"].as<int>();
	}

	// Arbitrary beams section
	if (config["arbitrary_beams"])
	{
		for (const auto &beam : config["arbitrary_beams"])
		{
			BeamDirection bd;
			if (beam["mu"]) bd.mu = beam["mu"].as<double>();
			if (beam["chi"]) bd.chi = beam["chi"].as<double>();
			cfg.arbitrary_beams.push_back(bd);
		}
	}

	return cfg;
}

void
writeConfigResume(const AppConfig &cfg, std::ostream &os)
{
	const int label_width = 36;
	auto	  print		  = [&](const std::string &label, const std::string &value)
	{ os << std::left << std::setw(label_width) << (label + ":") << value << std::endl; };

	os << "YAML Configuration Summary:" << std::endl;
	print("Input Directory", cfg.input_directory.string());
	print("Input File", cfg.input_file.string());
	print("Frequency File", cfg.frequency_file.string());
	print("Output Directory", cfg.output_directory.string());
	print("Output Enabled", (cfg.output ? "Yes" : "No"));
	print("Output Overwrite Prevention", (cfg.output_overwrite_prevention ? "Yes" : "No"));
	print("Emissivity Model", emissivity_model_to_string_long(cfg.emissivity_model));
	print("Use Magnetic Field", (cfg.use_B ? "Yes" : "No"));
	print("Enable Continuum", (cfg.enable_continuum ? "Yes" : "No"));
	print("Use 1.5D Approximation", (cfg.use_1_5D_approx ? "Yes" : "No"));
	print("Formal Solver", cfg.formal_solver);
	print("Angular Grid - N_theta", std::to_string(cfg.N_theta));
	print("Angular Grid - N_chi", std::to_string(cfg.N_chi));
	print("Horizontal Grid - N_x", std::to_string(cfg.N_x));
	print("Horizontal Grid - N_y", std::to_string(cfg.N_y));
	print("Horizontal Grid - L", std::to_string(cfg.L));
	print("Input CUL File", cfg.input_cul.string());
	print("Input QEL File", cfg.input_qel.string());
	print("Input LLP File", cfg.input_llp.string());
	print("Input BACK File", cfg.input_back.string());
	print("Use Preconditioner", (cfg.use_prec ? "Yes" : "No"));

	os << std::endl << "Solver Configuration:" << std::endl;
	print("KSP Solver Type", KSPTypeToString(cfg.solver.ksp_solver_type));
	print("KSP Relative Tolerance", std::to_string(cfg.solver.ksp_rtol));
	print("KSP Maximum Iterations", std::to_string(cfg.solver.ksp_max_it));
	print("GMRES Restart", std::to_string(cfg.solver.gmres_restart));

	os << std::endl << "Preconditioner Configuration:" << std::endl;
	print("PC Solver Type", KSPTypeToString(cfg.prec.pc_solver_type));
	print("PC Relative Tolerance", std::to_string(cfg.prec.pc_rtol));
	print("PC Maximum Iterations", std::to_string(cfg.prec.pc_max_it));

	os << std::endl;
}

// std::string get_arg(const std::string &input, const std::string &word) {
//   std::regex pattern("^\\s+" + word + ":\\s+(.*)$", std::regex::multiline);
//   std::smatch match;

//   if (std::regex_search(input, match, pattern)) {
//     return match.str(1);
//   }

//   return std::string();
// }

std::string
get_arg(const std::string &input, const std::string &word)
{
	std::regex		   pattern("^\\s*" + word + ":\\s*([^\\s]*)\\s*$");
	std::smatch		   match;
	std::istringstream stream(input);
	std::string		   line;

	while (std::getline(stream, line))
	{
		if (std::regex_search(line, match, pattern))
		{
			return match.str(1);
		}
	}

	return std::string();
}

std::map<std::string, std::string>
get_input_PORTA(const std::filesystem::path &config_file, int mpi_rank)
{
	std::map<std::string, std::string> input_PORTA_map;

	std::ifstream input_file(config_file.string());

	if (input_file)
	{
		// read the input file
		std::stringstream buffer;
		buffer << input_file.rdbuf();
		std::string input_string = buffer.str();

		std::string pmd_file  = get_arg(input_string, "pmd");
		std::string cul_file  = get_arg(input_string, "cul");
		std::string qel_file  = get_arg(input_string, "qel");
		std::string llp_file  = get_arg(input_string, "llp");
		std::string back_file = get_arg(input_string, "back");

		if (pmd_file.empty())
		{
			if (mpi_rank == 0) std::cerr << "Error reading PORTA input file: " << config_file << std::endl;
			exit(1);
		}

		input_PORTA_map["pmd"]	= pmd_file;
		input_PORTA_map["cul"]	= cul_file;
		input_PORTA_map["qel"]	= qel_file;
		input_PORTA_map["llp"]	= llp_file;
		input_PORTA_map["back"] = back_file;
	}
	else
	{
		if (mpi_rank == 0) std::cerr << "Error in opening PORTA input file: " << config_file << std::endl;
		exit(1);
	}

	input_file.close();
	return input_PORTA_map;
}

std::string
getCurrentDateTime()
{
	auto now	   = std::chrono::system_clock::now();
	auto in_time_t = std::chrono::system_clock::to_time_t(now);

	std::stringstream ss;
	ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %H:%M:%S");
	return ss.str();
}

// compile and run with:
// make -j 32 && srun -n 4 ./main

std::string
getOptionArgument(int argc, char *argv[], const std::string &option, const std::string &defaultValue)
{
	for (int i = 1; i < argc; ++i)
	{
		std::string arg = argv[i];
		if (arg == option && i + 1 < argc)
		{
			return argv[i + 1];
		}
	}
	return defaultValue; // Return default when not found
}

bool
getOptionFlag(int argc, char *argv[], const std::string &option)
{
	for (int i = 1; i < argc; ++i)
	{
		std::string arg = argv[i];
		if (arg == option)
		{
			return true;
		}
	}
	return false;
}

// TODO
void
print_help()
{
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

int
acc_devices_print_info(const int mpi_rank, const int mpi_size, std::ostream &os)
{
#if ACC_SOLAR_3D == _ON_

	for (int i = 0; i < LIMIT_OUT_DEVICE_MEMORY_USAGE; i++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		if (mpi_rank == i && RII_epsilon_contrib::RII_contrib_MPI_Is_Device_Handler())
		{
			os << RII_epsilon_contrib::RII_contrib_MPI_Device_Info() << std::endl;
		}
	}

	const int is_device_handler = int(RII_epsilon_contrib::RII_contrib_MPI_Is_Device_Handler());

	int devices_cnt = 0;
	MPI_Reduce(&is_device_handler, &devices_cnt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (mpi_rank == 0)
	{
		os << "Using ACC version with " << devices_cnt << " devices." << std::endl;
	}

	return devices_cnt;
#else
	return 0;
#endif
}
