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
	if (name == "KSPRICHARDSON") return KSPRICHARDSON;

	throw std::runtime_error("Unknown KSPType: " + name);
}

inline std::string
KSPTypeToString(KSPType type)
{
	PetscBool match;

	PetscStrcmp(type, KSPGMRES, &match);
	if (match) return "KSPGMRES";

	PetscStrcmp(type, KSPFGMRES, &match);
	if (match) return "KSPFGMRES";

	PetscStrcmp(type, KSPPGMRES, &match);
	if (match) return "KSPPGMRES";

	PetscStrcmp(type, KSPBCGS, &match);
	if (match) return "KSPBCGS";

	PetscStrcmp(type, KSPPREONLY, &match);
	if (match) return "KSPPREONLY";

	PetscStrcmp(type, KSPRICHARDSON, &match);
	if (match) return "KSPRICHARDSON";

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

	if (config["write_whole_3D_field_hdf5"])
		cfg.write_whole_3D_field_hdf5 = config["write_whole_3D_field_hdf5"].as<bool>();

	if (config["write_text_output"]) cfg.write_text_output = config["write_text_output"].as<bool>();

	// Optional string (converted to filesystem::path)
	if (config["output_directory"])
		cfg.output_directory = std::filesystem::path(config["output_directory"].as<std::string>());

	// Optional string for testing
	if (config["reference_sol_directory"])
		cfg.reference_sol_directory = std::filesystem::path(config["reference_sol_directory"].as<std::string>());

	// printing flag
	if (config["verbose"]) cfg.verbose = config["verbose"].as<bool>();

	// Emissivity model (required)
	cfg.emissivity_model = config["emissivity_model"].as<emissivity_model_t>();

	// Flags
	if (config["use_B"]) cfg.use_B = config["use_B"].as<bool>();
	if (config["use_1_5D_approx"]) cfg.use_1_5D_approx = config["use_1_5D_approx"].as<bool>();
	if (config["enable_continuum"]) cfg.enable_continuum = config["enable_continuum"].as<bool>();
	if (config["use_Vb"]) cfg.use_Vb = config["use_Vb"].as<bool>();

	if (config["use_D2_from_input"]) cfg.use_D2_from_input = config["use_D2_from_input"].as<bool>();

	// Formal solver
	if (config["formal_solver"])
	{
		std::string fs	  = config["formal_solver"].as<std::string>();
		cfg.formal_solver = validateFormalSolver(fs);
	}

	// Set constant B
	if (config["set_uniform_B"]) cfg.set_uniform_B = config["set_uniform_B"].as<bool>();
	if (config["B_field"])
	{
		auto B_field_node = config["B_field"];
		if (B_field_node.IsSequence() && B_field_node.size() == 3)
		{
			for (size_t i = 0; i < 3; ++i)
			{
				cfg.B_field[i] = B_field_node[i].as<double>();
			}
		}
		else
		{
			throw std::runtime_error("B_field must be a sequence of three numbers.");
		}
	}

	if (config["set_uniform_Vb"]) cfg.set_uniform_Vb = config["set_uniform_Vb"].as<bool>();
	if (config["Vb_field"])
	{
		auto Vb_field_node = config["Vb_field"];
		if (Vb_field_node.IsSequence() && Vb_field_node.size() == 3)
		{
			for (size_t i = 0; i < 3; ++i)
			{
				cfg.Vb_field[i] = Vb_field_node[i].as<double>();
			}
		}
		else
		{
			throw std::runtime_error("Vb_field must be a sequence of three numbers.");
		}
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
		if (s["ksp_use_J_KQ"]) cfg.solver.ksp_use_J_KQ = s["ksp_use_J_KQ"].as<bool>();
	}

	// Preconditioner section
	if (config["prec"])
	{
		auto p = config["prec"];

		if (p["pc_solver_type"]) cfg.prec.pc_solver_type = toKSPType(p["pc_solver_type"].as<std::string>());
		if (p["pc_rtol"]) cfg.prec.pc_rtol = p["pc_rtol"].as<double>();
		if (p["pc_max_it"]) cfg.prec.pc_max_it = p["pc_max_it"].as<int>();
		if (p["pc_use_J_KQ"]) cfg.prec.pc_use_J_KQ = p["pc_use_J_KQ"].as<bool>();
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
	const int label_width = 44;
	auto	  print		  = [&](const std::string &label, const std::string &value)
	{ os << std::left << std::setw(label_width) << (label + ":") << value << std::endl; };

	auto to_sci = [](double v, int prec = 2)
	{
		std::ostringstream oss;
		oss << std::scientific << std::setprecision(prec) << v;
		return oss.str();
	};

	os << "YAML Configuration Summary:" << std::endl;
	print("Input Directory", cfg.input_directory.string());
	print("Input File", cfg.input_file.string());
	print("Frequency File", cfg.frequency_file.string());
	print("Output Directory", cfg.output_directory.string());
	print("Output Enabled", (cfg.output ? "Yes" : "No"));
	print("Output Overwrite Prevention", (cfg.output_overwrite_prevention ? "Yes" : "No"));
	print("Write Whole 3D Field (HDF5)", (cfg.write_whole_3D_field_hdf5 ? "Yes" : "No"));
	print("Write Text Output", (cfg.write_text_output ? "Yes" : "No"));
	print("Reference Solution Directory", cfg.reference_sol_directory.string());
	print("Verbose", (cfg.verbose ? "Yes" : "No"));
	print("Emissivity Model", emissivity_model_to_string_long(cfg.emissivity_model));
	print("Use Magnetic Field", (cfg.use_B ? "Yes" : "No"));
	print("Use Bulk Velocity", (cfg.use_Vb ? "Yes" : "No"));
	print("Use D2 from Input", (cfg.use_D2_from_input ? "Yes" : "No"));
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
	print("Set Uniform Magnetic Field", (cfg.set_uniform_B ? "Yes" : "No"));
	if (cfg.set_uniform_B)
	{
		print("* Uniform Magnetic Field Value (Gauss)", std::to_string(cfg.B_field[0]));
		print("* Uniform Magnetic Field Theta (rad)", std::to_string(cfg.B_field[1]));
		print("* Uniform Magnetic Field Chi (rad)", std::to_string(cfg.B_field[2]));
	}

	print("Set Uniform Bulk Velocity", (cfg.set_uniform_Vb ? "Yes" : "No"));
	if (cfg.set_uniform_Vb)
	{
		print("* Uniform Bulk Velocity Vx (cm/s)", std::to_string(cfg.Vb_field[0]));
		print("* Uniform Bulk Velocity Vy (cm/s)", std::to_string(cfg.Vb_field[1]));
		print("* Uniform Bulk Velocity Vz (cm/s)", std::to_string(cfg.Vb_field[2]));
	}

	print("Arbitrary Beams Count", std::to_string(cfg.arbitrary_beams.size()));
	for (size_t i = 0; i < cfg.arbitrary_beams.size(); ++i)
	{
		const auto		 &beam		   = cfg.arbitrary_beams[i];
		const bool		  is_last_beam = (i + 1 == cfg.arbitrary_beams.size());
		const std::string beam_prefix  = is_last_beam ? "" : "";
		print(beam_prefix + "* Beam", std::to_string(i));
		print("   > mu", std::to_string(beam.mu));
		print("   > chi (rad)", std::to_string(beam.chi));
	}

	os << std::endl << "Solver Configuration:" << std::endl;
	print("KSP Solver Type", KSPTypeToString(cfg.solver.ksp_solver_type));
	print("KSP Relative Tolerance", to_sci(cfg.solver.ksp_rtol));
	print("KSP Maximum Iterations", std::to_string(cfg.solver.ksp_max_it));
	print("KSP Use J/KQ", (cfg.solver.ksp_use_J_KQ ? "Yes" : "No"));
	print("GMRES Restart", std::to_string(cfg.solver.gmres_restart));

	os << std::endl << "Preconditioner Configuration:" << std::endl;
	print("PC Solver Type", KSPTypeToString(cfg.prec.pc_solver_type));
	print("PC Relative Tolerance", to_sci(cfg.prec.pc_rtol));
	print("PC Maximum Iterations", std::to_string(cfg.prec.pc_max_it));
	print("PC Use J/KQ", (cfg.prec.pc_use_J_KQ ? "Yes" : "No"));

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

	for (int i = 0; i < 2 /* LIMIT_OUT_DEVICE_MEMORY_USAGE */; i++)
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

namespace TRIP_Comms
{
	std::shared_ptr<MPI_TRIP_Communicators> TRIP_comms_shared_ptr = nullptr;

	void
	initialize_MPI_TRIP_Communicators()
	{
		if (!TRIP_comms_shared_ptr)
		{
			TRIP_comms_shared_ptr = MPI_TRIP_Communicators::make_shared();
		}
	}

	void
	finalize_MPI_TRIP_Communicators()
	{
		TRIP_comms_shared_ptr.reset();
	}

	MPI_TRIP_Communicators::MPI_TRIP_Communicators_shared_ptr_t
	getTRIPCommunicators()
	{
		return TRIP_comms_shared_ptr;
	}
} // namespace TRIP_Comms

double
calculate_D2(const double a, const double b, const double nH, const double T, const double T_ref)
{
	// From: Collisional Depolarization of the Solar Ca, Mg, and Ba Levels. M. Derouich 2019. DOI 10.3847/1538-4357/ab26b4
	// eq. 6 in the paper.
	// In the paper, T_ref is set to 5000 K
	// a = 1.24 * 10^-8
	// b = 0.37
	// In our method a and b are given in the h5 input file.

	double D2 = a * nH * std::pow(T / T_ref, b);

	return D2;
}