#include "RT_arbitrary_beam.hpp"
#include "RT_solver.hpp"

//////////////////////////////////////////////////////////////////////////
// Main function
//////////////////////////////////////////////////////////////////////////
int
main(int argc, char *argv[])
{
	MPI_CHECK(MPI_Init(&argc, &argv));
	print_PETSc_mem();

	const double main_start_time = MPI_Wtime();

	int mpi_rank, mpi_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	// print help if requested by the user
	if (getOptionFlag(argc, argv, "--help") && mpi_rank == 0)
	{
		print_help();
		return 0;
	}

	// import input file config.yml
	const std::string yaml_config_file = getOptionArgument(argc, argv, "--yaml_config", "config.yml");
	AppConfig		  cfg;
	try
	{
		if (mpi_rank == 0) std::cout << "Loading config file: " << yaml_config_file << std::endl;
		cfg = loadConfig(yaml_config_file);
	}
	catch (const std::exception &e)
	{
		std::cerr << "Config error: " << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	std::stringstream	  ss_a, ss_b;
	std::filesystem::path output_info_file_name;

	// print some info
	if (mpi_rank == 0)
	{
		ss_a << "Date and time: " << getCurrentDateTime() << std::endl;
		ss_a << std::endl << "Command line arguments: " << std::endl;
		for (int i = 0; i < argc; ++i) ss_a << argv[i] << " ";
		ss_a << std::endl;

		ss_a << "PetscInt size: " << sizeof(PetscInt) << " bytes; " << (sizeof(PetscInt) * 8) << " bits." << std::endl
			 << std::endl;
		std::cout << ss_a.str();
	}

	PetscInitialize(&argc, &argv, (char *)0, NULL);
	// Kokkos::initialize(argc, argv);

#if ACC_SOLAR_3D == _ON_

#pragma message("ACC version ENABLED")

	const int devices_per_node = 4; // for MareNostrum 5 ACC & Daint Alps.
	int		  mpi_error		   = RII_epsilon_contrib::RII_contrib_MPI_Init(devices_per_node, MPI_COMM_WORLD);

	set_RII_contrib_block_size(RII_CONTRIB_BLOCK_SIZE);

	if (RII_epsilon_contrib::RII_contrib_MPI_Is_Device_Handler()) //
		RII_epsilon_contrib::RII_contrib_MPI_Init_Memory_Pool();  //

	const int devices_cnt =									   //
		acc_devices_print_info(mpi_rank, mpi_size, std::cout); //

#else

	set_RII_contrib_block_size(1);

#endif // ACC_SOLAR_3D

	if (mpi_rank == 0)
	{
		writeConfigResume(cfg, ss_b);
	}

	{ // start scope for RT_problem and RT_solver

		// create problem object
		auto rt_problem_ptr = std::make_shared<RT_problem>(cfg);
		rt_problem_ptr->space_grid_->dumpDecompositionCSV("decomposition.csv");

		// create solver object
		RT_solver rt_solver(rt_problem_ptr);

		//////////////////////////////////////////////////////////////////////////
		// Prepare output directory
		// If the output directory does not exist, create it. If it exists, abort !
		const std::filesystem::path output_subdir =
			std::filesystem::path(cfg.input_file.string() + "." + emissivity_model_to_string_long(cfg.emissivity_model));
		const std::filesystem::path output_path =
			std::filesystem::path(cfg.output_directory) / std::filesystem::path(output_subdir);
		std::string output_file = (output_path / "profiles").string();

		if (mpi_rank == 0)
		{
			std::cout << "Output path: " << output_path << std::endl;
		}

		if (cfg.output) // Set up output
		{
			if (mpi_rank == 0)
			{
				if (not std::filesystem::exists(output_path))
				{
					std::filesystem::create_directories(output_path);
				}
				else if (cfg.output_overwrite_prevention)
				{
					std::cerr << "File: " << __FILE__ << " Line: " << __LINE__ << " Directory: " << output_path
							  << " already exists" << std::endl;
					std::cerr << "Use different output directory" << std::endl;
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
			}

			if (mpi_rank == 0)
			{
				ss_b << "Output directory: " << output_path << std::endl;
				std::cout << ss_b.str();
				std::cout.flush();

				output_info_file_name = output_path / "info.txt"; // <-- Remove 'const auto', assign to outer variable
				std::ofstream output_file_info(output_info_file_name);

				output_file_info << ss_a.str();
				output_file_info << ss_b.str();
				output_file_info.close();
			}
		} // end if (output)

		// print some timing info
		const double main_setup_time = MPI_Wtime();
		if (mpi_rank == 0)
		{
			std::stringstream ss_mem;
			ss_mem << std::fixed << std::setprecision(2);
			ss_mem << "Setup time: " << (main_setup_time - main_start_time) << " seconds." << std::endl;

			std::cout << ss_mem.str();
			if (cfg.output)
			{
				std::ofstream output_file_info(output_info_file_name, std::ios_base::app);
				output_file_info << ss_mem.str();
				output_file_info.close();
			}
		}

		///////////////////////////////////////////////////
		// solve //////////////////////////////////////////
		rt_solver.solve(); 
		// rt_solver.apply_formal_solver();

		const double main_solve_end_time = MPI_Wtime();

		// print some memory info
		print_PETSc_mem();

		double time_output_csv;
		double time_output_hdf5_emergent;
		double time_output_hdf5_JKQ;
		if (cfg.output)
		{
			// write output Surface profiles for all surface points
			const int N_x = rt_problem_ptr->N_x_;
			const int N_y = rt_problem_ptr->N_y_;

			// Old school output ....
			// Where the output is written in separate files for each emergent spatial point
			////////////////////////////////////////////////////////////////////////////////////////////////
			time_output_csv = MPI_Wtime();
			for (int i = 0; i < N_x; ++i)
			{
				for (int j = 0; j < N_y; ++j)
				{
					// const std::string output_file_Omega_mu	   = output_file + "_mu" + mu_str + "_chi" + chi_str;

					rt_problem_ptr->write_surface_point_profiles(output_file, i, j);

					// TESTING purpose: write the CSV files too
					rt_problem_ptr->write_surface_point_profiles_csv(output_file, i, j);

					if (rt_problem_ptr->mpi_rank_ == 0 and i == 0 and j == 0)
					{
						const std::string output_file_frequencies  = (output_path / "frequencies_grid_Hz").string();
						const std::string output_file_angular_grid = (output_path / "angular_grid").string();

						rt_problem_ptr->write_angular_grid_csv(output_file_angular_grid, i, j);
						rt_problem_ptr->write_frequencies_grid_csv(output_file_frequencies, i, j);
					}
				}
			}
			time_output_csv = MPI_Wtime() - time_output_csv;

			// New HDF5 output .... TESTING PASSED
			////////////////////////////////////////////////////////////////////////////////////////////////
			time_output_hdf5_emergent = MPI_Wtime();
			write_emergent_field_hdf5(*rt_problem_ptr,											  //
									  (output_path / "emergent_field_angular_grid.h5").string()); //
			time_output_hdf5_emergent = MPI_Wtime() - time_output_hdf5_emergent;
			time_output_hdf5_JKQ = MPI_Wtime();
			rt_problem_ptr->write_JKQ_field_hdf5((output_path / "JKQ_field.h5").string()); //
			time_output_hdf5_JKQ = MPI_Wtime() - time_output_hdf5_JKQ;

			////////////////////////////////////////////////////////////////////////////////////////////////
			// compute arbitrary beams
			////////////////////////////////////////////////////////////////////////////////////////////////
			/// Set arbitrary beam directions
			// free some memory
			rt_problem_ptr->free_fields_memory();
			rt_solver.free_fields_memory();

			process_arbitrary_beams_hdf(cfg.arbitrary_beams,										 //
										rt_solver,													 //
										rt_problem_ptr,												 //
										(output_path / "emergent_field_abitrary_Omega.h5").string(), //
										output_info_file_name);										 //

			const double main_end_time = MPI_Wtime();

			////////////////////////////////////////////////////////////////////////////////////////////////
#if ACC_SOLAR_3D == _ON_
			if (RII_epsilon_contrib::RII_contrib_MPI_Is_Device_Handler() and //
				mpi_rank < LIMIT_OUT_DEVICE_MEMORY_USAGE)					 //
			{
				std::string pool_info = RII_epsilon_contrib::RII_contrib_MPI_Generate_DPool_Report(true, true);

				const auto file_name_pool = output_path /													  //
											(boost::format("device_pool_info_rank_%d.txt") % mpi_rank).str(); //

				std::ofstream myfile;
				myfile.open(file_name_pool.string());
				if (myfile.is_open())
				{
					myfile << pool_info;
					myfile.close();
				}
				else
				{
					std::cerr << "Unable to open file: " << file_name_pool << std::endl;
				}
			}
#endif // ACC_SOLAR_3D

			// print memory usage
			const double byte_to_GB = 1.0 / (1000 * 1024 * 1024);

			unsigned long long vm_usage;
			unsigned long long resident_set;

			rii::process_mem_usage(vm_usage, resident_set);

			// print_PETSc_mem();

			if (mpi_rank == 0)
			{
				std::stringstream ss_mem;
				ss_mem << "Total memory usage (vm_usage):               " << byte_to_GB * vm_usage << " GB" << std::endl;
				ss_mem << "Total memory usage (resident_set):           " << byte_to_GB * resident_set << " GB"
					   << std::endl;

				// ss_mem << mem_petsc << std::endl << std::endl;

#if ACC_SOLAR_3D == _ON_
				ss_mem << "Total number of devices (accelerators) used: " << devices_cnt << std::endl;
#endif // ACC_SOLAR_3D

				ss_mem << format_execution_time_summary(rt_problem_ptr->mpi_size_, //
														main_start_time,		   //
														main_setup_time,		   //
														main_solve_end_time,	   //
														main_end_time,
														time_output_csv,
														time_output_hdf5_emergent,
														time_output_hdf5_JKQ);			   //

				std::cout << ss_mem.str();

				if (!output_info_file_name.empty()) // <-- Add this check
				{
					std::ofstream output_file_info(output_info_file_name, std::ios::app);
					output_file_info << ss_mem.str();
					output_file_info.close();
				} // end if !output_info_file_name.empty()
			} // end if mpi_rank_ == 0
		} // end if (cfg.output)
	} // end scope for RT_problem and RT_solver

#if ACC_SOLAR_3D == _ON_
	if (RII_epsilon_contrib::RII_contrib_MPI_Is_Device_Handler())
		RII_epsilon_contrib::RII_contrib_MPI_Finalize_Memory_Pool();

	RII_epsilon_contrib::RII_contrib_MPI_Finalize();
#endif // ACC_SOLAR_3D

	// Kokkos::finalize();
	PetscFinalize(); // CHKERRQ(ierr);

	MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank == 0) std::cout << "TRIP ended successfully at: " << getCurrentDateTime() << std::endl;

	MPI_CHECK(MPI_Finalize());

	return 0;
} // end main
//////////////////////////////////////////////////////////////////////////
