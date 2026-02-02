#include "RT_problem.hpp"

//////////////////////////////////////////////////////////////////////////
// Write surface point profiles to file
// write_surface_point_profiles
//////////////////////////////////////////////////////////////////////////
void 
RT_problem::write_surface_point_profiles(input_string file_name, const int i_space, const int j_space) const
{
	// // a single MPI rank writes output
	// if (mpi_rank_ == 0) std::cout << " Writing output in spatial point (" << i_space << ", " << j_space << ")" <<
	// std::endl;
	// indeces
	const int i_start = space_grid_->getGlobalStartX(); //g_dev.margin[0];
	const int j_start = space_grid_->getGlobalStartY(); //g_dev.margin[1];
	const int k_start = space_grid_->getGlobalStartZ(); //g_dev.margin[2];

	const int i_end = i_start + space_grid_->getLocalSizeX(); //g_dev.dim[0];
	const int j_end = j_start + space_grid_->getLocalSizeY(); //g_dev.dim[1];

	int i_global, j_global;

	double I, QUV;

	// write profiles
	if (space_grid_->local_to_global_coordinate(2, k_start) == 0)
	{
		for (int i = i_start; i < i_end; ++i)
		{
			i_global = space_grid_->local_to_global_coordinate(0, i);

			if (i_global == i_space)
			{
				for (int j = j_start; j < j_end; ++j)
				{
					j_global = space_grid_->local_to_global_coordinate(1, j);

					if (j_global == j_space)
					{
						// Create a new file
						input_string output_file =
							file_name + "_" + std::to_string(i_space) + "_" + std::to_string(j_space) + ".m";
						std::ofstream outputFile(output_file);

						if (outputFile.is_open())
						{
							// write grids
							outputFile << "\nnu_grid_ = [ ";
							for (int n = 0; n < N_nu_; ++n)
								outputFile << std::scientific << std::setprecision(15) << nu_grid_[n] << " ";
							outputFile << "];\n";

							outputFile << "\ntheta_grid = [ ";
							for (int j_theta = 0; j_theta < N_theta_; ++j_theta)
								outputFile << std::scientific << std::setprecision(15) << theta_grid_[j_theta] << " ";
							outputFile << "];\n";

							outputFile << "\nmu_grid = [ ";
							for (int j_theta = 0; j_theta < N_theta_; ++j_theta)
								outputFile << std::scientific << std::setprecision(15) << mu_grid_[j_theta] << " ";
							outputFile << "];\n";

							outputFile << "\nchi_grid = [ ";
							for (int k_chi = 0; k_chi < N_chi_; ++k_chi)
								outputFile << std::scientific << std::setprecision(15) << chi_grid_[k_chi] << " ";
							outputFile << "];\n";

							// create MATLAB data structure
							outputFile << std::scientific << std::setprecision(15) << "\nField = cell(4," << N_theta_
									   << "," << N_chi_ << ");" << std::endl;

							for (int j_theta = N_theta_ / 2; j_theta < N_theta_; ++j_theta)
							{
								for (int k_chi = 0; k_chi < N_chi_; ++k_chi)
								{
									const int b_start = I_field_->local_to_block(j_theta, k_chi, 0);

									for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
									{
										outputFile << "\nField{" << i_stokes + 1 << "," << j_theta + 1 << "," << k_chi + 1
												   << "} = [ ";

										for (int b = 0; b < 4 * N_nu_; b = b + 4)
										{
											I = I_field_->block(i, j, k_start)[b_start + b];

											if (i_stokes == 0)
											{
												outputFile << std::scientific << std::setprecision(15) << I << " ";
											}
											else
											{
												QUV = I_field_->block(i, j, k_start)[b_start + b + i_stokes];

												outputFile << std::scientific << std::setprecision(15) << 100.0 * QUV / I
														   << " ";
											}
										}

										outputFile << "];\n";
									}
								}
							}

							// Close the file
							outputFile.close();
						}
						else
						{
							if (mpi_rank_ == 0) std::cout << "\nERROR: failed to create the output file." << std::endl;
						}
					}
				}
			}
		}
	}

	if (mpi_rank_ == 0 and i_space == 0 and j_space == 0)
		std::cout << "Output written in " << file_name << "_i_j.m" << "\n" << std::endl;
}

//////////////////////////////////////////////////////////////////////////
// Write surface point angular grid to file
// write_angular_grid_csv
//////////////////////////////////////////////////////////////////////////
void 
RT_problem::write_angular_grid_csv(input_string file_name, const int i_space, const int j_space,
								   const unsigned int precision) const
{
	// indeces
	const int i_start = space_grid_->getGlobalStartX(); //g_dev.margin[0];
	const int j_start = space_grid_->getGlobalStartY(); //g_dev.margin[1];
	const int k_start = space_grid_->getGlobalStartZ(); //g_dev.margin[2];

	const int i_end = i_start + space_grid_->getLocalSizeX(); //g_dev.dim[0];
	const int j_end = j_start + space_grid_->getLocalSizeY(); //g_dev.dim[1];


	int i_global, j_global;

	double I, QUV;

	// write profiles
	if (space_grid_->local_to_global_coordinate(2, k_start) == 0)
	{
		for (int i = i_start; i < i_end; ++i)
		{
			i_global = space_grid_->local_to_global_coordinate(0, i);

			if (i_global == i_space)
			{
				for (int j = j_start; j < j_end; ++j)
				{
					j_global = space_grid_->local_to_global_coordinate(1, j);

					if (j_global == j_space)
					{
						// Create a new file
						input_string output_file =
							file_name + "_" + std::to_string(i_space) + "_" + std::to_string(j_space) + ".csv";
						std::ofstream outputFile(output_file);

						if (outputFile.is_open())
						{
							outputFile << "theta_i,chi_i,theta,mu,chi" << std::endl;
							// write grids
							for (int j_theta = N_theta_ / 2; j_theta < N_theta_; ++j_theta)
							{
								for (int k_chi = 0; k_chi < N_chi_; ++k_chi)
								{
									outputFile << std::scientific << std::setprecision(precision) << j_theta << ","
											   << k_chi << "," << theta_grid_[j_theta] << "," << mu_grid_[j_theta] << ","
											   << chi_grid_[k_chi] << std::endl;
								}
							}
						}
						else
						{
							if (mpi_rank_ == 0) std::cout << "\nERROR: failed to create the output file." << std::endl;
						}
					}
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////
// Write surface point frequency grid to file
// write_frequencies_grid_csv
//////////////////////////////////////////////////////////////////////////
void 
RT_problem::write_frequencies_grid_csv(input_string file_name, const int i_space, const int j_space,
									   const unsigned int precision) const
{
	// indeces
	const int i_start = space_grid_->getGlobalStartX(); //g_dev.margin[0];
	const int j_start = space_grid_->getGlobalStartY(); //g_dev.margin[1];
	const int k_start = space_grid_->getGlobalStartZ(); //g_dev.margin[2];

	const int i_end = i_start + space_grid_->getLocalSizeX(); //g_dev.dim[0];
	const int j_end = j_start + space_grid_->getLocalSizeY(); //g_dev.dim[1];

	int i_global, j_global;

	double I, QUV;

	// write profiles
	if (space_grid_->local_to_global_coordinate(2, k_start) == 0)
	{
		for (int i = i_start; i < i_end; ++i)
		{
			i_global = space_grid_->local_to_global_coordinate(0, i);

			if (i_global == i_space)
			{
				for (int j = j_start; j < j_end; ++j)
				{
					j_global = space_grid_->local_to_global_coordinate(1, j);

					if (j_global == j_space)
					{
						// Create a new file
						input_string output_file =
							file_name + "_" + std::to_string(i_space) + "_" + std::to_string(j_space) + ".csv";

						std::ofstream outputFile(output_file);

						if (outputFile.is_open())
						{
							// write grids
							outputFile << "";
							for (int n = 0; n < N_nu_; ++n)
							{
								std::string sep = n < (N_nu_ - 1) ? "," : "";
								outputFile << std::scientific << std::setprecision(15) << nu_grid_[n] << sep;
							}
							outputFile << std::endl;
						}
					}
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////
// Write surface point profiles to file in CSV format
// write_surface_point_profiles_csv
//////////////////////////////////////////////////////////////////////////
void 
RT_problem::write_surface_point_profiles_csv(input_string file_name, const int i_space, const int j_space,
											 const unsigned int precision) const
{
	// indeces
	const int i_start = space_grid_->getGlobalStartX(); //g_dev.margin[0];
	const int j_start = space_grid_->getGlobalStartY(); //g_dev.margin[1];
	const int k_start = space_grid_->getGlobalStartZ(); //g_dev.margin[2];

	const int i_end = i_start + space_grid_->getLocalSizeX(); //g_dev.dim[0];
	const int j_end = j_start + space_grid_->getLocalSizeY(); //g_dev.dim[1];

	int i_global, j_global;

	double I, QUV;

	// write profiles
	if (space_grid_->local_to_global_coordinate(2, k_start) == 0)
	{
		for (int i = i_start; i < i_end; ++i)
		{
			i_global = space_grid_->local_to_global_coordinate(0, i);

			if (i_global == i_space)
			{
				for (int j = j_start; j < j_end; ++j)
				{
					j_global = space_grid_->local_to_global_coordinate(1, j);

					if (j_global == j_space)
					{
						// Create a new file
						input_string output_file =
							file_name + "_" + std::to_string(i_space) + "_" + std::to_string(j_space) + ".csv";
						std::ofstream outputFile(output_file);

						if (outputFile.is_open())
						{
							for (int j_theta = N_theta_ / 2; j_theta < N_theta_; ++j_theta)
							{
								for (int k_chi = 0; k_chi < N_chi_; ++k_chi)
								{
									const int b_start = I_field_->local_to_block(j_theta, k_chi, 0);

									for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
									{
										// outputFile << "\nField{" << i_stokes + 1 << ","
										//            << j_theta + 1 << "," << k_chi + 1 << "} = [ ";

										for (int b = 0; b < 4 * N_nu_; b = b + 4)
										{
											I = I_field_->block(i, j, k_start)[b_start + b];

											const std::string sep = b < (4 * N_nu_ - 1) ? "," : "";

											if (i_stokes == 0)
											{
												outputFile << std::scientific << std::setprecision(precision) << I << sep;
											}
											else
											{
												QUV = I_field_->block(i, j, k_start)[b_start + b + i_stokes];

												outputFile << std::scientific << std::setprecision(precision)
														   << 100.0 * QUV / I << sep;
											}
										}

										outputFile << "\n";
									}
								}
							}

							// Close the file
							outputFile.close();
						}
						else
						{
							if (mpi_rank_ == 0) std::cout << "\nERROR: failed to create the output file." << std::endl;
						}
					}
				}
			}
		}
	}

	if (mpi_rank_ == 0 and i_space == 0 and j_space == 0)
		std::cout << "Output written in " << file_name << "_i_j.csv" << "\n" << std::endl;
}

////////////////////////////////////////////////////////////////////////////
// write_surface_profiles
////////////////////////////////////////////////////////////////////////////
void 
RT_problem::write_surface_profiles(input_string file_name) const
{
	if (mpi_size_x_ > 1 or mpi_size_y_ > 1)
	{
		std::cerr << "\nWARNING: write_surface_profiles not supported for hotizontal decomposition!" << std::endl;
	}

	// indeces
	const int i_start = space_grid_->getGlobalStartX(); //g_dev.margin[0];
	const int j_start = space_grid_->getGlobalStartY(); //g_dev.margin[1];
	const int k_start = space_grid_->getGlobalStartZ(); //g_dev.margin[2];

	const int i_end = i_start + space_grid_->getLocalSizeX(); //g_dev.dim[0];
	const int j_end = j_start + space_grid_->getLocalSizeY(); //g_dev.dim[1];

	int i_global, j_global;

	double I, QUV;

	// write profiles
	if (space_grid_->local_to_global_coordinate(2, k_start) == 0)
	{
		// Create a new file
		input_string  output_file = file_name + ".m";
		std::ofstream outputFile(output_file);

		if (outputFile.is_open())
		{
			if (mpi_rank_ == 0) std::cout << "\nWriting surface output in " << output_file << "\n" << std::endl;

			// write grids
			outputFile << "\nnu_grid_ = [ ";
			for (int n = 0; n < N_nu_; ++n) outputFile << nu_grid_[n] << " ";
			outputFile << "];\n";

			outputFile << "\ntheta_grid = [ ";
			for (int j_theta = 0; j_theta < N_theta_; ++j_theta) outputFile << theta_grid_[j_theta] << " ";
			outputFile << "];\n";

			outputFile << "\nmu_grid = [ ";
			for (int j_theta = 0; j_theta < N_theta_; ++j_theta) outputFile << mu_grid_[j_theta] << " ";
			outputFile << "];\n";

			outputFile << "\nchi_grid = [ ";
			for (int k_chi = 0; k_chi < N_chi_; ++k_chi) outputFile << chi_grid_[k_chi] << " ";
			outputFile << "];\n";

			// create MATLAB data structure
			outputFile << "\nField = cell(4," << N_x_ << "," << N_y_ << "," << N_theta_ << "," << N_chi_ << ");"
					   << std::endl;
		}
		else
		{
			if (mpi_rank_ == 0) std::cout << "\nERROR: failed to create the output file." << std::endl;
		}

		// loop over spatial ppoints and directions
		for (int i = i_start; i < i_end; ++i)
		{
			i_global = space_grid_->local_to_global_coordinate(0, i);

			for (int j = j_start; j < j_end; ++j)
			{
				j_global = space_grid_->local_to_global_coordinate(1, j);

				for (int j_theta = N_theta_ / 2; j_theta < N_theta_; ++j_theta)
				{
					for (int k_chi = 0; k_chi < N_chi_; ++k_chi)
					{
						const int b_start = I_field_->local_to_block(j_theta, k_chi, 0);

						for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
						{
							outputFile << "\nField{" << i_stokes + 1 << "," << i_global + 1 << "," << j_global + 1 << ","
									   << j_theta + 1 << "," << k_chi + 1 << "} = [ ";

							for (int b = 0; b < 4 * N_nu_; b = b + 4)
							{
								I = I_field_->block(i, j, k_start)[b_start + b];

								if (i_stokes == 0)
								{
									outputFile << I << " ";
								}
								else
								{
									QUV = I_field_->block(i, j, k_start)[b_start + b + i_stokes];

									outputFile << 100.0 * QUV / I << " ";
								}
							}

							outputFile << "];\n";
						}
					}
				}
			}
		}

		// Close the file
		outputFile.close();
	}
}

////////////////////////////////////////////////////////////////////////////
// Write surface point profiles to file in CSV format
// write_surface_point_profiles_Omega_csv
////////////////////////////////////////////////////////////////////////////
void 
RT_problem::write_surface_point_profiles_Omega_csv(input_string file_name, const int i_space, const int j_space,
												   const unsigned int precision) const
{
	const int block_size = 4 * N_nu_;

	// indeces
	const int i_start = space_grid_->getGlobalStartX(); //g_dev.margin[0];
	const int j_start = space_grid_->getGlobalStartY(); //g_dev.margin[1];
	const int k_start = space_grid_->getGlobalStartZ(); //g_dev.margin[2];

	const int i_end = i_start + space_grid_->getLocalSizeX(); //g_dev.dim[0];
	const int j_end = j_start + space_grid_->getLocalSizeY(); //g_dev.dim[1];

	int i_global, j_global;

	double I, QUV;

	// write profiles
	if (space_grid_->local_to_global_coordinate(2, k_start) == 0)
	{
		for (int i = i_start; i < i_end; ++i)
		{
			i_global = space_grid_->local_to_global_coordinate(0, i);

			if (i_global == i_space)
			{
				for (int j = j_start; j < j_end; ++j)
				{
					j_global = space_grid_->local_to_global_coordinate(1, j);

					if (j_global == j_space)
					{
						// Create a new file
						input_string output_file =
							file_name + "_" + std::to_string(i_space) + "_" + std::to_string(j_space) + ".csv";
						std::ofstream outputFile(output_file);

						if (outputFile.is_open())
						{
							// Write CSV header
							outputFile << "I,Q/I(%),U/I(%),V/I(%)" << std::endl;

							// Write data rows
							for (int b = 0; b < block_size; b = b + 4)
							{
								I = I_field_Omega_->block(i, j, k_start)[b];

								outputFile << std::scientific << std::setprecision(precision) << I << ",";

								// Q/I, U/I, V/I
								for (int i_stokes = 1; i_stokes < 4; ++i_stokes)
								{
									QUV = I_field_Omega_->block(i, j, k_start)[b + i_stokes];

									std::string sep = (i_stokes < 3) ? "," : "";
									outputFile << std::scientific << std::setprecision(precision) << 100.0 * QUV / I
											   << sep;
								}

								outputFile << std::endl;
							}

							// Close the file
							outputFile.close();

							if (mpi_rank_ == 0) std::cout << "Output written in " << output_file << "\n" << std::endl;
						}
						else
						{
							if (mpi_rank_ == 0) std::cout << "\nERROR: failed to create the output file." << std::endl;
						}
					}
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////
// write surface profile in one single point
// write_surface_point_profiles_Omega
////////////////////////////////////////////////////////////////////////////
void 
RT_problem::write_surface_point_profiles_Omega(input_string file_name, const int i_space, const int j_space) const
{
	// // a single MPI rank writes output
	// if (mpi_rank_ == 0) std::cout << " Writing output in spatial point (" << i_space << ", " << j_space << ")" <<
	// std::endl;

	const int block_size = 4 * N_nu_;

	// indeces
	const int i_start = space_grid_->getGlobalStartX(); //g_dev.margin[0];
	const int j_start = space_grid_->getGlobalStartY(); //g_dev.margin[1];
	const int k_start = space_grid_->getGlobalStartZ(); //g_dev.margin[2];

	const int i_end = i_start + space_grid_->getLocalSizeX(); //g_dev.dim[0];
	const int j_end = j_start + space_grid_->getLocalSizeY(); //g_dev.dim[1];

	int i_global, j_global;

	double I, QUV;

	// write profiles
	if (space_grid_->local_to_global_coordinate(2, k_start) == 0)
	{
		for (int i = i_start; i < i_end; ++i)
		{
			i_global = space_grid_->local_to_global_coordinate(0, i);

			if (i_global == i_space)
			{
				for (int j = j_start; j < j_end; ++j)
				{
					j_global = space_grid_->local_to_global_coordinate(1, j);

					if (j_global == j_space)
					{
						// Create a new file
						input_string output_file =
							file_name + "_" + std::to_string(i_space) + "_" + std::to_string(j_space) + ".m";
						std::ofstream outputFile(output_file);

						if (outputFile.is_open())
						{
							// create MATLAB data structure
							outputFile << "\nField_Omega = cell(4,1);" << std::endl;

							for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
							{
								outputFile << "\nField_Omega{" << i_stokes + 1 << "} = [ ";

								for (int b = 0; b < block_size; b = b + 4)
								{
									I = I_field_Omega_->block(i, j, k_start)[b];

									if (i_stokes == 0)
									{
										outputFile << I << " ";
									}
									else
									{
										QUV = I_field_Omega_->block(i, j, k_start)[b + i_stokes];

										outputFile << 100.0 * QUV / I << " ";
									}
								}

								outputFile << "];\n";
							}

							// Close the file
							outputFile.close();

							if (mpi_rank_ == 0) std::cout << "Output written in " << output_file << "\n" << std::endl;
						}
						else
						{
							if (mpi_rank_ == 0) std::cout << "\nERROR: failed to create the output file." << std::endl;
						}
					}
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////
// write surface profile for arbitrary direction
// write_surface_profiles_Omega
////////////////////////////////////////////////////////////////////////////
void 
RT_problem::write_surface_profiles_Omega(input_string file_name) const
{
	if (mpi_size_x_ > 1 or mpi_size_y_ > 1)
	{
		std::cerr << "\nWARNING: write_surface_profiles_Omega not supported for hotizontal decomposition!" << std::endl;
	}

	const int block_size = 4 * N_nu_;

	// indeces
	const int i_start = space_grid_->getGlobalStartX(); //g_dev.margin[0];
	const int j_start = space_grid_->getGlobalStartY(); //g_dev.margin[1];
	const int k_start = space_grid_->getGlobalStartZ(); //g_dev.margin[2];

	const int i_end = i_start + space_grid_->getLocalSizeX(); //g_dev.dim[0];
	const int j_end = j_start + space_grid_->getLocalSizeY(); //g_dev.dim[1];

	int i_global, j_global;

	double I, QUV;

	// write profiles
	if (space_grid_->local_to_global_coordinate(2, k_start) == 0)
	{
		// Create a new file
		input_string  output_file = file_name + ".m";
		std::ofstream outputFile(output_file);

		if (outputFile.is_open())
		{
			if (mpi_rank_ == 0) std::cout << "\nWriting surface output in " << output_file << "\n" << std::endl;

			// create MATLAB data structure
			outputFile << "\nField = cell(4," << N_x_ << "," << N_y_ << ");" << std::endl;
		}
		else
		{
			if (mpi_rank_ == 0) std::cout << "\nERROR: failed to create the output file." << std::endl;
		}

		for (int i = i_start; i < i_end; ++i)
		{
			i_global = space_grid_->local_to_global_coordinate(0, i);

			for (int j = j_start; j < j_end; ++j)
			{
				j_global = space_grid_->local_to_global_coordinate(1, j);

				for (int i_stokes = 0; i_stokes < 4; ++i_stokes)
				{
					outputFile << "\nField{" << i_stokes + 1 << "," << i_global + 1 << "," << j_global + 1 << "} = [ ";

					for (int b = 0; b < block_size; b = b + 4)
					{
						I = I_field_Omega_->block(i, j, k_start)[b];

						if (i_stokes == 0)
						{
							outputFile << I << " ";
						}
						else
						{
							QUV = I_field_Omega_->block(i, j, k_start)[b + i_stokes];

							outputFile << 100.0 * QUV / I << " ";
						}
					}

					outputFile << "];\n";
				}
			}
		}

		// Close the file
		outputFile.close();
	}
}