#include "RT_arbitrary_beam.hpp"
#include <mpi.h>	 // for MPI_Wtime
#include <algorithm> // for std::remove
#include <cmath>	 // for std::acos
#include <cstdio>	 // for std::snprintf
#include <filesystem>
#include <fstream>	// for std::ofstream
#include <iomanip>	// for std::setw, std::setprecision
#include <iostream> // for std::cout
#include <sstream>	// for std::stringstream

////////////////////////////////////////////////////////
// Compute a single arbitrary beam
// compute_arbitrary_beam
////////////////////////////////////////////////////////
void
compute_arbitrary_beam(RT_solver				   &rt_solver,		//
					   std::shared_ptr<RT_problem> &rt_problem_ptr, //
					   Real							mu,				//
					   Real							chi,			//
					   const std::string		   &output_file)
{
	char mu_charv[40];
	char chi_charv[40];

	std::snprintf(&mu_charv[0], 40, "%.4f", mu);
	std::snprintf(&chi_charv[0], 40, "%.4f", chi);

	std::string mu_str(mu_charv);
	std::string chi_str(chi_charv);

	// remove the dot from the output file name
	mu_str.erase(std::remove(mu_str.begin(), mu_str.end(), '.'), mu_str.end());
	chi_str.erase(std::remove(chi_str.begin(), chi_str.end(), '.'), chi_str.end());

	const std::string output_file_Omega_mu = output_file + "_mu" + mu_str + "_chi" + chi_str;

	const Real theta = std::acos(mu);
	rt_solver.apply_formal_solver_Omega(theta, chi);

	const int Nx = rt_problem_ptr->N_x_;
	const int Ny = rt_problem_ptr->N_y_;

	for (int i = 0; i < Nx; ++i)
	{
		for (int j = 0; j < Ny; ++j)
		{
			rt_problem_ptr->write_surface_point_profiles_Omega(output_file_Omega_mu, i, j);
			rt_problem_ptr->write_surface_point_profiles_Omega_csv(output_file_Omega_mu, i, j, 14);
		}
	}
}

////////////////////////////////////////////////////////
// Compute a single arbitrary beam
// compute_arbitrary_beam_hdf
////////////////////////////////////////////////////////
void
compute_arbitrary_beam_hdf(RT_solver						&rt_solver,			//
						   bool								 write_text_output, //
						   std::shared_ptr<RT_problem>		&rt_problem_ptr,	//
						   const std::vector<BeamDirection> &beams,				//
						   const int						 beam_index,		//
						   const std::string				&output_file)						//
{
	const Real mu  = beams[beam_index].mu;
	const Real chi = beams[beam_index].chi;

	const Real theta = std::acos(mu);
	rt_solver.apply_formal_solver_Omega(theta, chi);

	// Accumulate surface profiles data across MPI ranks
	std::vector<Real> surface_data_I;
	std::vector<Real> surface_data_Q;
	std::vector<Real> surface_data_U;
	std::vector<Real> surface_data_V;

	MPI_Comm write_comm;
	const auto [color, mpi_status] = rt_problem_ptr->make_write_surface_MPI_Comm(MPI_COMM_WORLD, write_comm);

	if (color == 0)
	{
		MPI_Barrier(write_comm);
		MPI_Comm_free(&write_comm);
		return;
	}

	rt_problem_ptr->accumulate_surface_profiles_Omega_domain_data(surface_data_I,  //
																  surface_data_Q,  //
																  surface_data_U,  //
																  surface_data_V); //

	rt_problem_ptr->write_emergent_field_Omega_hdf5(output_file,	   //
													write_comm,		   //
													beams, beam_index, //
													surface_data_I,	   //
													surface_data_Q,	   //
													surface_data_U,	   //
													surface_data_V);   //

	MPI_Barrier(write_comm);
	MPI_Comm_free(&write_comm);

	if (write_text_output)
	{ // Write individual surface profiles too
		// For TESTING and comparison with old output

		char mu_charv[40];
		char chi_charv[40];

		std::snprintf(&mu_charv[0], 40, "%.4f", mu);
		std::snprintf(&chi_charv[0], 40, "%.4f", chi);

		std::string mu_str(mu_charv);
		std::string chi_str(chi_charv);

		// remove the dot from the output file name
		mu_str.erase(std::remove(mu_str.begin(), mu_str.end(), '.'), mu_str.end());
		chi_str.erase(std::remove(chi_str.begin(), chi_str.end(), '.'), chi_str.end());

		const std::filesystem::path output_file_fs = output_file;
		const std::filesystem::path head		   = output_file_fs.parent_path();
		const std::string			file_name_base = "emergent_field_abitrary_Omega";

		const std::string output_file_Omega_mu = (head / (file_name_base + "_mu" + mu_str + "_chi" + chi_str)).string();

		const int Nx = rt_problem_ptr->N_x_;
		const int Ny = rt_problem_ptr->N_y_;

		for (int i = 0; i < Nx; ++i)
		{
			for (int j = 0; j < Ny; ++j)
			{
				rt_problem_ptr->write_surface_point_profiles_Omega(output_file_Omega_mu, i, j);
				rt_problem_ptr->write_surface_point_profiles_Omega_csv(output_file_Omega_mu, i, j, 14);
			}
		}
	} // END if (write_text_output)
}

////////////////////////////////////////////////////////
// Format arbitrary beams info string for printing
// format_arbitrary_beams_info
////////////////////////////////////////////////////////
std::string
format_arbitrary_beams_info(const std::vector<BeamDirection> &beams)
{
	std::stringstream ss;
	ss << "\n" << std::string(50, '=') << std::endl;
	ss << "Arbitrary Beam Directions" << std::endl;
	ss << std::string(50, '=') << std::endl;
	ss << std::setw(8) << "Beam #" << std::setw(15) << "mu" << std::setw(15) << "chi" << std::setw(12) << "theta (rad)"
	   << std::endl;
	ss << std::string(50, '-') << std::endl;

	int beam_count = 0;
	for (const auto &beam : beams)
	{
		Real theta = std::acos(beam.mu);
		ss << std::setw(8) << beam_count++ << std::setw(15) << std::fixed << std::setprecision(6) << beam.mu
		   << std::setw(15) << std::fixed << std::setprecision(6) << beam.chi << std::setw(12) << std::fixed
		   << std::setprecision(6) << theta << std::endl;
	}
	ss << std::string(50, '-') << std::endl;
	ss << "Total number of beams: " << beam_count << std::endl;
	ss << std::string(50, '=') << std::endl << std::endl;

	return ss.str();
}

////////////////////////////////////////////////////////
// Process arbitrary beams from config
// process_arbitrary_beams
////////////////////////////////////////////////////////
int
process_arbitrary_beams(const std::vector<BeamDirection> &beams, RT_solver &rt_solver,
						std::shared_ptr<RT_problem> &rt_problem_ptr, const std::string &output_file,
						const std::filesystem::path &output_info_file_name) // No default value here
{
	const int mpi_rank = rt_problem_ptr->mpi_rank_;

	if (beams.empty())
	{
		if (mpi_rank == 0) std::cout << "WARNING: no arbitrary beams defined in config" << std::endl;
		return 0;
	}

	// Print beam info
	if (mpi_rank == 0)
	{
		std::string beam_info = format_arbitrary_beams_info(beams);
		std::cout << beam_info;

		if (!output_info_file_name.empty())
		{
			std::ofstream output_file_info(output_info_file_name, std::ios::app);
			output_file_info << beam_info;
			output_file_info.close();
		}
	}

	std::cout.flush();

	const double tick		   = MPI_Wtime();
	int			 beams_counter = 0;

	for (const auto &beam : beams)
	{
		if (mpi_rank == 0)
		{
			std::cout << getCurrentDateTime() << " - Computing arbitrary beam: mu = " << beam.mu << ", chi = " << beam.chi
					  << std::endl;
		}

		compute_arbitrary_beam(rt_solver, rt_problem_ptr, beam.mu, beam.chi, output_file);
		beams_counter++;
	}

	const double tock = MPI_Wtime();

	if (mpi_rank == 0)
	{
		std::cout << "Arbitrary beams time (s) = " << tock - tick << std::endl;
		std::cout << "Number of beams          = " << beams_counter << std::endl;
		std::cout << "Mean time per beam (s)   = " << (tock - tick) / double(beams_counter) << std::endl;
	}

	return beams_counter;
}

int
process_arbitrary_beams_hdf(const std::vector<BeamDirection> &beams,			 //
							bool							  write_text_output, //
							RT_solver						 &rt_solver,		 //
							std::shared_ptr<RT_problem>		 &rt_problem_ptr,	 //
							const std::string				 &output_file,		 //
							const std::filesystem::path		 &output_info_file_name)	 //
{
	const int mpi_rank = rt_problem_ptr->mpi_rank_;

	if (beams.empty())
	{
		if (mpi_rank == 0) std::cout << "WARNING: no arbitrary beams defined in config" << std::endl;
		return 0;
	}

	// Print beam info
	if (mpi_rank == 0)
	{
		std::string beam_info = format_arbitrary_beams_info(beams);
		std::cout << beam_info;

		if (!output_info_file_name.empty())
		{
			std::ofstream output_file_info(output_info_file_name, std::ios::app);
			output_file_info << beam_info;
			output_file_info.close();
		}
	}

	std::cout.flush();

	const double tick		   = MPI_Wtime();
	int			 beams_counter = 0;

	rt_problem_ptr->write_beams_frequency_grids_Omega_hdf5(beams, output_file);

	for (int beam_index = 0; beam_index < beams.size(); ++beam_index)
	{
		const auto &beam = beams[beam_index];
		if (mpi_rank == 0)
		{
			std::cout << getCurrentDateTime() << " - Computing arbitrary beam: mu = " << beam.mu << ", chi = " << beam.chi
					  << std::endl;
		}

		compute_arbitrary_beam_hdf(rt_solver, write_text_output, rt_problem_ptr, beams, beam_index, output_file);
		beams_counter++;
	}

	const double tock = MPI_Wtime();

	if (mpi_rank == 0)
	{
		std::cout << "Arbitrary beams time (s) = " << tock - tick << std::endl;
		std::cout << "Number of beams          = " << beams_counter << std::endl;
		std::cout << "Mean time per beam (s)   = " << (tock - tick) / double(beams_counter) << std::endl;
	}

	return beams_counter;
}
