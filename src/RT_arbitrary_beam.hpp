#ifndef RT_ARBITRARY_BEAM_HPP
#define RT_ARBITRARY_BEAM_HPP

#include <filesystem>
#include <memory>
#include <string>
#include <vector>
#include "RT_problem.hpp"
#include "RT_solver.hpp"
#include "RT_utility.hpp" // for BeamDirection and getCurrentDateTime

// Compute a single arbitrary beam
void
compute_arbitrary_beam(RT_solver				   &rt_solver,		//
					   std::shared_ptr<RT_problem> &rt_problem_ptr, //
					   Real							mu,				//
					   Real							chi,			//
					   const std::string		   &output_file);				//

// Compute a single arbitrary beam
void
compute_arbitrary_beam_hdf(RT_solver						&rt_solver,			//
						   bool								 write_text_output, //
						   std::shared_ptr<RT_problem>		&rt_problem_ptr,	//
						   const std::vector<BeamDirection> &beams,				//
						   const int						 beam_index,		//
						   const std::string				&output_file);						//

// Format arbitrary beams info as a string for printing//
std::string
format_arbitrary_beams_info(const std::vector<BeamDirection> &beams); //

// Process all arbitrary beams and return the number of beams processed
int
process_arbitrary_beams(const std::vector<BeamDirection> &beams,				  //
						RT_solver						 &rt_solver,			  //
						std::shared_ptr<RT_problem>		 &rt_problem_ptr,		  //
						const std::string				 &output_file,			  //
						const std::filesystem::path		 &output_info_file_name = {}); //

int
process_arbitrary_beams_hdf(const std::vector<BeamDirection> &beams,				  //
							bool							  write_text_output,	  //
							RT_solver						 &rt_solver,			  //
							std::shared_ptr<RT_problem>		 &rt_problem_ptr,		  //
							const std::string				 &output_file,			  //
							const std::filesystem::path		 &output_info_file_name = {}); //

#endif // RT_ARBITRARY_BEAM_HPP