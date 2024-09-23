#ifndef __RT_UTILITY_HPP__
#define __RT_UTILITY_HPP__

#include <string>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <filesystem>

// Function declarations
std::string getCurrentDateTime();

std::map<std::string, std::string> get_input_PORTA(const std::filesystem::path &config_file, int mpi_rank);

std::string getOptionArgument(int argc, char *argv[], const std::string &option);

bool getOptionFlag(int argc, char *argv[], const std::string &option);

void print_help();

#endif