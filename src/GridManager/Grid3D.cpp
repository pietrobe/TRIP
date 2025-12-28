#include <petsc.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <array>
#include <string>
#include <stdexcept>
#include <vector>
#include <iostream>
#include "Grid3D.hpp"

void Grid3D::print_info()
{
    int size = local_sizes[0] * local_sizes[1] * local_sizes[2];
    std::cout <<  "[Rank " << rank << "] " << name << " has " << size << " cells out of " << Nr   
                << " (start_indices = {" << start_indices[0] << ", " << start_indices[1] << ", " << start_indices[2] 
                << "}, local_indices = {" << local_sizes[0] << ", " << local_sizes[1] << ", " << local_sizes[2] << "})\n" << std::flush; 
}


DM Grid3D::getDA() const { return da; }
std::string Grid3D::getName() { return name; }
MPI_Comm Grid3D::getComm() const { return comm; }
int Grid3D::getRank() const { return rank; }
int Grid3D::getMPISize() const { return mpi_size; }
PetscInt Grid3D::getNx() const { return Nx; }
PetscInt Grid3D::getNy() const { return Ny; }
PetscInt Grid3D::getNz() const { return Nz; }
PetscInt Grid3D::getGlobalNumNodes() const { return Nr; } 
PetscInt Grid3D::getLocalNumNodes() const { return local_sizes[0] * local_sizes[1] * local_sizes[2]; } 
PetscInt Grid3D::getLocalSizeX() const { return local_sizes[0]; }
PetscInt Grid3D::getLocalSizeY() const { return local_sizes[1]; }
PetscInt Grid3D::getLocalSizeZ() const { return local_sizes[2]; }
std::array<PetscInt,3> Grid3D::getLocalSizes() const { return local_sizes; }
PetscInt Grid3D::getGlobalStartX() const { return start_indices[0]; }
PetscInt Grid3D::getGlobalStartY() const { return start_indices[1]; }
PetscInt Grid3D::getGlobalStartZ() const { return start_indices[2]; }
std::array<PetscInt,3> Grid3D::getGlobalStarts() const { return start_indices; }

PetscInt Grid3D::local_to_global_coordinate(const int d, const PetscInt local_coord) const { 
    return start_indices[d] + local_coord; // no margin since we do not have ghost layers
}

void Grid3D::create(
    std::array<PetscInt, 3> p, 
    std::array<DMBoundaryType,3> boundary_type)
{
    Px = p[0]; Py = p[1]; Pz = p[2];
    // create DMDA
    PetscErrorCode ierr = DMDACreate3d(
        comm,
        boundary_type[0], boundary_type[1], boundary_type[2], 
        DMDA_STENCIL_BOX,
        Nx, Ny, Nz,
        Px, Py, Pz, // MPI decomposition
        /*dof*/ 1, 
        0,
        NULL, NULL, NULL, // uniform partitioning
        &da
    );
    if (ierr) throw std::runtime_error("DMDACreate3d failed in Grid3D::create.");
    ierr = DMSetUp(da);
    if (ierr) throw std::runtime_error("DMSetUp failed in Grid3D::create.");
    // get local sizes
    ierr = DMDAGetCorners(
        da, 
        &start_indices[0], &start_indices[1], &start_indices[2], 
        &local_sizes[0], &local_sizes[1], &local_sizes[2]
    );
    if (ierr) throw std::runtime_error("Failed DMDAGetCorners in Grid3D::create.");
}