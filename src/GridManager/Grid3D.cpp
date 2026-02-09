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
    PetscInt local_owned = getLocalNumNodes();
    PetscInt local_ghosted = getLocalNumNodesWithGhosts();

    std::cout << "[Rank " << rank << "] " << name << ":\n"
              << "  Global size       = (" << Nx << ", " << Ny << ", " << Nz << ") = " << Nr << " nodes\n"
              << "  Local owned size  = (" << local_sizes[0] << ", " << local_sizes[1] << ", " << local_sizes[2] 
              << ") = " << local_owned << " nodes\n"
              << "  Local ghosted size= (" << local_ghosted_sizes[0] << ", " << local_ghosted_sizes[1] << ", " << local_ghosted_sizes[2] 
              << ") = " << local_ghosted << " nodes\n"
              << "  Ghost margins     = start: (" << margin_start[0] << ", " << margin_start[1] << ", " << margin_start[2] 
              << "), end: (" << margin_end[0] << ", " << margin_end[1] << ", " << margin_end[2] << ")\n"
              << "  Global start idx  = owned: (" << global_start[0] << ", " << global_start[1] << ", " << global_start[2] 
              << "), ghosted: (" << global_ghosted_start[0] << ", " << global_ghosted_start[1] << ", " << global_ghosted_start[2] << ")\n"
              << std::flush;
}


DM Grid3D::getDA() const { return da; }
std::string Grid3D::getName() { return name; }
MPI_Comm Grid3D::getComm() const { return comm; }
int Grid3D::getRank() const { return rank; }
int Grid3D::getMPISize() const { return mpi_size; }
int Grid3D::getNx() const { return Nx; }
int Grid3D::getNy() const { return Ny; }
int Grid3D::getNz() const { return Nz; }
PetscInt Grid3D::getGlobalNumNodes() const { return Nr; } 
PetscInt Grid3D::getLocalNumNodes() const { 
    return local_sizes[0] * local_sizes[1] * local_sizes[2]; 
} 
PetscInt Grid3D::getLocalNumNodesWithGhosts() const {
    return getLocalSizeWithGhostX() * getLocalSizeWithGhostY() * getLocalSizeWithGhostZ(); 
}
PetscInt Grid3D::getLocalSizeX() const { return local_sizes[0]; }
PetscInt Grid3D::getLocalSizeY() const { return local_sizes[1]; }
PetscInt Grid3D::getLocalSizeZ() const { return local_sizes[2]; }
std::array<PetscInt,3> Grid3D::getLocalSizes() const { return local_sizes; }
PetscInt Grid3D::getLocalSizeWithGhostX() const { return local_ghosted_sizes[0]; }
PetscInt Grid3D::getLocalSizeWithGhostY() const { return local_ghosted_sizes[1]; }
PetscInt Grid3D::getLocalSizeWithGhostZ() const { return local_ghosted_sizes[2]; }
std::array<PetscInt,3> Grid3D::getLocalSizesWithGhost() const { return local_ghosted_sizes; }
PetscInt Grid3D::getGlobalStartX() const { return global_start[0]; }
PetscInt Grid3D::getGlobalStartY() const { return global_start[1]; }
PetscInt Grid3D::getGlobalStartZ() const { return global_start[2]; }
std::array<PetscInt,3> Grid3D::getGlobalStarts() const { return global_start; }
PetscInt Grid3D::getGhostMarginX() const { return margin_start[0]; }
PetscInt Grid3D::getGhostMarginY() const { return margin_start[1]; }
PetscInt Grid3D::getGhostMarginZ() const { return margin_start[2]; }
std::array<PetscInt,3> Grid3D::getGhostMargins() const { return margin_start; }

PetscInt Grid3D::local_to_global_coordinate(const int d, const PetscInt local_coord) const { 
    return global_ghosted_start[d] + local_coord;
}

void Grid3D::create(
    std::array<PetscInt, 3> p, 
    std::array<DMBoundaryType,3> boundary_type,
    int ghost_width)
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
        ghost_width,
        NULL, NULL, NULL, // uniform partitioning
        &da
    );
    if (ierr) throw std::runtime_error("DMDACreate3d failed in Grid3D::create.");
    ierr = DMSetUp(da);
    if (ierr) throw std::runtime_error("DMSetUp failed in Grid3D::create.");
    // get local sizes without ghost layers
    ierr = DMDAGetCorners(
        da, 
        &global_start[0], &global_start[1], &global_start[2], 
        &local_sizes[0], &local_sizes[1], &local_sizes[2]
    );
    ierr = DMDAGetGhostCorners(
        da, 
        &global_ghosted_start[0], &global_ghosted_start[1], &global_ghosted_start[2], 
        &local_ghosted_sizes[0], &local_ghosted_sizes[1], &local_ghosted_sizes[2]
    );
    for (int d = 0; d < 3; ++d) {
        // ghost cells before owned region
        margin_start[d] = global_start[d] - global_ghosted_start[d];

        // ghost cells after owned region
        margin_end[d] =
            (global_ghosted_start[d] + local_ghosted_sizes[d]) -
            (global_start[d] + local_sizes[d]);
    }
    if (ierr) throw std::runtime_error("Failed DMDAGetCorners in Grid3D::create.");
}