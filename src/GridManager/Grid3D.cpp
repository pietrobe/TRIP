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
    Int local_owned = getLocalNumNodes();
    Int local_ghosted = getLocalNumNodesWithGhosts();

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
Int Grid3D::getNx() const { return Nx; }
Int Grid3D::getNy() const { return Ny; }
Int Grid3D::getNz() const { return Nz; }
Int Grid3D::getGlobalNumNodes() const { return Nr; } 
Int Grid3D::getLocalNumNodes() const { 
    return local_sizes[0] * local_sizes[1] * local_sizes[2]; 
} 
Int Grid3D::getLocalNumNodesWithGhosts() const {
    return getLocalSizeWithGhostX() * getLocalSizeWithGhostY() * getLocalSizeWithGhostZ(); 
}
Int Grid3D::getLocalSizeX() const { return local_sizes[0]; }
Int Grid3D::getLocalSizeY() const { return local_sizes[1]; }
Int Grid3D::getLocalSizeZ() const { return local_sizes[2]; }
std::array<Int,3> Grid3D::getLocalSizes() const { return local_sizes; }
Int Grid3D::getLocalSizeWithGhostX() const { return local_ghosted_sizes[0]; }
Int Grid3D::getLocalSizeWithGhostY() const { return local_ghosted_sizes[1]; }
Int Grid3D::getLocalSizeWithGhostZ() const { return local_ghosted_sizes[2]; }
std::array<Int,3> Grid3D::getLocalSizesWithGhost() const { return local_ghosted_sizes; }
Int Grid3D::getGlobalStartX() const { return global_start[0]; }
Int Grid3D::getGlobalStartY() const { return global_start[1]; }
Int Grid3D::getGlobalStartZ() const { return global_start[2]; }
std::array<Int,3> Grid3D::getGlobalStarts() const { return global_start; }
Int Grid3D::getGhostMarginX() const { return margin_start[0]; }
Int Grid3D::getGhostMarginY() const { return margin_start[1]; }
Int Grid3D::getGhostMarginZ() const { return margin_start[2]; }
std::array<Int,3> Grid3D::getGhostMargins() const { return margin_start; }
std::array<Int,3> Grid3D::getDecomposition() const { return {Px, Py, Pz}; }

Int Grid3D::local_to_global_coordinate(const int d, const Int local_coord) const { 
    return global_ghosted_start[d] + local_coord;
}

void Grid3D::create(
    std::array<Int, 3> p, 
    std::array<DMBoundaryType,3> boundary_type,
    Int ghost_width)
{
    // create DMDA
    PetscErrorCode ierr = DMDACreate3d(
        comm,
        boundary_type[0], boundary_type[1], boundary_type[2], 
        DMDA_STENCIL_BOX,
        Nx, Ny, Nz,
        p[0], p[1], p[2], // MPI decomposition
        /*dof*/ 1, 
        ghost_width,
        NULL, NULL, NULL, // uniform partitioning
        &da
    );
    if (ierr) throw std::runtime_error("DMDACreate3d failed in Grid3D::create.");
    ierr = DMSetUp(da);
    if (ierr) throw std::runtime_error("DMSetUp failed in Grid3D::create.");

    // retrieve actual decomposition chosen by petsc
    ierr = DMDAGetInfo(
        da, NULL, NULL, NULL, NULL,
        &Px, &Py, &Pz,
        NULL, NULL, NULL, NULL, NULL, NULL
    );
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