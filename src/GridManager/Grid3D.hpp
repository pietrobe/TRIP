#pragma once
#include <petsc.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <array>
#include <string>
#include <stdexcept>
#include <vector>
#include <iostream>


/** @brief wrapper around PETSc DMDA describing a 3D spatial grid. */
class Grid3D {

public:
    /** @brief constructors that takes the following inputs to create a distributed PETSc DMDA 3D grid.
     * @param nx, ny, nz spatial grid dimensions.
     * @param p (optional) array of 3 integer, the i-th int is the number of MPI process in the i-th direction (all set to PETSC_DECIDE by default).
     * @param boundary_type (optional) array of 3 PETSc boundary type (https://petsc.org/main/manualpages/DM/DMBoundaryType/), 
     *          one for each dimension.  By default, x and y are period and z has no ghosts nodes.
     * @param name_ (optional) name of the grid
     */
    Grid3D(
        MPI_Comm c, 
        PetscInt nx, PetscInt ny, PetscInt nz,
        std::array<PetscInt, 3> p = {PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE}, 
        std::array<DMBoundaryType,3> boundary_type = {DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_NONE},
        std::string name_ = "Grid")
        : comm(c), Nx(nx), Ny(ny), Nz(nz)
    {
        Nr = nx * ny * nz;
        MPI_Comm_rank(c, &rank);
        MPI_Comm_size(c, &mpi_size);
        if (name_.compare("Grid") != 0) name = "Grid " + name_;
        create(p, boundary_type);
    }

    /** @brief destructor that destroys the DMDA object if it exists. */
    ~Grid3D() { if (da) (void)DMDestroy(&da); }

    void print_info();

    // --------- getters ------------
    DM getDA() const;
    std::string getName();
    MPI_Comm getComm() const;
    int getRank() const;
    int getMPISize() const;
    PetscInt getNx() const;
    PetscInt getNy() const;
    PetscInt getNz() const;
    PetscInt getGlobalNumNodes() const;
    PetscInt getLocalNumNodes() const;
    PetscInt getLocalSizeX() const;
    PetscInt getLocalSizeY() const;
    PetscInt getLocalSizeZ() const;
    std::array<PetscInt,3> getLocalSizes() const;
    PetscInt getGlobalStartX() const;
    PetscInt getGlobalStartY() const;
    PetscInt getGlobalStartZ() const;
    std::array<PetscInt,3> getGlobalStarts() const;

    /** @brief method that given a dimension and the local coordinate of a cell, returns the global coordinate of the 
     * cell in that direction.
     * @param d the chosen dimension.
     * @param local_coord the local index of the cell.
     */
    PetscInt local_to_global_coordinate(const int d, const PetscInt local_coord) const;

    /** @brief iterates over the local domain owned by this process and 
     * applies the passed lambda function f. 
     */
    template<typename Function>
    void parallel_for(Function&& f) const 
    {
        #pragma omp parallel for collapse(3)
        for (PetscInt k = 0; k < local_sizes[2]; k++) {
            for (PetscInt j = 0; j < local_sizes[1]; j++) {
                for (PetscInt i = 0; i < local_sizes[0]; i++) {
                    f(i,j,k);
                }
            }
        }
    }

    //TODO iterator for (i,j,k)

private:
    DM da = nullptr; // PETSc DMDA grid (only topology)
    PetscInt Nx = 0, Ny = 0, Nz = 0, Nr = 0; // Spatial grid dimensions
    PetscInt Px = PETSC_DECIDE, Py = PETSC_DECIDE, Pz = PETSC_DECIDE; // number of processors in each direction
    
    MPI_Comm comm = PETSC_COMM_WORLD; // MPI communicator 
    int rank, mpi_size;
    std::array<PetscInt,3> local_sizes = {0,0,0}; // local sizes owned by this MPI process
    std::array<PetscInt,3> start_indices = {0,0,0};    // Global start indices of this process
    std::string name = "Grid";

    /** @brief create the 3D DMDA grid. */
    void create(
        std::array<PetscInt, 3> p, 
        std::array<DMBoundaryType,3> boundary_type);

};
