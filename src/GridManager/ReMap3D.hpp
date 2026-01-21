#pragma once
#include <petsc.h>
#include <vector>
#include <array>
#include <cassert>
#include <iostream>
#include <omp.h>
#include <iomanip>
#include <memory>
#include <fstream>
#include <type_traits>
#include "Grid3D.hpp"
#include "Field.hpp"

using Field_ptr = std::shared_ptr<Field>;

/** @brief
 * This ReMap3D class is based on the homonymous class in sgrid library
 * (https://bitbucket.org/zulianp/sgrid/src/main/core/sgrid_ReMap.hpp).
 * Parallel: [space (distributed)]×[full block]
 * Serial: [full space]×[block slice (distributed)]
 */
class ReMap3D {
public:
    using Real = double;

    bool enable_timing = false; // default: false

    /** @brief default constructor for ReMap3D.
     * @note This requires a call to the init method afterwards.
     */
    ReMap3D() = default;

    /** @brief constructor for ReMap3D that takes as input the following parameter:
     * @param space_grid spatially-distributed grid
     * @param block_grid block-distributed grid
     * @param block_size_ size of the block of the fields defined in space_grid
     * @param tile_size_ size of the block of the fields defined in block_grid
     * @param num_tiles_per_block_ (optional) number of pieces in which each block is split. By default, 
     *          it is set as block_size_ / tile_size_.
     * @note The same ReMap3D object can be used for all spatially-distributed fields and block-distributed fields
     *          having the same layout as the one declared in this constructor.
     */
    ReMap3D(
        Grid_ptr space_grid_,
        Grid_ptr block_grid_,
        int block_size_,
        int tile_size_,
        int num_tiles_per_block_ = -1) 
        : space_grid(space_grid_),
          block_grid(block_grid_),
          comm(space_grid_->getComm()),
          block_size(block_size_),
          tile_size(tile_size_)
    {
        init(block_size_, tile_size_, num_tiles_per_block_);
    }

    /** @brief convert from spatially-distributed Field to block-distributed Field (serial is on rank 0).
     * @param space_distr spatially-distributed field
     * @param block_distr block-distributed field
     * @param num_tiles number of tiles in which the redistribution happens
     */
    void from_space_to_block_distributed(Field_ptr space_distr, Field_ptr block_distr, int num_tiles = 0);


    /** @brief convert from block-distributed field to spatially-distributed field. */
    void from_block_to_space_distributed(Field_ptr block_distr, Field_ptr space_distr, PetscInt num_tiles = 0);


    /** @brief initialise the values using the space-distributed and block-distributed grid.
     * @param space_grid spatially-distributed grid
     * @param block_grid block-distributed grid
     * @param block_size size of the block of the fields defined in space_grid
     * @param tile_size size of the block of the fields defined in block_grid
     * @param num_tiles_per_block_ number of pieces in which each block is split
     * If num_tiles_per_block_ == -1, it is set as space_field.getBlockSize() / block_field.getBlockSize().
     */
    void init(
        Grid_ptr space_grid_, Grid_ptr block_grid_,
        PetscInt block_size_, PetscInt tile_size_,
        PetscInt num_tiles_per_block_ = -1);
    
private:
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0, size = 1;

    Grid_ptr space_grid;
    Grid_ptr block_grid;

    PetscInt local_nodes = 0; // local size for rank rank
    PetscInt total_global_nodes = 0; // global number of grid nodes
    PetscInt tile_size = 0;
    PetscInt block_size = 0;
    PetscInt num_tiles_per_block = -1;

    std::vector<std::array<PetscInt,3>> local_starts = {}; // start indices of each rank
    std::vector<std::array<PetscInt,3>> local_sizes = {}; // local sizes of each rank
    std::vector<PetscInt> per_rank_nodes = {}; // number of nodes per rank

    bool uniform = true;
    std::vector<PetscMPIInt> send_counts = {}, send_offsets = {};
    std::vector<PetscMPIInt> recv_counts = {}, recv_offsets = {};

    std::vector<Real> space_buffer = {};   
    std::vector<Real> block_buffer = {}; 

    /** @brief initialise the values using the space-distributed and block-distributed grid.
     * @param block_size size of the block of the fields defined in space_grid
     * @param tile_size size of the block of the fields defined in block_grid
     * @param num_tiles_per_block_ number of pieces in which each block is split
     * If num_tiles_per_block_ == -1, it is set as space_field.getBlockSize() / block_field.getBlockSize().
     */
    void init(
        PetscInt block_size_, PetscInt tile_size_,
        PetscInt num_tiles_per_block_);

    /** @brief pack the space-distributed field into the space_buffer communicator buffer. */ 
    void pack_space_distr_field(Field_ptr field, PetscInt num_tiles);

    /** @brief unpack the space_buffer communicator buffer the space-distributed field. */ 
    void unpack_space_distr_field(Field_ptr space_field, PetscInt num_tiles);

    /** @brief pack the local block-distributed field into the received block_buffer communicator buffer. */ 
    void pack_block_distr_field(Field_ptr field);

    /** @brief unpack the received block_buffer communicator buffer into the local block-distributed field. */ 
    void unpack_block_distr_field(Field_ptr field);

    template<class T>
    MPI_Datatype mpi_type() const {
        if (std::is_same<T, double>::value) return MPI_DOUBLE;
        if (std::is_same<T, float>::value) return MPI_FLOAT;
        if (std::is_same<T, int>::value) return MPI_INT;
        if (std::is_same<T, long>::value) return MPI_LONG;
        if (std::is_same<T, PetscInt>::value) return MPIU_INT;
        return MPI_BYTE; 
    }

};
