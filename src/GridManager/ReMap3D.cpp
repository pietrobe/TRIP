#include <petsc.h>
#include <vector>
#include <array>
#include <cassert>
#include <iostream>
#include <omp.h>
#include <iomanip>
#include <fstream>
#include <type_traits>
#include "Grid3D.hpp"
#include "Field.hpp"
#include "ReMap3D.hpp"

/** MPI utils */
void check_mpi(int err) {
    if (err != MPI_SUCCESS) {
        char msg[MPI_MAX_ERROR_STRING];
        int len = 0;
        MPI_Error_string(err, msg, &len);
        std::cerr << "MPI error: " << std::string(msg, len) << std::endl;
        MPI_Abort(MPI_COMM_WORLD, err);
    }
}

void ReMap3D::from_space_to_block_distributed(Field_ptr space_distr, Field_ptr block_distr, int num_tiles)
{
    if (space_distr->getBlockSize() != block_size) 
        throw std::runtime_error(
            "In ReMap3D::from_space_to_block_distributed, provided spatially-distributed field does not have the right block size! (" 
            + std::to_string(block_size) + "-" 
            + std::to_string(space_distr->getBlockSize()) + ")"
        );
    if (block_distr->getBlockSize() != tile_size) 
        throw std::runtime_error(
            "In ReMap3D::from_space_to_block_distributed, provided block-distributed field does not have the right tile size! (" 
            + std::to_string(tile_size) + "-" 
            + std::to_string(block_distr->getBlockSize()) + ")"
        );

    // pack local data into communication buffer
    double t_pack = 0.0, t_mpi = 0.0, t_unpack = 0.0;
    if(enable_timing) t_pack = -MPI_Wtime();
    pack_space_distr_field(space_distr, num_tiles);
    if(enable_timing) t_pack += MPI_Wtime();

    // MPI all to all
    if(enable_timing) t_mpi = -MPI_Wtime();
    if (uniform) {
        int sendcount = space_buffer.size() / size;
        int recvcount = block_buffer.size() / size;
        check_mpi(
            MPI_Alltoall(
                space_buffer.data(), sendcount, mpi_type<Real>(),
                block_buffer.data(), recvcount, mpi_type<Real>(),
                comm
            )
        );
    } else {
        check_mpi(
            MPI_Alltoallv(
                space_buffer.data(), send_counts.data(), send_offsets.data(), mpi_type<Real>(),
                block_buffer.data(), recv_counts.data(), recv_offsets.data(), mpi_type<Real>(),
                comm
            )
        );
    }
    if(enable_timing) t_mpi += MPI_Wtime();

    // unpack into local block-distributed field 
    if(enable_timing) t_unpack = -MPI_Wtime();
    unpack_block_distr_field(block_distr);
    if(enable_timing) {
        t_unpack += MPI_Wtime();
        double max_mpi, max_pack, max_unpack;
        MPI_Reduce(&t_mpi, &max_mpi, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(&t_pack, &max_pack, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(&t_unpack, &max_unpack, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        double bytes_sent = static_cast<double>(space_buffer.size() * sizeof(Real)); 
        double bandwidth = 0.0;
        MPI_Reduce(&bytes_sent, &bandwidth, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (rank == 0) {
            bandwidth /= (max_mpi * 1.0e6); // MB/s 
            std::cout << "[ReMap3D::from_block_to_space_distributed] Max: pack = " << max_pack << " s, MPI = " << max_mpi 
                << " s, unpack = " << max_unpack << " s, total = " << (max_pack + max_mpi + max_unpack) << " s. "
                << "MPI aggregate effective BW = " << bandwidth << " MB/s\n" << std::flush;
        }
    }
}


void ReMap3D::from_block_to_space_distributed(Field_ptr block_distr, Field_ptr space_distr, PetscInt num_tiles)
{
    if (space_distr->getBlockSize() != block_size) 
        throw std::runtime_error("In ReMap3D::from_block_to_space_distributed, provided spatially-distributed field does not have the right block size!");
    if (block_distr->getBlockSize() != tile_size) 
        throw std::runtime_error("In ReMap3D::from_block_to_space_distributed, provided block-distributed field does not have the right tile size!");
        
    // pack local data into communication buffer
    double t_pack = 0.0, t_mpi = 0.0, t_unpack = 0.0;
    if(enable_timing) t_pack = -MPI_Wtime();
    pack_block_distr_field(block_distr);
    if(enable_timing) t_pack += MPI_Wtime();

    // MPI all to all 
    if(enable_timing) t_mpi = -MPI_Wtime();
    if (uniform) {
        int sendcount = static_cast<PetscInt>(block_buffer.size() / size);
        int recvcount = static_cast<PetscInt>(space_buffer.size() / size);
        check_mpi(
            MPI_Alltoall(
                block_buffer.data(), sendcount, mpi_type<Real>(),
                space_buffer.data(), recvcount, mpi_type<Real>(),
                comm
            )
        );
    } else {
        check_mpi(
            MPI_Alltoallv(
                block_buffer.data(), recv_counts.data(), recv_offsets.data(), mpi_type<Real>(),
                space_buffer.data(), send_counts.data(), send_offsets.data(), mpi_type<Real>(),
                comm
            )
        );
    }
    if(enable_timing) t_mpi += MPI_Wtime();

    // unpack into local space-distributed field 
    if(enable_timing) t_unpack = -MPI_Wtime();
    unpack_space_distr_field(space_distr, num_tiles);
    if(enable_timing) {
        t_unpack += MPI_Wtime();
        double max_mpi, max_pack, max_unpack;
        MPI_Reduce(&t_mpi, &max_mpi, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(&t_pack, &max_pack, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(&t_unpack, &max_unpack, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        double bytes_sent = static_cast<double>(space_buffer.size() * sizeof(Real)); 
        double bandwidth = 0.0;
        MPI_Reduce(&bytes_sent, &bandwidth, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (rank == 0) {
            bandwidth /= (max_mpi * 1.0e6); // MB/s 
            std::cout << "[ReMap3D::from_block_to_space_distributed] Max: pack = " << max_pack << " s, MPI = " << max_mpi 
                << " s, unpack = " << max_unpack << " s, total = " << (max_pack + max_mpi + max_unpack) << " s. "
                << "MPI aggregate effective BW = " << bandwidth << " MB/s\n" << std::flush;
        }
    }
}


void ReMap3D::init(
    PetscInt block_size_, PetscInt tile_size_,
    PetscInt num_tiles_per_block_)
{ 
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    assert((block_size_ / tile_size_) * tile_size_ == block_size_);

    if (num_tiles_per_block_ == -1) {
        num_tiles_per_block_ = size;
        if (num_tiles_per_block_ * tile_size_ > block_size_) {
            num_tiles_per_block_ = block_size_ / tile_size_;
        }
    }
    assert(num_tiles_per_block_ == block_size_ / tile_size_);
    assert(block_grid->getMPISize() == 1);
    assert(space_grid->getGlobalNumNodes() == block_grid->getLocalNumNodes());

    block_size = block_size_;
    tile_size = tile_size_;
    num_tiles_per_block = num_tiles_per_block_;

    // collect starts and dims from other ranks
    PetscInt* local_starts_ = new PetscInt[size * 3];
    PetscInt* local_sizes_ = new PetscInt[size * 3];
    std::array<PetscInt,3> my_start = space_grid->getGlobalStarts();
    std::array<PetscInt,3> my_dim   = space_grid->getLocalSizes();
    MPI_Allgather(
        my_start.data(), 3, MPIU_INT,
        &local_starts_[0], 3, MPIU_INT, comm
    );
    MPI_Allgather(
        my_dim.data(), 3, MPIU_INT,
        &local_sizes_[0], 3, MPIU_INT, comm
    );

    local_starts.resize(size);
    local_sizes.resize(size);
    for (int r = 0; r < size; r++) {
        local_sizes[r] = std::array<PetscInt,3>{
            local_sizes_[r*3 + 0], 
            local_sizes_[r*3 + 1], 
            local_sizes_[r*3 + 2]
        };
        local_starts[r] = std::array<PetscInt,3>{
            local_starts_[r*3 + 0], 
            local_starts_[r*3 + 1], 
            local_starts_[r*3 + 2]
        };
    }

    // compute rank node counts
    per_rank_nodes.resize(size);
    total_global_nodes = space_grid->getGlobalNumNodes();
    for (int r = 0; r < size; ++r) 
        per_rank_nodes[r] = local_sizes[r][0] * local_sizes[r][1] * local_sizes[r][2];
    local_nodes = per_rank_nodes[rank];

    send_counts.resize(size);
    send_offsets.resize(size + 1);
    recv_counts.resize(size);
    recv_offsets.resize(size + 1);

    // each rank sends (local_nodes * tile_size) values to each rank
    // each rank receives (their_local_nodes * tile_size) from every other rank 
    for (int r = 0; r < size; ++r) {
        send_counts[r] = local_nodes * tile_size; 
        recv_counts[r] = per_rank_nodes[r] * tile_size;
        
    }
    recv_offsets[0] = 0;
    send_offsets[0] = 0;
    for (int r = 0; r < size; ++r) {
        send_offsets[r+1] = send_offsets[r] + send_counts[r]; 
        recv_offsets[r+1] = recv_offsets[r] + recv_counts[r];
    }

    // allocate communication buffers
    // std::cout << "[Rank " << rank << "] init: space_buffer(" 
    //     << num_tiles_per_block * local_nodes * tile_size << ", " << send_offsets.back() + send_counts.back()
    //     << "), block_buffer(" << total_global_nodes * tile_size << ", " << recv_offsets.back() + recv_counts.back() 
    //     << ")\n" << std::flush;
    space_buffer.resize(size * local_nodes * tile_size); //(send_offsets.back() + send_counts.back());   // num_tiles_per_block * local_nodes * tile_size
    block_buffer.resize(total_global_nodes * tile_size); //(recv_offsets.back() + recv_counts.back()); // global_nodes * tile_size

    // uniform check
    uniform = true;
    for (int r = 1; r < size; ++r) {
        if (recv_counts[r] != recv_counts[0]) { 
            uniform = false; 
            break; 
        }
    }
}



void ReMap3D::init(
    Grid_ptr space_grid_, Grid_ptr block_grid_,
    PetscInt block_size_, PetscInt tile_size_,
    PetscInt num_tiles_per_block_)
{
    space_grid = space_grid_;
    block_grid = block_grid_;
    init(block_size_, tile_size_, num_tiles_per_block_);
}


void ReMap3D::pack_space_distr_field(Field_ptr field, PetscInt num_tiles)
{
    // auto grid = field->getGrid();
    const PetscInt block_size = field->getBlockSize();
    const PetscInt block_part = block_size / num_tiles_per_block;
    const PetscInt rank_offset = local_nodes * tile_size; // offset per destination rank in space_buffer
    // const PetscInt xm = space_grid->getLocalSizeX();
    const PetscInt ym = space_grid->getLocalSizeY();
    const PetscInt zm = space_grid->getLocalSizeZ();

    if(num_tiles_per_block < size) {
        assert(num_tiles == 0);
        if(num_tiles > 0) MPI_Abort(comm, -1);
    }

    space_grid->parallel_for(
        [&](PetscInt i, PetscInt j, PetscInt k) {
            Real* block_ptr = field->block(i,j,k);
            PetscInt node_idx = (i * ym + j) * zm + k; 
            PetscInt data_offset = node_idx * tile_size;
            // copy the slice for each partition r into the buffer segment for dest rank r
            for (PetscInt r = 0; r < num_tiles_per_block; ++r) {
                PetscInt block_start = /*tile_size * (r + num_tiles);*/ (block_part * r) + (num_tiles * tile_size); // TOFIX
                PetscInt out_offset = r * rank_offset + data_offset;
                for (PetscInt t = 0; t < tile_size; ++t) {
                    space_buffer[out_offset + t] = block_ptr[block_start + t];
                }
            }
        }
    );
}

void ReMap3D::unpack_space_distr_field(Field_ptr space_field, PetscInt num_tiles)
{
    // auto grid = space_field->getGrid();
    const PetscInt block_size = space_field->getBlockSize();
    const PetscInt block_part = block_size / num_tiles_per_block;
    const PetscInt xm = space_grid->getLocalSizeX();
    const PetscInt ym = space_grid->getLocalSizeY();
    const PetscInt zm = space_grid->getLocalSizeZ();
    const PetscInt rank_offset = local_nodes * tile_size;

    space_grid->parallel_for(
        [&](PetscInt i, PetscInt j, PetscInt k) {
            Real* block_ptr = space_field->block(i,j,k);
            PetscInt local_node = (i * ym + j) * zm + k; 
            PetscInt data_offset = local_node * tile_size;
            // copy the slice for each partition r
            for (PetscInt r = 0; r < num_tiles_per_block; ++r) {
                PetscInt block_start = (block_part * r) + (num_tiles * tile_size);
                PetscInt src_offset = r * rank_offset + data_offset;
                for (PetscInt t = 0; t < tile_size; ++t) {
                    block_ptr[block_start + t] = space_buffer[src_offset + t];
                }
            }
        }
    );
}


void ReMap3D::pack_block_distr_field(Field_ptr field)
{
    // if (size >= num_tiles_per_block) return;

    // #pragma omp parallel for
    for (int r = 0; r < size; ++r) {
        const auto &s = local_starts[r];
        const auto &d = local_sizes[r];
        PetscInt xs = s[0], ys = s[1], zs = s[2];
        PetscInt xm = d[0], ym = d[1], zm = d[2];
        PetscInt rank_offset = recv_offsets[r];

        // iterate over local block in global indices
        // #pragma omp parallel for collapse(3)
        for (PetscInt k = 0; k < zm; ++k) {
            for (PetscInt j = 0; j < ym; ++j) {
                for (PetscInt i = 0; i < xm; ++i) {
                    Real* block_ptr = field->block(xs + i, ys + j, zs + k); 
                    PetscInt node_local = (i * ym + j) * zm + k;
                    PetscInt out_offset = rank_offset + node_local * tile_size;
                    for (PetscInt t = 0; t < tile_size; ++t) {
                        block_buffer[out_offset + t] = block_ptr[t];
                    }
                }
            }
        }
    }
}

void ReMap3D::unpack_block_distr_field(Field_ptr field)
{
    // if (size >= num_tiles_per_block) return;

    // For each rank r, take the slice at recv_offsets[r] and write into serial at global starts[r]
    // #pragma omp parallel for
    for (PetscInt r = 0; r < size; ++r) {
        const auto &s = local_starts[r];
        const auto &d = local_sizes[r];
        PetscInt xs = s[0], ys = s[1], zs = s[2];
        PetscInt xm = d[0], ym = d[1], zm = d[2];
        PetscInt rank_offset = recv_offsets[r];

        // #pragma omp parallel for collapse(3)
        for (PetscInt k = 0; k < zm; ++k) {
            for (PetscInt j = 0; j < ym; ++j) {
                for (PetscInt i = 0; i < xm; ++i) {
                    Real* block_ptr = field->block(xs + i, ys + j, zs + k);
                    PetscInt node_local = (i * ym + j) * zm + k;
                    PetscInt src = rank_offset + node_local * tile_size;
                    for (PetscInt t = 0; t < tile_size; ++t) {
                        block_ptr[t] = block_buffer[src + t];
                    }
                }
            }
        }
    }
}