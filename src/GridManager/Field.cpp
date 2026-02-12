#include "Grid3D.hpp"
#include "Field.hpp"
#include <algorithm>
#include <vector>
#include <array>
#include <stdexcept>
#include <iostream>
#include <petscerror.h>
#include <cassert>


void Field::print_info()
{
    int rank;
    MPI_Comm_rank(grid->getComm(), &rank);

    std::cout << "[Rank " << rank << "] Field " << name
                << "\n  owned dofs:   (" << start_index << ", " << end_index << ")"
                << "\n  ghosted dofs: (" << ghosted_start_index << ", " << ghosted_end_index << ")"
                << "\n  local ghosted grid size = ("
                << grid->getLocalSizeWithGhostX() << ","
                << grid->getLocalSizeWithGhostY() << ","
                << grid->getLocalSizeWithGhostZ() << ")\n"
                << std::flush;
}

Grid_ptr Field::getGrid() const { return grid; }

Int Field::getBlockSize() const { return block_size; }

Int Field::getNumLocalBlocks() const { return n_local_blocks; }

Vec Field::getVec() { 
    if (!is_Vec_allocated) throw std::runtime_error("getVec(): PETSc Vec not allocated for field " + name);
    return vec; 
}

Real* Field::getData() { return data_host; }

Int Field::getStartIndex() const { return start_index; }
Int Field::getEndIndex() const { return end_index; }

std::string Field::getName() const { return name; }
Int Field::getNPhysical() const { return Nphysical; } 
Int Field::getNparam() const { return Nparam; }
std::vector<Int> Field::getNParamSizes() const { return Nparam_size; }
Int Field::getParamSize(int p) const { 
    if (Nparam_size.size() == 0 || Nparam_size.size() < p)
        throw std::runtime_error("In Field::getParamSize, p is out of Nparam size.");
    return Nparam_size[p];
}


Real* Field::block(Int i,Int j,Int k)
{
    auto [xm, ym, zm] = grid->getLocalSizesWithGhost();
    if (i < 0 || i >= xm) throw std::out_of_range("i index out of bounds in Field::block() for " + name);
    if (j < 0 || j >= ym) throw std::out_of_range("j index out of bounds in Field::block() for " + name);
    if (k < 0 || k >= zm) throw std::out_of_range("k index out of bounds in Field::block() for " + name);
    if (!data_host) throw std::runtime_error("Data not allocated in Field::block() for " + name);

    Int linear_block = (i * ym + j) * zm + k; // z-fastest
    Int idx = linear_block * block_size;
    return &data_host[idx];
}
Real* Field::operator()(Int i,Int j,Int k) { return block(i,j,k); }

Real& Field::ref(Int i,Int j,Int k)
{
    return block(i,j,k)[0];
}

std::vector<Int> Field::block_to_local(Int b) const 
{
    if (b < 0 || b >= n_local_blocks * block_size) throw std::out_of_range("Block index out of bounds in Field::block_to_local.");
    std::vector<Int> indices(Nparam + 1); // +1 for physical component
    Int block_index = b / Nphysical;  
    for (Int i = Nparam - 1; i > 0; --i) {
        indices[i] = block_index % Nparam_size[i];
        block_index /= Nparam_size[i];
    }
    indices[0] = block_index;
    indices.back() = b % Nphysical; 
    return indices;
}


Int Field::local_to_block(const std::vector<Int>& idx) const 
{
    if (idx.size() != ((Nparam == 0) ? 1 : Nparam)) {
        throw std::runtime_error("Number of indices does not match layout of Field.");
    }
    if (Nparam == 0 && (idx[0] < 0 || idx[0] >= Nphysical)) throw std::out_of_range("Index out of bounds in local_to_block.");
    for (Int k = 0; k < Nparam; ++k) {
        if (idx[k] < 0 || idx[k] >= Nparam_size[k]) {
            throw std::out_of_range("Index out of bounds in local_to_block.");
        }
    }

    Int block_index = 0;
    Int stride = 1;
    if(Nparam == 0) return idx[0];

    for (Int i = idx.size() - 1; i >= 0; --i) {
        block_index += idx[i] * stride;
        stride *= Nparam_size[i];
    }
    return Nphysical * block_index;
}

void Field::set_to_zero() 
{
    Int n = grid->getLocalNumNodesWithGhosts() * block_size;
    std::fill(data_host, data_host + n, Real(0));
}


void Field::allocate(bool allocate_vec) // TODO : fix in cause ghost layers will actually be used
{
    MPI_Comm comm = grid->getComm();
    PetscErrorCode ierr;
    Int ghosted_local_dofs = grid->getLocalNumNodesWithGhosts() * block_size;
    Int local_dofs = grid->getLocalNumNodes() * block_size;
    n_local_blocks = local_dofs / block_size;
    ghosted_start_index = 0;
    ghosted_end_index   = ghosted_local_dofs;
    if (allocate_vec) {
        ierr = VecCreate(comm, &vec);  CHKERRABORT(comm, ierr); 
        ierr = VecSetSizes(vec, local_dofs, PETSC_DECIDE);  CHKERRABORT(comm, ierr); 
        ierr = VecSetFromOptions(vec);  CHKERRABORT(comm, ierr);
        // Bind host array to Vec, so that data_host is the Vec’s internal storage
        ierr = PetscCalloc(ghosted_local_dofs * sizeof(Real), &data_host); CHKERRABORT(comm, ierr);
        ierr = VecReplaceArray(vec, data_host);  CHKERRABORT(comm, ierr);
        ierr = VecGetOwnershipRange(vec, &start_index, &end_index); CHKERRABORT(comm, ierr);
    } else {
        data_host = new Real[ghosted_local_dofs];
        std::fill(data_host, data_host + ghosted_local_dofs, 0.0);
        Int offset = 0;
        MPI_Exscan(&ghosted_local_dofs, &offset, 1, MPIU_INT, MPI_SUM, comm); 
        start_index = offset; 
        end_index = start_index + local_dofs; 
    }
}