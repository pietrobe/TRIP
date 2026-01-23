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
    
    std::cout << "[Rank " << rank << "] Field " << name << " has " << n_local_blocks 
                << " blocks of size " << block_size << " (" << start_index << "-" << end_index << ") with local grid size " 
                << "(" << xm << "," << ym << "," << zm << ")\n" << std::flush; 
}

Grid_ptr Field::getGrid() const { return grid; }

const PetscInt Field::getBlockSize() const { return block_size; }

const PetscInt Field::getNumLocalBlocks() const { return n_local_blocks; }

Vec Field::getVec() { 
    if (!is_Vec_allocated) throw std::runtime_error("getVec(): PETSc Vec not allocated for field " + name);
    return vec; 
}

Real* Field::getData() { return data_host; }

const PetscInt Field::getStartIndex() const { return start_index; }
const PetscInt Field::getEndIndex() const { return end_index; }

std::string Field::getName() const { return name; }
PetscInt Field::getNPhysical() const { return Nphysical; } 
PetscInt Field::getNparam() const { return Nparam; }
std::vector<PetscInt> Field::getNParamSizes() const { return Nparam_size; }
PetscInt Field::getParamSize(int p) const { 
    if (Nparam_size.size() == 0 || Nparam_size.size() < p)
        throw std::runtime_error("In Field::getParamSize, p is out of Nparam size.");
    return Nparam_size[p];
}


Real* Field::block(PetscInt i,PetscInt j,PetscInt k)
{
    assert(i >= 0 && i < xm);
    assert(j >= 0 && j < ym);
    assert(k >= 0 && k < zm);
    assert(data_host != nullptr);

    PetscInt linear_block = (i * ym + j) * zm + k;
    PetscInt idx = linear_block * block_size;
    return &data_host[idx];
}
Real* Field::operator()(PetscInt i,PetscInt j,PetscInt k) { return block(i,j,k); }

Real& Field::ref(PetscInt i,PetscInt j,PetscInt k)
{
    return block(i,j,k)[0];
}

std::vector<PetscInt> Field::block_to_local(PetscInt b) const 
{
    if (b < 0 || b >= n_local_blocks * block_size) throw std::out_of_range("Block index out of bounds in Field::block_to_local.");
    std::vector<PetscInt> indices(Nparam + 1); // +1 for physical component
    PetscInt block_index = b / Nphysical;  
    for (int i = Nparam - 1; i > 0; --i) {
        indices[i] = block_index % Nparam_size[i];
        block_index /= Nparam_size[i];
    }
    indices[0] = block_index;
    indices.back() = b % Nphysical; 
    return indices;
}


PetscInt Field::local_to_block(const std::vector<PetscInt>& idx) const 
{
    if (idx.size() != ((Nparam == 0) ? 1 : Nparam)) {
        throw std::runtime_error("Number of indices does not match layout of Field.");
    }
    if (Nparam == 0 && (idx[0] < 0 || idx[0] >= Nphysical)) throw std::out_of_range("Index out of bounds in local_to_block.");
    for (int k = 0; k < Nparam; ++k) {
        if (idx[k] < 0 || idx[k] >= Nparam_size[k]) {
            throw std::out_of_range("Index out of bounds in local_to_block.");
        }
    }

    PetscInt block_index = 0;
    int stride = 1;
    if(Nparam == 0) return idx[0];

    for (int i = idx.size() - 1; i >= 0; --i) {
        block_index += idx[i] * stride;
        stride *= Nparam_size[i];
    }
    return Nphysical * block_index;
}

void Field::allocate(bool allocate_vec)
{
    MPI_Comm comm = grid->getComm();
    PetscErrorCode ierr;
    PetscInt local_dofs = xm * ym * zm * block_size;
    n_local_blocks = local_dofs / block_size;
    if (allocate_vec) {
        ierr = VecCreate(comm, &vec);  CHKERRABORT(comm, ierr); // TODO check memory, flag per abiliate/disabilitare
        ierr = VecSetSizes(vec, local_dofs, PETSC_DECIDE);  CHKERRABORT(comm, ierr); // TOFIX check if it allocates memory
        ierr = VecSetFromOptions(vec);  CHKERRABORT(comm, ierr);
        // Bind host array to Vec, so that data_host is the Vec’s internal storage
        ierr = PetscCalloc(local_dofs * sizeof(Real), &data_host); CHKERRABORT(comm, ierr);
        ierr = VecReplaceArray(vec, data_host);  CHKERRABORT(comm, ierr);
        ierr = VecGetOwnershipRange(vec, &start_index, &end_index); CHKERRABORT(comm, ierr);
    } else {
        data_host = new Real[local_dofs];
        std::fill(data_host, data_host + local_dofs, 0.0);
        PetscInt offset = 0;
        MPI_Exscan(&local_dofs, &offset, 1, MPIU_INT, MPI_SUM, comm); // check
        start_index = offset; 
        end_index = start_index + local_dofs; 
    }
}