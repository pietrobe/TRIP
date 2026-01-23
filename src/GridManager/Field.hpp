#ifndef FIELD_HH
#define FIELD_HH 

#include "Grid3D.hpp"
#include <algorithm>
#include <memory>
#include <vector>
#include <array>
#include <stdexcept>
#include <iostream>
#include <petscerror.h>

using Real = double;
using Grid_ptr  = std::shared_ptr<Grid3D>;

/**
 * @brief Field class represents a physical quantity defined on a Grid3D.
 * Each spatial cell (i,j,k) on the grid contains a block of block_size values.
 * Each MPI process has a Field object, the PETSc Vecs are distributed so that
 * each rank only stores its local chunck of the full field.
 * @note that Field uses a z-fastest mapping.
 */
class Field {

public:
    /** @brief constructor for a field with a specific layout. The constructor already allocates the memory required by the field.
     * For example, for Stokes parameters in the polarised case, Nphysical = 4 and Nparam_size={Ntheta, Nchi, Nnu};
     * for Stokes parameters in the unpolarised case, Nphysical = 1 and Nparam_size={Ntheta, Nchi, Nnu}.
     * @note the order of the Nparam_size matters!
     * @param fieldName name of the field.
     * @param grid_ pointer to associated local Grid3D.
     * @param Nphysical_ number of physical quantities per spatial cell.
     * @param Nparam_size vector of sizes of each independent parameters. 
     * @param allocate_PETSc_vec flag indicating whether PETSc Vec is needed (set to true) or not (set to false). 
     */
    Field(
        const std::string &fieldName, 
        Grid_ptr grid_,
        const PetscInt Nphysical_, 
        const std::vector<PetscInt>& Nparam_size_,
        const bool allocate_PETSc_vec)
        : name(fieldName), 
          grid(grid_), 
          Nphysical(Nphysical_), 
          Nparam(Nparam_size_.size()), 
          Nparam_size(Nparam_size_),
          xm(grid_->getLocalSizeX()), 
          ym(grid_->getLocalSizeY()), 
          zm(grid_->getLocalSizeZ()),
          is_Vec_allocated(allocate_PETSc_vec)
    {
        block_size = Nphysical_;
        for (auto s: Nparam_size_) block_size *= s;
        allocate(allocate_PETSc_vec);
    }

    /** @brief construct for a 'scalar' field. 
     * For 'scalar' we mean a field with Nphysical = 1 and Nparam = 0 
     * (for example the atmospheric quantities that depends only on the spatial cell).
     * In such case, block_size = 1.
     * @param fieldName name of the field.
     * @param grid_ pointer to associated local Grid3D.
     * @param allocate_PETSc_vec flag indicating whether PETSc Vec is needed (set to true) or not (set to false). 
     */
    Field(
        const std::string &fieldName, 
        Grid_ptr grid_, 
        const bool allocate_PETSc_vec = false)
        : name(fieldName), grid(grid_), 
          xm(grid_->getLocalSizeX()), 
          ym(grid_->getLocalSizeY()), 
          zm(grid_->getLocalSizeZ()),
          is_Vec_allocated(allocate_PETSc_vec)
    {
        allocate(allocate_PETSc_vec);
    }

    /** @brief construct for a 'scalar multi-physical' field. 
     * For 'scalar multi-physical' we mean a field with and arbitrary number of Nphysical and Nparam = 0
     * (for example the magnetic field which is defined in polar coordinates and the bulk velocity).
     * In such case, block_size = Nphysical.
     * @param fieldName name of the field.
     * @param grid_ pointer to associated local Grid3D.
     * @param Nphysical_ number of physical quantities per spatial cell.
     * @param allocate_PETSc_vec flag indicating whether PETSc Vec is needed (set to true) or not (set to false). 
     */
    Field(
        const std::string &fieldName, 
        Grid_ptr grid_, 
        const PetscInt Nphysical_,
        const bool allocate_PETSc_vec)
        : name(fieldName), grid(grid_), 
          Nphysical(Nphysical_), 
          block_size(Nphysical_),
          xm(grid_->getLocalSizeX()), 
          ym(grid_->getLocalSizeY()), 
          zm(grid_->getLocalSizeZ()),
          is_Vec_allocated(allocate_PETSc_vec)
    {
        allocate(allocate_PETSc_vec);
    }

    /** @brief destructor that release PETSc Vecs. */
    ~Field() 
    {
        if (is_Vec_allocated) {
            // unbind host array from PETSc Vec
            // PetscErrorCode ierr = VecReplaceArray(vec, nullptr); CHKERRABORT(grid->getComm(), ierr);  
            // destroy PETSc Vec
            PetscErrorCode ierr = VecRestoreArray(vec, &data_host); 
            ierr = VecDestroy(&vec); CHKERRABORT(grid->getComm(), ierr);
            // ierr = PetscFree(data_host);CHKERRABORT(grid->getComm(), ierr);
            vec = nullptr;
        } else {
            delete[] data_host; 
        }
        data_host = nullptr;
    }

    void print_info();
    
    /** @brief return local grid. */ 
    Grid_ptr getGrid() const;

    /** @brief return number of values per block. */ 
    const PetscInt getBlockSize() const;

    /** @brief return number of blocks in local vector. */ 
    const PetscInt getNumLocalBlocks() const;

    /** @brief return local PETSc vector. 
     * @note any changes done to PETSc vector are automatically synched to the 'standard' vector.
     */ 
    Vec getVec();

    /** @brief return pointer to local data. 
     * @note any changes done to this vector are automatically synched to the PETSc Vec (if allocated).
     */ 
    Real* getData();

    /** @brief return start and end indices of the local vector. */ 
    const PetscInt getStartIndex() const;
    const PetscInt getEndIndex() const;

    std::string getName() const;
    PetscInt getNPhysical() const;
    PetscInt getNparam() const;
    std::vector<PetscInt> getNParamSizes() const;
    PetscInt getParamSize(int p) const;
    
    /** @brief returns pointer to block (i,j,k) */
    Real* block(PetscInt i,PetscInt j,PetscInt k);

    /** @brief returns reference to first element in block (i,j,k) */
    Real& ref(PetscInt i,PetscInt j,PetscInt k);

    /** @brief same as block(i,j,k) but with less verbose syntax (field(i,j,k) instead of field.block(i,j,k)) */
    Real* operator()(PetscInt i,PetscInt j,PetscInt k);

    /** @brief convert linear block index to multi-dimensional indices */
    std::vector<PetscInt> block_to_local(PetscInt b) const;

    /** @brief convert multi-dimensional indices to linear block index */
    PetscInt local_to_block(const std::vector<PetscInt>& idx) const;

    template<typename... Indices>
    PetscInt local_to_block(Indices... indices) const 
    {
        std::array<int, sizeof...(Indices)> idx{indices...};
        if (idx.size() != ((Nparam == 0) ? 1 : Nparam)) {
            throw std::runtime_error("Number of indices does not match layout of Field.");
        }
        // out of bound indices
        if (Nparam == 0 && (idx[0] < 0 || idx[0] >= Nphysical)) throw std::out_of_range("Index out of bounds in local_to_block.");
        for (int k = 0; k < Nparam; ++k) {
            if (idx[k] < 0 || idx[k] >= Nparam_size[k]) {
                throw std::out_of_range("Index out of bounds in local_to_block.");
            }
        }

        int block_index = 0;
        int stride = 1;
        if(Nparam == 0) return idx[0];

        for (int i = idx.size() - 1; i >= 0; --i) {
            block_index += idx[i] * stride;
            stride *= Nparam_size[i];
        }
        return Nphysical * block_index;
    }

private:
    std::string name;
    Grid_ptr grid; // pointer to associated Grid3D
    Vec vec = nullptr; // PETSc distributed global vector
    bool is_Vec_allocated = false; 
    Real* data_host = nullptr; // host array bound to Vec
    PetscInt xm = 0, ym = 0, zm = 0; // local grid sizes
    PetscInt start_index = 0; // start index of local PETSc vector
    PetscInt end_index = 0; // end index of local PETSc vector

    // generic block layout
    PetscInt block_size = 1; 
    PetscInt n_local_blocks = 0; // number of blocks in the local PETSc vector
    PetscInt Nphysical = 1; // number of physical quantities per spatial cell 
    PetscInt Nparam = 0; // number of independent values for each physical quantities 
    std::vector<PetscInt> Nparam_size = {}; // size of each independent values. The order matters!


    /** @brief allocate PETSc Vec. */
    void allocate(bool allocate_vec);
    
};

#endif