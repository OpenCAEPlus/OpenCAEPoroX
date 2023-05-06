/*! \file    PetscSolver.hpp
 *  \brief   API of Petsc Solver
 *  \author  Shizhe Li
 *  \date    May/04/2023
 *
 *  \note    use Petsc as linear solver
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


#ifdef WITH_PETSC

#ifndef __PETSCSOLVER_HEADER__
#define __PETSCSOLVER_HEADER__

#include "LinearSolver.hpp"

#include <vector>

using namespace std;

/// API for Petsc Solver
// Note: the index is 0-baesd
class PetscSolver : public LinearSolver
{
public:
    PetscSolver() = default;
    /// Set parameters.
    void SetupParam(const string& dir, const string& file) override {};

    /// Initialize the Params for linear solver.
    void InitParam() override {};

    /// Allocate memoery for pardiso solver
    void Allocate(const OCP_USI& max_nnz, const OCP_USI& maxDim) override;

    /// Calculate terms used in communication
    void CalCommTerm(const USI& actWellNum, const Domain* domain) override;

    /// Assemble coefficient matrix.
    void AssembleMat(const vector<vector<USI>>& colId,
                     const vector<vector<OCP_DBL>>& val,
                     const OCP_USI& dim,
                     vector<OCP_DBL>& rhs,
                     vector<OCP_DBL>& u) override; 

    /// Get number of iterations used by iterative solver.
    USI GetNumIters() const override { return 1; }

protected:

    // CSR/BSR   
    int                     blockdim;   ///< block dim
    int                     dim_local;  ///< local dim               
    vector<double>          A;          ///< value
    vector<int>             iA;         ///< row ptr
    vector<int>             jA;         ///< col index  
    double*                 b;          ///< rhs
    double*                 x;          ///< solution

    // Communication
    int                     dim_global;  ///< total dim
    int                     row_begin;  ///< global row begin
    const vector<OCP_USI>*  global_index; ///< global index
};


class ScalarPetscSolver : public PetscSolver
{
public:
    ScalarPetscSolver() { blockdim = 1; }
    /// Solve the linear system.
    OCP_INT Solve() override {};

};

class VectorPetscSolver : public PetscSolver
{
public:
    VectorPetscSolver(const USI& blockDim) { blockdim = blockDim; }
    /// Solve the linear system.
    OCP_INT Solve() override {};

};

#endif 

#endif // WITH_PETSC


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           May/04/2023      Create file                          */
/*----------------------------------------------------------------------------*/