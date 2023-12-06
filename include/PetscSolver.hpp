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


#ifdef WITH_PETSCSOLVER

#ifndef __PETSCSOLVER_HEADER__
#define __PETSCSOLVER_HEADER__

#include "LinearSolver.hpp"
#include "PETScBSolverPS.h"
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
    void Allocate(const OCPMatrix& mat) override;

    /// Calculate terms used in communication
    void CalCommTerm(const USI& actWellNum, const Domain* domain) override;

    /// Assemble coefficient matrix.
    void AssembleMat(OCPMatrix& mat) override;

    /// Get number of iterations used by iterative solver.
    USI GetNumIters() const override { return 1; }

protected:

    // CSR/BSR   
    OCP_INT                 nb;         ///< block dim              
    vector<OCP_DBL>         A;          ///< value
    vector<OCP_SLL>         iA;         ///< row ptr
    vector<OCP_SLL>         jA;         ///< col index  
    OCP_DBL*                b;          ///< rhs
    OCP_DBL*                x;          ///< solution

    // Communication
    MPI_Comm                myComm;       ///< Communicator
    OCP_INT                 numproc;      ///< num of process
    OCP_INT                 myrank;       ///< current rank
    const vector<OCP_ULL>*  global_index; ///< global index
    vector<OCP_SLL>         allBegin;     ///< begin for all process(self-include) 
    vector<OCP_SLL>         allEnd;       ///< end for all process(self-include)
    vector<OCP_INT>         allEle;       ///< num of elements for each proc
};


class ScalarPetscSolver : public PetscSolver
{
public:
    ScalarPetscSolver() {};
    /// Solve the linear system.
    OCP_INT Solve() override;

};

class VectorPetscSolver : public PetscSolver
{
public:
    VectorPetscSolver(const Domain* domain) { 
        myComm  = domain->myComm;
        numproc = domain->numproc;
        myrank  = domain->myrank;
        allBegin.resize(numproc);
        allEnd.resize(numproc);
        allEle.resize(numproc);
    }
    /// Solve the linear system.
    OCP_INT Solve() override;

};

#endif 

#endif // WITH_PETSCSOLVER


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           May/04/2023      Create file                          */
/*----------------------------------------------------------------------------*/