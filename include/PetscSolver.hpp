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
#include <vector>

#if WITH_AS
#include "as.h"
#else
#include "PETScBSolverPS.h"
#endif


using namespace std;

/// API for Petsc Solver
// Note: the index is 0-baesd
class PetscSolver : public LinearSolver
{
public:
    PetscSolver() = default;

    /// Assemble coefficient matrix.
    void AssembleMat(OCPMatrix& mat, const Domain* domain) override;

    /// Get number of iterations used by iterative solver.
    USI GetNumIters() const override { return 1; }

protected:

    /// Allocate memoery for pardiso solver
    void Allocate(const OCPMatrix& mat, const Domain* domain);
    /// Calculate terms used in communication
    void CalCommTerm(const Domain* domain);

protected:

    // CSR/BSR   
    OCP_INT                 nb;         ///< block dim              
    vector<OCP_DBL>         A;          ///< value
    vector<OCP_SLL>         iA;         ///< row ptr
    vector<OCP_SLL>         jA;         ///< col index  
    OCP_DBL*                b;          ///< rhs
    OCP_DBL*                x;          ///< solution

    // Communication
    MPI_Comm                myComm{ MPI_COMM_NULL };  ///< Communicator
    const vector<OCP_ULL>*  global_index;             ///< global index
    vector<OCP_SLL>         allBegin;                 ///< begin for all process(self-include) 
    vector<OCP_SLL>         allEnd;                   ///< end for all process(self-include)
    vector<OCP_INT>         allEle;                   ///< num of elements for each proc
};


class ScalarPetscSolver : public PetscSolver
{
public:
    ScalarPetscSolver(const string& dir, const string& file, const OCPMatrix& mat, const Domain* domain) { OCP_ABORT("Inavailable!"); };
    /// Solve the linear system.
    OCP_INT Solve() override;

};

class VectorPetscSolver : public PetscSolver
{
public:
    VectorPetscSolver(const string& dir, const string& file, const OCPMatrix& mat, const Domain* domain);
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