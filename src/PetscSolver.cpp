/*! \file    PetscSolver.cpp
 *  \brief   API of Petsc Solver
 *  \author  Shizhe Li
 *  \date    May/04/2023
 *
 *  \note    use SAMG as linear solver
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifdef WITH_PETSCSOLVER

#include "PetscSolver.hpp"


 /// Allocate memoery for pardiso solver
void PetscSolver::Allocate(const OCP_USI& max_nnz, const OCP_USI& maxDim)
{
    A.resize(max_nnz * blockdim * blockdim);
    iA.resize(maxDim + 1);
    jA.resize(max_nnz);
}


/// Calculate terms used in communication
void PetscSolver::CalCommTerm(const USI& actWellNum, const Domain* domain)
{
    global_index  = domain->CalGlobalIndex(actWellNum);
    const OCP_INT numElementloc = actWellNum + domain->GetNumGridInterior();

    GetWallTime timer;
    timer.Start();

    MPI_Allgather(&numElementloc, 1, MPI_INT, &allEle[0], 1, MPI_INT, myComm);

    OCPTIME_COMM_COLLECTIVE += timer.Stop() / 1000;   
    
    allBegin[0] = 0;
    allEnd[0]   = allEle[0] - 1;
    for (OCP_USI p = 1; p < numproc; p++) {
        allBegin[p] = allEnd[p - 1] + 1;
        allEnd[p]   = allBegin[p] + allEle[p] - 1;
    }
}


/// Assemble coefficient matrix.
void PetscSolver::AssembleMat(const vector<vector<USI>>& colId,
    const vector<vector<OCP_DBL>>& val,
    const OCP_USI& dim,
    vector<OCP_DBL>& rhs,
    vector<OCP_DBL>& u)
{

    const USI blockSize = blockdim * blockdim;
    vector<OCP_INT> tmpJ;
    // Assemble iA, jA, A
    iA[0] = 0;
    for (OCP_USI i = 1; i < dim + 1; i++) {
        const USI nnzR = colId[i - 1].size();

        tmpJ.resize(nnzR);
        for (USI j = 0; j < nnzR; j++) {
            tmpJ[j] = global_index->at(colId[i - 1][j]);
        }

        iA[i] = iA[i - 1] + nnzR;
        copy(tmpJ.begin(), tmpJ.end(), &jA[iA[i - 1]]);
        copy(val[i - 1].begin(), val[i - 1].end(), &A[(iA[i - 1]) * blockSize]);
    }


    b = rhs.data();
    x = u.data();
}


OCP_INT ScalarPetscSolver::Solve()
{
    return IMPEC_solver_p(myrank, numproc, allBegin.data(), allEnd.data(), iA.data(), jA.data(), A.data(), b, x);
}


OCP_INT VectorPetscSolver::Solve()
{
    return FIM_solver_p(myrank, numproc, blockdim, allBegin.data(), allEnd.data(), iA.data(), jA.data(), A.data(), b, x);
}

#endif // WITH_PETSCSOLVER



 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           May/04/2023      Create file                          */
 /*----------------------------------------------------------------------------*/