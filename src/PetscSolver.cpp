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

#ifdef WITH_PETSC

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
    global_index = domain->CalGlobalIndex(actWellNum);

    const OCP_INT numElementLoc = actWellNum + domain->GetNumGridInterior();
    const OCP_INT global_end = global_index->at(numElementLoc - 1);

    // Get global row start
    row_begin = global_end - numElementLoc + 1;

    // Get global Dimension
    dim_global = global_end + 1;
    MPI_Bcast(&dim_global, 1, MPI_INT, domain->numproc - 1, domain->myComm);
}


/// Assemble coefficient matrix.
void PetscSolver::AssembleMat(const vector<vector<USI>>& colId,
    const vector<vector<OCP_DBL>>& val,
    const OCP_USI& dim,
    vector<OCP_DBL>& rhs,
    vector<OCP_DBL>& u)
{
    dim_local           = dim;
    const USI blockSize = blockdim * blockdim;
    vector<OCP_USI> tmpJ;
    // Assemble iA, jA, A
    iA[0] = 0;
    for (OCP_USI i = 1; i < dim + 1; i++) {
        const USI nnzR = colId[i - 1].size();

        tmpJ.resize(nnzR);
        for (USI j = 0; j < nnzR; j++) {
            tmpJ[j] = global_index->at(colId[i - 1][j]);
        }

        iA[i] = iA[i - 1] + nnzR;
        copy(&jA[iA[i - 1]], &jA[iA[i - 1]] + nnzR, tmpJ.data());
        const OCP_DBL* begin = &val[i - 1][0];
        const OCP_DBL* end = begin + nnzR * blockSize;
        copy(begin, end, &A[(iA[i - 1]) * blockSize]);
    }


    b = rhs.data();
    x = u.data();
}

#endif // WITH_PETSC



 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           May/04/2023      Create file                          */
 /*----------------------------------------------------------------------------*/