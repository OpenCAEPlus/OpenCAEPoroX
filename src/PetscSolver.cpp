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
void PetscSolver::Allocate(const OCPMatrix& mat, const Domain* domain)
{
    nb = mat.nb;
    A.resize(mat.max_nnz * nb * nb);
    iA.resize(mat.maxDim + 1);
    jA.resize(mat.max_nnz);

    // pre-allocation
    allBegin.resize(domain->global_numproc);
    allEnd.resize(domain->global_numproc);
    allEle.resize(domain->global_numproc);
}


/// Calculate terms used in communication
void PetscSolver::CalCommTerm(const Domain* domain)
{
    global_index  = domain->CalGlobalIndex();
    myComm        = domain->cs_comm;

    const OCP_INT numElementloc = domain->GetNumActElementForSolver(); 

    GetWallTime timer;
    timer.Start();

    MPI_Allgather(&numElementloc, 1, OCPMPI_INT, &allEle[0], 1, OCPMPI_INT, myComm);

    OCPTIME_COMM_COLLECTIVE += timer.Stop();
    
    allBegin[0] = 0;
    allEnd[0]   = allEle[0] - 1;
    for (OCP_USI p = 1; p < domain->cs_numproc; p++) {
        allBegin[p] = allEnd[p - 1] + 1;
        allEnd[p]   = allBegin[p] + allEle[p] - 1;
    }
}


/// Assemble coefficient matrix.
void PetscSolver::AssembleMat(OCPMatrix& mat, const Domain* domain)
{
    CalCommTerm(domain);
    nb = mat.nb;

    const USI blockSize = nb * nb;
    vector<OCP_SLL> tmpJ;
    // Assemble iA, jA, A
    iA[0] = 0;
    for (OCP_USI i = 1; i < mat.dim + 1; i++) {
        const USI nnzR = mat.colId[i - 1].size();

        tmpJ.resize(nnzR);
        for (USI j = 0; j < nnzR; j++) {
            tmpJ[j] = global_index->at(mat.colId[i - 1][j]);
        }

        iA[i] = iA[i - 1] + nnzR;
        copy(tmpJ.begin(), tmpJ.end(), &jA[iA[i - 1]]);
        copy(mat.val[i - 1].begin(), mat.val[i - 1].end(), &A[(iA[i - 1]) * blockSize]);
    }

    // fill(mat.u.begin(), mat.u.end(), 0.0);
    b = mat.b.data();
    x = mat.u.data();
}


OCP_INT ScalarPetscSolver::Solve()
{
    OCP_ABORT("Inavailable!");
}


VectorPetscSolver::VectorPetscSolver(const string& dir, const string& file, const OCPMatrix& mat, const Domain* domain)
{
    Allocate(mat, domain);
}


OCP_INT VectorPetscSolver::Solve()
{
    return FIM_solver_p(myComm, nb, allBegin.data(), allEnd.data(), iA.data(), jA.data(), A.data(), b, x);
}

#endif // WITH_PETSCSOLVER



 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           May/04/2023      Create file                          */
 /*----------------------------------------------------------------------------*/