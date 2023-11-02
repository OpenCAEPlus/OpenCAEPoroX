/*! \file    PardisoSolver.cpp
 *  \brief   PardisoSolver for OpenCAEPoroX simulator
 *  \author  Shizhe Li
 *  \date    Feb/21/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifdef WITH_PARDISO

#include "PardisoSolver.hpp"

#if    OCPFLOATTYPEWIDTH == 64


void PardisoSolver::SetupParam(const string& dir, const string& file)
{
    InitParam();
}


void PardisoSolver::InitParam()
{
    iparm[0]  = 1;  /* Solver default parameters overriden with provided by iparm */
    iparm[1]  = 2;  /* Use METIS for fill-in reordering */
    iparm[5]  = 0;  /* Write solution into x */
    iparm[7]  = 2;  /* Max number of iterative refinement steps */
    iparm[9]  = 13; /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
    iparm[12] = 1;  /* Switch on Maximum Weighted Matching algorithm (default for non-symmetric) */
    iparm[17] = 0;  /* (closed) Output: Number of nonzeros in the factor LU */ 
    iparm[18] = 0;  /* (closed) Output: Mflops for LU factorization */
    iparm[26] = 1;  /* Check input data for correctness */
    iparm[39] = 2;  /* Input: matrix/rhs/solution are distributed between MPI processes  */

    /*Zero-based indexing: columns and rows indexing in arrays ia, ja, and perm starts from 0 (C-style indexing).*/
    iparm[34] = 1;  
}


void PardisoSolver::Allocate(const OCP_USI& max_nnz, const OCP_USI& maxDim)
{
    if (blockdim > 1) {
        iparm[37] = blockdim * blockdim;
    }

    iA.resize(maxDim + 1);
    jA.resize(max_nnz);
    A.resize(max_nnz * blockdim * blockdim);
}


/// Calculate terms used in communication
void PardisoSolver::CalCommTerm(const USI& actWellNum, const Domain* domain)
{

    global_index = domain->CalGlobalIndex(actWellNum);

    const OCP_INT numGridInterior = domain->GetNumGridInterior();
    const OCP_INT numElementLoc   = actWellNum + numGridInterior;
    const OCP_INT global_end      = global_index->at(numElementLoc - 1);

    iparm[40] = global_end - numElementLoc + 1;  // global begin (included)
    iparm[41] = global_end;                      // global end   (included)

    // Get Dimension
    N = global_end + 1;

    GetWallTime timer;
    timer.Start();

    MPI_Bcast(&N, 1, MPI_INT, domain->numproc - 1, domain->myComm);

    OCPTIME_COMM_COLLECTIVE += timer.Stop() / TIME_S2MS;
}


/// Assemble coefficient matrix.
void PardisoSolver::AssembleMat(const vector<vector<USI>>& colId,
    const vector<vector<OCP_DBL>>& val,
    const OCP_USI& dim,
    vector<OCP_DBL>& rhs,
    vector<OCP_DBL>& u)
{
    const USI blockSize = blockdim * blockdim;
    vector<USI> tmp;  
    // Assemble iA, jA, A
    iA[0] = 0;
    for (OCP_USI i = 1; i < dim + 1; i++) {
        const USI nnzR = colId[i - 1].size();

        iA[i] = iA[i - 1] + nnzR;
        // reorder
        tmp = colId[i - 1];
        for (auto& t : tmp)    t = global_index->at(t);
        sort(tmp.begin(), tmp.end());
        for (USI j0 = 0; j0 < nnzR; j0++) {
            jA[iA[i - 1] + j0] = tmp[j0];
            for (USI j1 = 0; j1 < nnzR; j1++) {
                if (global_index->at(colId[i - 1][j1]) == tmp[j0]) {
                    const OCP_DBL* begin = &val[i - 1][0] + j1 * blockSize;
                    const OCP_DBL* end = begin + blockSize;
                    copy(begin, end, &A[(iA[i - 1] + j0) * blockSize]);
                    break;
                }
            }
        }           
    }


    b = rhs.data();
    x = u.data();
}


OCP_INT PardisoSolver::Solve()
{

    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates   */
    /* all memory that is necessary for the factorization.                  */
    /* -------------------------------------------------------------------- */


    phase = 11;
    cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase,
        &N, A.data(), iA.data(), jA.data(), &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &myComm, &error);

    if (error != 0)
        OCP_ABORT("ERROR during symbolic factorization: " + to_string((long long int)error));

    /* -------------------------------------------------------------------- */
    /* .. Numerical factorization.                                          */
    /* -------------------------------------------------------------------- */
    phase = 22;
    cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase,
        &N, A.data(), iA.data(), jA.data(), &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &myComm, &error);

    if (error != 0)
        OCP_ABORT("ERROR during numerical factorization: " + to_string((long long int)error));

    /* -------------------------------------------------------------------- */
    /* .. Back substitution and iterative refinement.                       */
    /* -------------------------------------------------------------------- */
    phase = 33;
    cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase,
        &N, A.data(), iA.data(), jA.data(), &idum, &nrhs, iparm, &msglvl, b, x, &myComm, &error);

    if (error != 0)
        OCP_ABORT("ERROR during numerical solution: " + to_string((long long int)error));

    /* -------------------------------------------------------------------- */
    /* .. Termination and release of memory. */
    /* -------------------------------------------------------------------- */
    phase = -1; /* Release internal memory. */
    cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase,
        &N, &ddum, iA.data(), jA.data(), &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &myComm, &error);

    if (error != 0)
        OCP_ABORT("ERROR during release memory: " + to_string((long long int)error));

    return 1;
}


/// Allocate memoery for pardiso solver
void VectorPardisoSolver::Allocate(const OCP_USI& max_nnz, const OCP_USI& maxDim)
{
    iA.resize(maxDim * blockdim + 1);
    jA.resize(max_nnz * blockdim * blockdim);
    A.resize(max_nnz * blockdim * blockdim);
}


/// Calculate terms used in communication
void VectorPardisoSolver::CalCommTerm(const USI& actWellNum, const Domain* domain)
{
    global_index = domain->CalGlobalIndex(actWellNum);

    const OCP_INT numGridInterior = domain->GetNumGridInterior();
    const OCP_INT numElementLoc   = actWellNum + numGridInterior;
    const OCP_INT global_end      = global_index->at(numElementLoc - 1);

    iparm[40] = (global_end - numElementLoc + 1) * blockdim;  // global begin (included)
    iparm[41] = (global_end + 1) * blockdim - 1;              // global end   (included)

    // Get Dimension
    N = (global_end + 1) * blockdim;

    GetWallTime timer;
    timer.Start();

    MPI_Bcast(&N, 1, MPI_INT, domain->numproc - 1, domain->myComm);

    OCPTIME_COMM_COLLECTIVE += timer.Stop() / TIME_S2MS;
}


/// Assemble coefficient matrix.
void VectorPardisoSolver::AssembleMat(const vector<vector<USI>>& colId,
    const vector<vector<OCP_DBL>>& val,
    const OCP_USI& dim,
    vector<OCP_DBL>& rhs,
    vector<OCP_DBL>& u)
{
    const USI blockSize = blockdim * blockdim;
    vector<USI> tmp;
    // Assemble iA, jA, A
    iA[0] = 0;
    for (OCP_USI i = 1; i < dim + 1; i++) {
        const USI nnzR = colId[i - 1].size();
        const OCP_USI bId = (i - 1) * blockdim;

        for (USI c = 0; c < blockdim; c++)
            iA[bId + c + 1] = iA[bId + c] + nnzR * blockdim;
       
        // reorder
        tmp = colId[i - 1];
        for (auto& t : tmp)    t = global_index->at(t);
        sort(tmp.begin(), tmp.end());
        for (USI j0 = 0; j0 < nnzR; j0++) {
            for (USI j1 = 0; j1 < nnzR; j1++) {
                if (global_index->at(colId[i - 1][j1]) == tmp[j0]) {

                    for (USI c = 0; c < blockdim; c++) {
                        for (USI c1 = 0; c1 < blockdim; c1++)
                            jA[iA[bId + c] + j0 * blockdim + c1] = tmp[j0] * blockdim + c1;

                        const OCP_DBL* begin = &val[i - 1][0] + j1 * blockSize + c * blockdim;
                        const OCP_DBL* end = begin + blockdim;
                        copy(begin, end, &A[iA[bId + c] + j0 * blockdim]);
                    }
                        
                    break;
                }
            }
        }
    }

    b = rhs.data();
    x = u.data();
}

#endif // OCPFLOATTYPEWIDTH == 64
#endif // WITH_PARDISO

 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Mar/30/2023      Create file                          */
 /*----------------------------------------------------------------------------*/