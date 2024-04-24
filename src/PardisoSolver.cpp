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

#if OCPFLOATTYPEWIDTH == 64


PardisoSolver::PardisoSolver(const string& dir, const string& file, const OCPMatrix& mat)
{
    SetupParam(dir, file);
    Allocate(mat);
}


/// Assemble coefficient matrix.
void PardisoSolver::AssembleMat(OCPMatrix& mat, const Domain* domain)
{
    CalCommTerm(domain);

    const USI blockSize = nb * nb;
    vector<USI> tmp;  
    // Assemble iA, jA, A
    iA[0] = 0;
    for (OCP_USI i = 1; i < mat.dim + 1; i++) {
        const USI nnzR = mat.colId[i - 1].size();

        iA[i] = iA[i - 1] + nnzR;
        // reorder
        tmp = mat.colId[i - 1];
        for (auto& t : tmp)    t = global_index->at(t);
        sort(tmp.begin(), tmp.end());
        for (USI j0 = 0; j0 < nnzR; j0++) {
            jA[iA[i - 1] + j0] = tmp[j0];
            for (USI j1 = 0; j1 < nnzR; j1++) {
                if (global_index->at(mat.colId[i - 1][j1]) == tmp[j0]) {
                    const OCP_DBL* begin = &mat.val[i - 1][0] + j1 * blockSize;
                    const OCP_DBL* end = begin + blockSize;
                    copy(begin, end, &A[(iA[i - 1] + j0) * blockSize]);
                    break;
                }
            }
        }           
    }


    b = mat.b.data();
    x = mat.u.data();
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


void PardisoSolver::SetupParam(const string& dir, const string& file)
{
    InitParam();
}


void PardisoSolver::InitParam()
{
    iparm[0] = 1;  /* Solver default parameters overriden with provided by iparm */
    iparm[1] = 2;  /* Use METIS for fill-in reordering */
    iparm[5] = 0;  /* Write solution into x */
    iparm[7] = 2;  /* Max number of iterative refinement steps */
    iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
    iparm[12] = 1;  /* Switch on Maximum Weighted Matching algorithm (default for non-symmetric) */
    iparm[17] = 0;  /* (closed) Output: Number of nonzeros in the factor LU */
    iparm[18] = 0;  /* (closed) Output: Mflops for LU factorization */
    iparm[26] = 1;  /* Check input data for correctness */
    iparm[39] = 2;  /* Input: matrix/rhs/solution are distributed between MPI processes  */

    /*Zero-based indexing: columns and rows indexing in arrays ia, ja, and perm starts from 0 (C-style indexing).*/
    iparm[34] = 1;
}


void PardisoSolver::Allocate(const OCPMatrix& mat)
{
    nb = mat.nb;
    if (nb > 1) {
        iparm[37] = nb * nb;
    }

    iA.resize(mat.maxDim + 1);
    jA.resize(mat.max_nnz);
    A.resize(mat.max_nnz * nb * nb);
}


/// Calculate terms used in communication
void PardisoSolver::CalCommTerm(const Domain* domain)
{

    global_index = domain->CalGlobalIndex();
    myComm       = MPI_Comm_c2f(domain->cs_comm);

    const OCP_INT numElementLoc = domain->GetNumActElementForSolver();
    const OCP_INT global_end = global_index->at(numElementLoc - 1);

    iparm[40] = global_end - numElementLoc + 1;  // global begin (included)
    iparm[41] = global_end;                      // global end   (included)

    // Get Dimension
    N = global_end + 1;

    GetWallTime timer;
    timer.Start();

    MPI_Bcast(&N, 1, MPI_INT, *domain->cs_group_local_rank.rbegin(), domain->cs_comm);

    OCPTIME_COMM_COLLECTIVE += timer.Stop();
}


VectorPardisoSolver::VectorPardisoSolver(const string& dir, const string& file, const OCPMatrix& mat)
{
    SetupParam(dir, file);
    Allocate(mat);
}


/// Assemble coefficient matrix.
void VectorPardisoSolver::AssembleMat(OCPMatrix& mat, const Domain* domain)
{
    CalCommTerm(domain);

    const USI blockSize = nb * nb;
    vector<USI> tmp;
    // Assemble iA, jA, A
    iA[0] = 0;
    for (OCP_USI i = 1; i < mat.dim + 1; i++) {
        const USI nnzR = mat.colId[i - 1].size();
        const OCP_USI bId = (i - 1) * nb;

        for (USI c = 0; c < nb; c++)
            iA[bId + c + 1] = iA[bId + c] + nnzR * nb;
       
        // reorder
        tmp = mat.colId[i - 1];
        for (auto& t : tmp)    t = global_index->at(t);
        sort(tmp.begin(), tmp.end());
        for (USI j0 = 0; j0 < nnzR; j0++) {
            for (USI j1 = 0; j1 < nnzR; j1++) {
                if (global_index->at(mat.colId[i - 1][j1]) == tmp[j0]) {

                    for (USI c = 0; c < nb; c++) {
                        for (USI c1 = 0; c1 < nb; c1++)
                            jA[iA[bId + c] + j0 * nb + c1] = tmp[j0] * nb + c1;

                        const OCP_DBL* begin = &mat.val[i - 1][0] + j1 * blockSize + c * nb;
                        const OCP_DBL* end = begin + nb;
                        copy(begin, end, &A[iA[bId + c] + j0 * nb]);
                    }
                        
                    break;
                }
            }
        }
    }

    b = mat.b.data();
    x = mat.u.data();
}


/// Allocate memoery for pardiso solver
void VectorPardisoSolver::Allocate(const OCPMatrix& mat)
{
    nb = mat.nb;
    iA.resize(mat.maxDim * nb + 1);
    jA.resize(mat.max_nnz * nb * nb);
    A.resize(mat.max_nnz * nb * nb);
}


/// Calculate terms used in communication
void VectorPardisoSolver::CalCommTerm(const Domain* domain)
{
    global_index = domain->CalGlobalIndex();
    myComm       = MPI_Comm_c2f(domain->cs_comm);

    const OCP_INT numElementLoc = domain->GetNumActElementForSolver();
    const OCP_INT global_end = global_index->at(numElementLoc - 1);

    iparm[40] = (global_end - numElementLoc + 1) * nb;  // global begin (included)
    iparm[41] = (global_end + 1) * nb - 1;              // global end   (included)

    // Get Dimension
    N = (global_end + 1) * nb;

    GetWallTime timer;
    timer.Start();

    MPI_Bcast(&N, 1, MPI_INT, *domain->cs_group_local_rank.rbegin(), domain->cs_comm);

    OCPTIME_COMM_COLLECTIVE += timer.Stop();
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