/*! \file    SamgSolver.hpp
 *  \brief   API of SANG Solver
 *  \author  Shizhe Li
 *  \date    Apr/07/2023
 *
 *  \note    use SAMG as linear solver
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "SamgSolver.hpp"


/// Set parameters.
void SamgSolver::SetupParam(const string& dir, const string& file)
{

}


/// Initialize the Params for linear solver.
void SamgSolver::InitParam()
{

}


/// Allocate memoery for pardiso solver
void SamgSolver::Allocate(const OCP_USI& max_nnz,
    const OCP_USI& maxDim,
    const USI& blockDim)
{
    blockdim = blockDim;
    iA.resize(maxDim * blockdim + 1);
    jA.resize(max_nnz * blockdim * blockdim);
    A.resize(max_nnz * blockdim * blockdim);
}


/// Calculate terms used in communication
void SamgSolver::CalCommTerm(const USI& actWellNum, const Domain* domain)
{

}


/// Assemble coefficient matrix.
void ScalarSamgSolver::AssembleMat(const vector<vector<USI>>& colId,
    const vector<vector<OCP_DBL>>& val,
    const OCP_USI& dim,
    vector<OCP_DBL>& rhs,
    vector<OCP_DBL>& u)
{

    b = rhs.data();
    x = u.data();

    OCP_USI bId = 0;
    iA[0] = 0;
    for (OCP_USI i = 1; i < dim + 1; i++) {
        USI nnz_Row = colId[i - 1].size();
        iA[i] = iA[i - 1] + nnz_Row;

        copy(colId[i - 1].begin(), colId[i - 1].end(), &jA[bId]);
        // the first entry is in diagnal line all right
        if (val[i - 1][0] > 0) {           
            copy(val[i - 1].begin(), val[i - 1].end(), &A[bId]);
            bId += nnz_Row;
        }
        else {
            for (USI j = 0; j < nnz_Row; j++) {
                A[bId]  = -val[i - 1][j];
                bId++;
            }
            b[i - 1] = -b[i - 1];
        }      
    }
}


/// Solve the linear system.
OCP_INT SamgSolver::Solve()
{

}


/// Assemble coefficient matrix.
void VectorSamgSolver::AssembleMat(const vector<vector<USI>>& colId,
    const vector<vector<OCP_DBL>>& val,
    const OCP_USI& dim,
    vector<OCP_DBL>& rhs,
    vector<OCP_DBL>& u)
{

    b = rhs.data();
    x = u.data();

    const USI blockSize = blockdim * blockdim;
    // Assemble iA, jA, A
    iA[0] = 0;
    for (OCP_USI i = 1; i < dim + 1; i++) {
        const USI nnzR    = colId[i - 1].size();
        const OCP_USI bId = (i - 1) * blockdim;

        // iA
        for (USI c = 0; c < blockdim; c++)
            iA[bId + c + 1] = iA[bId + c] + nnzR * blockdim;

        // iA and A
        // copy each block
        for (USI j = 0; j < nnzR; j++) {
            for (USI c = 0; c < blockdim; c++) {
                for (USI c1 = 0; c1 < blockdim; c1++)
                    jA[iA[bId + c] + j * blockdim + c1] = colId[i - 1][j] * blockdim + c1;

                const OCP_DBL* begin = &val[i - 1][0] + j * blockSize + c * blockdim;
                const OCP_DBL* end = begin + blockdim;
                copy(begin, end, &A[iA[bId + c] + j * blockdim]);
            }
        }
        // the first entry should be in diagnal line and positive
        for (USI c = 0; c < blockdim; c++) {
            swap(jA[iA[bId + c]], jA[iA[bId + c] + c]);
            swap(A[iA[bId + c]], A[iA[bId + c] + c]);
            if (A[iA[bId + c]] < 0) {
                for (OCP_USI j = iA[bId + c]; j < iA[bId + c + 1]; j++)
                    A[j] = -A[j];
                b[bId + c] = -b[bId + c];
            }
        }
    }

    if (false) {
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        ofstream myFile;
        myFile.open("test/test" + to_string(myrank) + ".out");
        ios::sync_with_stdio(false);
        myFile.tie(0);

        myFile << dim * blockdim << endl;
        for (OCP_USI i = 0; i <= dim * blockdim; i++) {
            myFile << iA[i] << endl;
        }
        for (OCP_USI i = 0; i < dim * blockdim; i++) {
            for (OCP_USI j = iA[i]; j < iA[i + 1]; j++) {
                myFile << jA[j] << endl;
            }
        }
        for (OCP_USI i = 0; i < dim * blockdim; i++) {
            for (OCP_USI j = iA[i]; j < iA[i + 1]; j++) {
                myFile << A[j] << endl;
            }
        }

        myFile.close();
    }
}




 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Apr/07/2023      Create file                          */
 /*----------------------------------------------------------------------------*/