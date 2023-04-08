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
        // the first entry is in diagnal line
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

}




 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Apr/07/2023      Create file                          */
 /*----------------------------------------------------------------------------*/