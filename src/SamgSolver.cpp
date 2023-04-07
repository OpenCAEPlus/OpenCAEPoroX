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

}


/// Calculate terms used in communication
void SamgSolver::CalCommTerm(const USI& actWellNum, const Domain* domain)
{

}


/// Assemble coefficient matrix.
void SamgSolver::AssembleMat(const vector<vector<USI>>& colId,
    const vector<vector<OCP_DBL>>& val,
    const OCP_USI& dim,
    vector<OCP_DBL>& rhs,
    vector<OCP_DBL>& u)
{

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