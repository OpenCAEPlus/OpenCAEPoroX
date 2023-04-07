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

#ifndef __SAMGSOLVER_HEADER__
#define __SAMGSOLVER_HEADER__

#include "LinearSolver.hpp"

using namespace std;

class SamgSolver : public LinearSolver
{

public:
    /// Set parameters.
    void SetupParam(const string& dir, const string& file) override;

    /// Initialize the Params for linear solver.
    void InitParam() override;

    /// Allocate memoery for pardiso solver
    void Allocate(const OCP_USI& max_nnz,
        const OCP_USI& maxDim,
        const USI& blockDim) override;

    /// Calculate terms used in communication
    void CalCommTerm(const USI& actWellNum, const Domain* domain) override;

    /// Assemble coefficient matrix.
    void AssembleMat(const vector<vector<USI>>& colId,
        const vector<vector<OCP_DBL>>& val,
        const OCP_USI& dim,
        vector<OCP_DBL>& rhs,
        vector<OCP_DBL>& u) override;

    /// Solve the linear system.
    OCP_INT Solve() override;

    /// Get number of iterations used by iterative solver.
    USI GetNumIters() const override { return 1; }

protected:

    
};


class VectorSamgSolver : public SamgSolver
{
public:

    /// Assemble coefficient matrix.
    void AssembleMat(const vector<vector<USI>>& colId,
        const vector<vector<OCP_DBL>>& val,
        const OCP_USI& dim,
        vector<OCP_DBL>& rhs,
        vector<OCP_DBL>& u) override;
};



#endif


 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Apr/07/2023      Create file                          */
 /*----------------------------------------------------------------------------*/