/*! \file    LinearSystem.hpp
 *  \brief   Linear solver class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __LINEARSYSTEM_HEADER__
#define __LINEARSYSTEM_HEADER__

// Standard header files
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "OCPControlMethod.hpp"
#include "OCPMatrix.hpp"
#include "PardisoSolver.hpp"
#include "SamgSolver.hpp"
#include "FaspSolver.hpp"
#include "PetscSolver.hpp"
#include "LinearSolver.hpp"

using namespace std;


/// Linear solvers for discrete systems.
//  Note: The matrix is stored in the form of row-segmented CSR internally
class LinearSystem
{
public:
    /// Setup linear system
    USI Setup(const OCPModel& model, const string& dir, const string& file, const Domain& d, const USI& nb);
    /// Set work LS
    void SetWorkLS(const USI& i);
    /// Clear the internal matrix data for scalar-value problems.
    void ClearData() { mat.ClearData(); }
    /// Calculate Global start
    void CalCommTerm(const USI& actWellNum) { LS[wIndex]->CalCommTerm(actWellNum, domain); }
    /// Assemble Mat for Linear Solver.
    void AssembleMatLinearSolver();
    /// Solve the Linear System.
    OCP_INT Solve();

protected:
    /// Setup LinearSolver.
    void SetupLinearSolver(const OCPModel& model, const string& lsFile);

public:
    /// Setup dimensions.
    OCP_USI AddDim(const OCP_USI& n) { return mat.AddDim(n); }
    /// Push back a diagonal val, which is always at the first location.
    void NewDiag(const OCP_USI& n, const OCP_DBL& v) { mat.NewDiag(n, v); }
    void NewDiag(const OCP_USI& n, const vector<OCP_DBL>& v) { mat.NewDiag(n, v); }
    /// Add a value at diagonal value.
    void AddDiag(const OCP_USI& n, const OCP_DBL& v) { mat.AddDiag(n, v); }
    void AddDiag(const OCP_USI& n, const vector<OCP_DBL>& v) { mat.AddDiag(n, v); }
    /// Push back a off-diagonal value.
    void NewOffDiag(const OCP_USI& bId, const OCP_USI& eId, const OCP_DBL& v) { mat.NewOffDiag(bId, eId, v); }
    void NewOffDiag(const OCP_USI& bId, const OCP_USI& eId, const vector<OCP_DBL>& v) { mat.NewOffDiag(bId, eId, v); }
    /// Add a value at b[n].
    // template <typename T>
    void AddRhs(const OCP_USI& n, const OCP_DBL& v) { mat.AddRhs(n, v); }
    void AddRhs(const OCP_USI& n, const vector<OCP_DBL>& v) { mat.AddRhs(n, v); }
    /// copy rhs to b
    void CopyRhs(const vector<OCP_DBL>& rhs) { mat.CopyRhs(rhs); }
    /// Assign an initial value at u[n].
    void SetGuess(const OCP_USI& n, const OCP_DBL& v) { mat.SetGuess(n, v); }
    /// Return the solution.
    vector<OCP_DBL>& GetSolution() { return mat.GetSolution(); }


public:
    /// Output the mat and rhs to fileA and fileb.
    void OutputLinearSystem(const string& fileA, const string& fileb) const { mat.OutputLinearSystem(solveDir, fileA, fileb); }
    /// Output the solution to file
    void OutputSolution(const string& fileX) const { mat.OutputSolution(solveDir, fileX); }

protected:
    /// working index of LS
    USI                   wIndex;  
    /// Current workdir.
    string                solveDir; 
    /// domain decomposition information
    const Domain*         domain;
    /// internal matrix
    OCPMatrix             mat;
    /// block dimension set
    vector<USI>           nbs;
    /// LStype sets
    vector<OCPLStype>     LStype;
    /// LS sets
    vector<LinearSolver*> LS;

};

#endif /* end if __LINEARSOLVER_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Nov/09/2021      Remove decoupling methods            */
/*  Chensong Zhang      Nov/22/2021      renamed to LinearSystem              */
/*----------------------------------------------------------------------------*/