/*! \file    OCPMatrix.hpp
 *  \brief   OCPMatrix class declaration
 *  \author  Shizhe Li
 *  \date    Dec/06/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPMATRIX_HEADER__
#define __OCPMATRIX_HEADER__

 // Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "Domain.hpp"

using namespace std;


/// Internal maxtirx format in OCP
class OCPMatrix
{

public:
    /// allocate memory
    void Allocate(const Domain& domain, const USI& blockdim);
    /// clear data not free memory
    void ClearData();
    /// Set block dim and block size
    void SetBlockDim(const USI& blockdim);

protected:
    /// allocate memory for row
    void AllocateRowMem(const Domain& domain);
    /// allocate memory for column
    void AllocateColMem(const Domain& domain);

public:
    /// Setup dimensions.
    OCP_USI AddDim(const OCP_USI& n);
    /// add a new value in diagnal line for the first time
    void NewDiag(const OCP_USI& n, const OCP_DBL& v);
    void NewDiag(const OCP_USI& n, const vector<OCP_DBL>& v);
    /// add a vaule in diagnal line 
    void AddDiag(const OCP_USI& n, const OCP_DBL& v);
    void AddDiag(const OCP_USI& n, const vector<OCP_DBL>& v);
    /// Push back a off-diagonal value.
    void NewOffDiag(const OCP_USI& bId, const OCP_USI& eId, const OCP_DBL& v);
    void NewOffDiag(const OCP_USI& bId, const OCP_USI& eId, const vector<OCP_DBL>& v);
    /// Add a value at rhs
    void AddRhs(const OCP_USI& n, const OCP_DBL& v);
    void AddRhs(const OCP_USI& n, const vector<OCP_DBL>& v);
    /// copy rhs to b
    void CopyRhs(const vector<OCP_DBL>& rhs);
    /// set a guess
    void SetGuess(const OCP_USI& n, const OCP_DBL& v) { u[n] = v; }
    /// return the solution
    auto& GetSolution() { return u; }

public:
    /// output A and b to files
    void OutputLinearSystem(const string& dir, const string& fileA, const string& fileb) const;
    /// output solution to file
    void OutputSolution(const string& dir, const string& fileU) const;

public:
    /// check if there are nans in linear system
    void CheckLinearSystem() const;
    /// check if there are nans in solution
    void CheckSolution() const;

public:
    /// Dimens of small block matrix.
    USI nb;
    /// Size of small block matrix.
    USI nb2;

    /// Maximal possible dimension of matrix.
    OCP_USI maxDim;
    /// Actual dimension of matrix.
    OCP_USI dim;
    /// maximum nnz
    OCP_USI max_nnz;

    /// Column indices of nonzero entry.
    vector<vector<OCP_USI>> colId;
    /// Nonzero values.
    vector<vector<OCP_DBL>> val;
    /// Right-hand side of linear system.
    vector<OCP_DBL>         b;
    /// Solution of linear system.
    vector<OCP_DBL>         u;
};







#endif /* end if __OCPMATRIX_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/06/2023      Create file                          */
/*----------------------------------------------------------------------------*/