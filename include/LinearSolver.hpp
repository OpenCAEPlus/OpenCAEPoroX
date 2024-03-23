/*! \file    LinearSolver.hpp
 *  \brief   LinearSolver class declaration
 *  \author  Shizhe Li
 *  \date    Nov/22/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __LINEARSOLVER_HEADER__
#define __LINEARSOLVER_HEADER__

// Standard header files
#include <string>
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "Domain.hpp"
#include "OCPMatrix.hpp"

using namespace std;


enum class OCPLStype : USI
{
    none,
    fasp,
    pardiso,
    petsc,
    samg
};


/// Virtual base class for linear solvers.
class LinearSolver
{
public:
    /// Assemble matrix for linear solver from the internal matrix data.
    virtual void AssembleMat(OCPMatrix& mat, const Domain* domain) = 0;

    /// Solve the linear system and return the number of iterations.
    virtual OCP_INT Solve() = 0;

    /// Get number of iterations.
    virtual USI GetNumIters() const = 0;
};

#endif // __LINEARSOLVER_HEADER__

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/22/2021      Create file                          */
/*  Chensong Zhang      Jan/18/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/