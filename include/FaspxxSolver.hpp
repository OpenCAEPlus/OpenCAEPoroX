/*! \file    FaspxxSolver.hpp
 *  \brief   FaspxxSolver class declaration
 *  \author  Li Zhao
 *  \date    April/27/2024
 *
 *  \note    The params used in OpenCAEPoroX is mostly compatible with Eclipse by SLB,
 *           but it has some own rules for easy to use. It is extensible and friendly.
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


#ifdef OCP_USE_FASPXX


#ifndef __FASPXX_HEADER__
#define __FASPXX_HEADER__

#include "LinearSolver.hpp"
#include <vector>

// faspxx header files
extern "C" {

#include "faspxx.h"

}


/// API for Faspxx Solver
// Note: the index is 0-baesd
class FaspxxSolver : public LinearSolver
{
public:
    FaspxxSolver() = default;

    /// Assemble coefficient matrix.
    void AssembleMat(OCPMatrix& mat, const Domain* domain) override;

    /// Get number of iterations used by iterative solver.
    USI GetNumIters() const override { return 1; }

protected:

    /// Allocate memoery for pardiso solver
    void Allocate(const OCPMatrix& mat, const Domain* domain);
    /// Calculate terms used in communication
    void CalCommTerm(const Domain* domain);

protected:

    // CSR/BSR   
    OCP_INT                 nb;         ///< block dim              
    vector<OCP_DBL>         A;          ///< value
    vector<OCP_SLL>         iA;         ///< row ptr
    vector<OCP_SLL>         jA;         ///< col index  
    OCP_DBL*                b;          ///< rhs
    OCP_DBL*                x;          ///< solution

    // Communication
    MPI_Comm                myComm{ MPI_COMM_NULL };  ///< Communicator
    const vector<OCP_ULL>*  global_index;             ///< global index
    vector<OCP_SLL>         allBegin;                 ///< begin for all process(self-include) 
    vector<OCP_SLL>         allEnd;                   ///< end for all process(self-include)
    vector<OCP_INT>         allEle;                   ///< num of elements for each proc
};


class ScalarFaspxxSolver : public FaspxxSolver
{
public:
    ScalarFaspxxSolver(const string& dir, const string& file, const OCPMatrix& mat, const Domain* domain);
    /// Solve the linear system.
    OCP_INT Solve() override;

};

class VectorFaspxxSolver : public FaspxxSolver
{
public:
    VectorFaspxxSolver(const string& dir, const string& file, const OCPMatrix& mat, const Domain* domain);
    /// Solve the linear system.
    OCP_INT Solve() override;
};

#endif // __FASPXX_HEADER__ 
#endif // OCP_USE_FASPXX

 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Li Zhao           April/27/2024      Create file                          */
 /*----------------------------------------------------------------------------*/