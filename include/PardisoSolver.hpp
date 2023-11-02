/*! \file    PardisoSolver.hpp
 *  \brief   PardisoSolver class declaration
 *  \author  Shizhe Li
 *  \date    Mar/30/2023
 *
 *  \note    The params used in OpenCAEPoroX is mostly compatible with Eclipse by SLB,
 *           but it has some own rules for easy to use. It is extensible and friendly.
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifdef WITH_PARDISO


#ifndef __PARDISOSOLVER_HEADER__
#define __PARDISOSOLVER_HEADER__

#include <math.h>
#include <mpi.h>
#include <mkl.h>
#include <mkl_cluster_sparse_solver.h>
#include <algorithm>

#include "LinearSolver.hpp"

#if    OCPFLOATTYPEWIDTH == 64


#ifdef MKL_ILP64
#define MPI_DT MPI_LONG
#else
#define MPI_DT MPI_INT
#endif
#define MPI_REDUCE_AND_BCAST \
        MPI_Reduce(&err_mem, &error, 1, MPI_DT, MPI_SUM, 0, MPI_COMM_WORLD); \
        MPI_Bcast(&error, 1, MPI_DT, 0, MPI_COMM_WORLD);

using namespace std;

// A unify API for csr and bsr matrix(bsr doesn't work now)
class PardisoSolver : public LinearSolver
{
public:
    PardisoSolver() = default;
    PardisoSolver(const USI& blockDim) { blockdim = blockDim; }

    /// Set parameters.
    void SetupParam(const string& dir, const string& file) override;

    /// Initialize the Params for linear solver.
    void InitParam() override;

    /// Allocate memoery for pardiso solver
    void Allocate(const OCP_USI& max_nnz, const OCP_USI& maxDim) override;

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
    void*     pt[64]    = { 0 };  ///< Internal solver memory pointer pt, don't modify it.
    MKL_INT   iparm[64] = { 0 };  ///< Cluster Sparse Solver control parameters.
    MKL_INT   nrhs      = 1;      ///< Number of right hand sides.
    MKL_INT   maxfct    = 1;      ///< Maximum number of numerical factorizations.
    MKL_INT   mnum      = 1;      ///< Which factorization to use.
    MKL_INT   mtype     = 11;     ///< Defines the matrix type, which influences the pivoting method.
    MKL_INT   phase;              ///< Controls the execution of the solver.
    MKL_INT   N;                  ///< Dimension of Matrix
    MKL_INT   msglvl    = 0;      ///< Message level information.
    MKL_INT   error     = 0;
    MKL_INT   err_mem   = 0;
    double    ddum;               ///< Double dummy
    MKL_INT   idum;               ///< Integer dummy

    // CSR/BSR Matrix
    MKL_INT          blockdim;
    vector<MKL_INT>  iA;
    vector<MKL_INT>  jA;
    vector<double>   A;
    double*          b = nullptr;
    double*          x = nullptr;

    // conmunication
    int              myComm = MPI_Comm_c2f(MPI_COMM_WORLD);
    const vector<OCP_ULL>* global_index;
};

// For unknown reasons, pardiso's BSR version doesn't work, so convert a BSR mat to CSR first
class VectorPardisoSolver : public PardisoSolver
{
public:

    VectorPardisoSolver(const USI& blockDim) { blockdim = blockDim; }

    /// Allocate memoery for pardiso solver
    void Allocate(const OCP_USI& max_nnz, const OCP_USI& maxDim) override;

    void CalCommTerm(const USI& actWellNum, const Domain* domain) override;

    /// Assemble coefficient matrix.
    void AssembleMat(const vector<vector<USI>>& colId,
        const vector<vector<OCP_DBL>>& val,
        const OCP_USI& dim,
        vector<OCP_DBL>& rhs,
        vector<OCP_DBL>& u) override;

};


#endif // OCPFLOATTYPEWIDTH == 64
#endif // __PARDISOSOLVER_HEADER__ 
#endif // WITH_PARDISO

 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Mar/30/2023      Create file                          */
 /*----------------------------------------------------------------------------*/