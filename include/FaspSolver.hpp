/*! \file    FaspSolver.hpp
 *  \brief   Declaration of classes interfacing to the FASP solvers
 *  \author  Shizhe Li
 *  \date    Nov/22/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "../config/config.hpp"


#ifdef OCP_USE_FASP

#ifndef __FASPSOLVER_HEADER__
#define __FASPSOLVER_HEADER__


// Standard header files
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// OpenCAEPoroX header files
#include "LinearSolver.hpp"

// faspsolver header files
extern "C" {

#if OCPFLOATTYPEWIDTH == 64

#define SHORT  short
#define INT    int   
#define LONG   long  
#define REAL   double

#endif // OCPFLOATTYPEWIDTH == 64

#if OCPFLOATTYPEWIDTH == 128

#define SHORT  short
#define INT    int   
#define LONG   long  
#define REAL   long double

#endif // OCPFLOATTYPEWIDTH == 128

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

#undef SHORT
#undef INT 
#undef LONG
#undef REAL
}

// faspcpr header files
#if WITH_FASPCPR
extern "C" {
#include "faspcpr.h"
#include "faspcpr_functs.h"
}
#endif

// fasp4blkoil header files
#ifdef OCP_USE_FASP4BLKOIL
extern "C" {

#if OCPFLOATTYPEWIDTH == 64

#define SHORT  short
#define INT    int   
#define LONG   long  
#define REAL   double

#endif // OCPFLOATTYPEWIDTH == 64

#if OCPFLOATTYPEWIDTH == 128

#define SHORT  short
#define INT    int   
#define LONG   long  
#define REAL   long double

#endif // OCPFLOATTYPEWIDTH == 128

#include "fasp4blkoil.h"
#include "fasp4blkoil_functs.h"

#undef SHORT
#undef INT 
#undef LONG
#undef REAL
}
#endif

// fasp4cuda header files
// Note: It should not inside extern "C" {} !
#if WITH_FASP4CUDA
#include "fasp4cuda.h"
#include "fasp4cuda_functs.h"
#endif

using namespace std;

// Standard preconditioner types
#define PC_NULL  60 ///< None:  no preconditioner
#define PC_FASP1 61 ///< FASP1: MSP, default for FIM from 2020
#define PC_FASP2 62 ///< FASP2: MSP, experimental only
#define PC_FASP3 63 ///< FASP3: MSP, monolithic preconditioner
#define PC_FASP4 64 ///< FASP4: MSP, default for FIM from 2015
#define PC_FASP5 65 ///< FASP5: MSP, experimental only
#define PC_DIAG  68 ///< DIAG:  diagonal preconditioner
#define PC_BILU  69 ///< BILU:  block ILU preconditioner

// Sharing-setup preconditioner types
#define PC_FASP1_SHARE 71 ///< Sharing setup stage for PC_FASP1, use with caution
#define PC_FASP4_SHARE 74 ///< Sharing setup stage for PC_FASP4, use with caution
#define RESET_CONST    35 ///< Sharing threshold for PC_FASP1_SHARE, PC_FASP4_SHARE

/// Basic FASP solver class.
class FaspSolver : public LinearSolver
{
public:

    /// Get number of iterations used by iterative solver.
    USI GetNumIters() const override { return itsParam.maxit; }

protected:
    /// Set FASP parameters.
    void SetupParam(const string& dir, const string& file);

    virtual void InitParam() = 0;

public:
    string      solveDir;  ///< Current work dir
    string      solveFile; ///< Relative path of fasp file
    input_param inParam;   ///< Parameters from input files
    ITS_param   itsParam;  ///< Parameters for iterative method
    AMG_param   amgParam;  ///< Parameters for AMG method
    ILU_param   iluParam;  ///< Parameters for ILU method
    SWZ_param   swzParam;  ///< Parameters for Schwarz method
};

/// Scalar solvers in CSR format from FASP.
class ScalarFaspSolver : public FaspSolver
{
public:
    ScalarFaspSolver(const string& dir, const string& file, const OCPMatrix& mat);

    /// Assemble coefficient matrix.
    void AssembleMat(OCPMatrix& mat, const Domain* domain) override;

    /// Solve the linear system.
    OCP_INT Solve() override;

protected:

    /// Initialize the Params for linear solver.
    void InitParam() override;

    /// Allocate memory for the linear system.
    void Allocate(const OCPMatrix& mat);

protected:
    dCSRmat A; ///< Matrix for scalar-value problems
    dvector b; ///< Right-hand side for scalar-value problems
    dvector x; ///< Solution for scalar-value problems
};

/// Vector solvers in BSR format from FASP.
class VectorFaspSolver : public FaspSolver
{
public:
    VectorFaspSolver(const string& dir, const string& file, const OCPMatrix& mat);

    /// Assemble coefficient matrix.
    void AssembleMat(OCPMatrix& mat, const Domain* domain) override;

    /// Solve the linear system.
    OCP_INT Solve() override;

protected:

    /// Initialize the Params for linear solver.
    void InitParam() override;

    /// Allocate memory for the linear system.
    void Allocate(const OCPMatrix& mat);

    /// Apply decoupling to the linear system.
    void Decoupling(dBSRmat* Absr,
                    dvector* b,
                    dBSRmat* Asc,
                    dvector* fsc,
                    ivector* order,
                    OCP_DBL* Dmatvec,
                    int      decouple_type);

protected:
    dBSRmat A; ///< Matrix for vector-value problems
    dvector b; ///< Right-hand side for vector-value problems
    dvector x; ///< Solution for vector-value problems

    dBSRmat Asc;   ///< Scaled matrix for vector-value problems
    dvector fsc;   ///< Scaled right-hand side for vector-value problems
    ivector order; ///< User-defined ordering for smoothing process

    vector<OCP_DBL> Dmat; ///< Decoupling matrices
};

#endif // __FASPSOLVER_HEADER__

#endif // OCP_USE_FASP

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/22/2021      Create file                          */
/*  Chensong Zhang      Jan/08/2022      Update Doxygen                       */
/*  Chensong Zhang      Jan/19/2022      Set FASP4BLKOIL as optional          */
/*  Li Zhao             Apr/04/2022      Set FASP4CUDA   as optional          */
/*----------------------------------------------------------------------------*/
