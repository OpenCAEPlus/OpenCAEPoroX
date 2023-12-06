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


#ifdef WITH_SAMG

#ifndef __SAMGSOLVER_HEADER__
#define __SAMGSOLVER_HEADER__

#include "LinearSolver.hpp"
#include "samg.h"

using namespace std;

/// API for SAMG Solver
// Note: the index is 1-baesd
class SamgSolver : public LinearSolver
{
public:

    /// Set parameters.
    void SetupParam(const string& dir, const string& file) override;

    /// Initialize the Params for linear solver.
    void InitParam() override;

    /// Allocate memoery for pardiso solver
    void Allocate(const OCPMatrix& mat) override;

    /// Calculate terms used in communication
    void CalCommTerm(const USI& actWellNum, const Domain* domain) override;

    /// Solve the linear system.
    OCP_INT Solve() override;

    /// Get number of iterations used by iterative solver.
    USI GetNumIters() const override { return 1; }

protected:

    // CSR Mat
    SAMG_INT            nb;
    vector<SAMG_INT>    iA;
    vector<SAMG_INT>    jA;
    vector<SAMG_REAL>   A;

    SAMG_REAL*          b = nullptr;
    SAMG_REAL*          x = nullptr;

    // Physical Info
    SAMG_INT           nsys;         ///< Number of unknowns(physical variable type)
    SAMG_INT           nnu;          ///< number of (local) variables
    SAMG_INT           nna;          ///< number of matrix entries stored in the (local) vector a(local nnz)
    SAMG_INT           ndiu;         ///< size of iu
    SAMG_INT           ndip;         ///< size of ip
    vector<SAMG_INT>   iu;           ///< variable-to-unknown pointer: nnu
    vector<SAMG_INT>   ip;           ///< variable-to-point pointer: nnu
    vector<SAMG_INT>   iu_tmp;       ///< template of iu: nb
    vector<SAMG_INT>   nunknown_description;
    
    // comunication
    SAMG_INT            myComm = (SAMG_INT)MPI_Comm_c2f(MPI_COMM_WORLD);
    SAMG_INT            npsnd;       ///< Total number of neighboring processors to send
    vector<SAMG_INT>    iranksnd;    ///< Rank of neighboring processors to send
    SAMG_INT            nshalo;      ///< Total number of variables to send
    vector<SAMG_INT>    ipts;        ///< Range of variables to send in isndlist
    vector<SAMG_INT>    isndlist;    ///< Index of variables to send
    SAMG_INT            nprec;       ///< Total number of neighboring processors to receive
    vector<SAMG_INT>    irankrec;    ///< Rank of neighboring processors to receive
    SAMG_INT            nrhalo;      ///< Total number of variables to receive
    vector<SAMG_INT>    iptr;        ///< Range of variables to receive in isndlist
    vector<SAMG_INT>    ireclist;    ///< Index of variables to recv

    const vector<OCP_USI>* global_index;

protected:
    // Samg params
    SAMG_INT          noil_approach            = 12;
    SAMG_INT          noil_cyc                 = 19;
    SAMG_REAL         noil_preparation         = 19.4;
    SAMG_INT          ierr                     = 0;      ///< Code number indicating errors or warnings.
    SAMG_INT          samg_matrix              = 220;    ///< Type of the matrix A
    SAMG_INT          ncyc_done                = 0;      ///< Total number of cycles (iterations) performed.
    SAMG_REAL         res_in                   = 0;      ///< Residual of first guess.
    SAMG_REAL         res_out                  = 0;      ///< Residual of final approximation.
    
    SAMG_INT          nsolve                   = 2;      ///< Specifies SAMG¡¯s solution strategy. Suggested choice.
    SAMG_INT          ncyc                     = 13050;  ///< Cycling and acceleration strategy. Suggested choice.
    SAMG_INT          ifirst                   = 1;      ///< Selects first approximation for guess
    SAMG_INT          iswtch                   = 51;     ///< Memory extension switch
    SAMG_INT          iout                     = -1;      ///< Print output during the solution phase
    SAMG_INT          idump                    = -1;      ///< Print output during the setup phase
    SAMG_REAL         eps                      = 1.0e-3; ///< Standard stopping criterion for the AMG iteration
    SAMG_REAL         chktol                   = -1.0;   ///< Checking of input matrix
    SAMG_REAL         a_cmplx                  = 2.5;    ///< used to allocate SAMG¡¯s initial memory
    SAMG_REAL         g_cmplx                  = 1.8;    ///< used to allocate SAMG¡¯s initial memory
    SAMG_REAL         p_cmplx                  = 2.0;    ///< used to allocate SAMG¡¯s initial memory
    SAMG_REAL         w_avrge                  = 2.5;    ///< used to allocate SAMG¡¯s initial memory

};

// Convert Internal mat(csr-like) to CSR mat 
class ScalarSamgSolver : public SamgSolver
{
public:
    ScalarSamgSolver(const OCPModel& model) {
        nb       = 1;
        nsys     = 1;
        ndiu     = 1;
        ndip     = 1;
        ifirst   = 0;   // last solution as initial guess
    }

    /// Assemble coefficient matrix.
    void AssembleMat(OCPMatrix& mat) override;
};

// Convert Internal mat(bsr-like) to CSR mat 
class VectorSamgSolver : public SamgSolver
{
public:
    VectorSamgSolver(const USI& blockDim, const OCPModel& model) {
        nb = blockDim;
        nsys     = nb;
        iu_tmp.resize(nb);
        ifirst   = 1;   // zero solution as initial guess
        for (USI i = 0; i < nb; i++)  iu_tmp[i] = i + 1;
        if (model == OCPModel::isothermal) {
            nunknown_description.resize(nsys, 2);
            nunknown_description[0] = 0;          // Pressure
        }
        else if (model == OCPModel::thermal) {
            nunknown_description.resize(nsys, 2); // Concentration(init)
            nunknown_description.front() = 0;     // Pressure
            nunknown_description.back()  = 100;   // Temperature           
        }
        else                            OCP_ABORT("Wrong Model for SAMG Solver!");              
    }
    /// Assemble coefficient matrix.
    void AssembleMat(OCPMatrix& mat) override;
};

#endif 

#endif // WITH_SAMG


 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Apr/07/2023      Create file                          */
 /*----------------------------------------------------------------------------*/