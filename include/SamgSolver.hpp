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
    void Allocate(const OCP_USI& max_nnz, const OCP_USI& maxDim) override;

    /// Calculate terms used in communication
    void CalCommTerm(const USI& actWellNum, const Domain* domain) override;

    /// Solve the linear system.
    OCP_INT Solve() override;

    /// Get number of iterations used by iterative solver.
    USI GetNumIters() const override { return 1; }

protected:

    // CSR Mat
    int            blockdim;
    vector<int>    iA;
    vector<int>    jA;
    vector<double> A;

    double*        b = nullptr;
    double*        x = nullptr;

    // Physical Info
    int           nsys;         ///< num of systems(physical variable type)
    int           nnu;          ///< number of (local) variables
    int           nna;          ///< number of matrix entries stored in the (local) vector a(local nnz)
    int           ndiu;         ///< size of iu
    int           ndip;         ///< size of ip
    vector<int>   iu;           ///< variable-to-unknown pointer: dim
    vector<int>   ip;           ///< variable-to-point pointer: dim
    vector<int>   iu_tmp;       ///< template of iu: blockdim
    vector<int>   nunknown_description;
    

    // comunication
    int            myComm = MPI_Comm_c2f(MPI_COMM_WORLD);
    int            npsnd;       ///< Total number of neighboring processors to send
    vector<int>    iranksnd;    ///< Rank of neighboring processors to send
    int            nshalo;      ///< Total number of variables to send
    vector<int>    ipts;        ///< Range of variables to send in isndlist
    vector<int>    isndlist;    ///< Index of variables to send
    int            nprec;       ///< Total number of neighboring processors to receive
    vector<int>    irankrec;    ///< Rank of neighboring processors to receive
    int            nrhalo;      ///< Total number of variables to receive
    vector<int>    iptr;        ///< Range of variables to receive in isndlist
    vector<int>    ireclist;    ///< Index of variables to recv

protected:
    // Samg params
    int            noil_approach            = 12;
    int            noil_cyc                 = 19;
    double         noil_preparation         = 19.4;   
    int            ierr                     = 0;
    int            samg_matrix              = 220;
    int            ncyc_done                = 0;
    double         res_in                   = 0;
    double         res_out                  = 0;
    
    int            nsolve                   = 2;      ///< AMG approach
    int            ncyc                     = 11030;  ///< type of cycling. Here: CG/V-Cycle, maixmally 30 iterations
    int            ifirst                   = 1;      ///< which first guess (0=input vector u; 1=zero)
    int            iswtch                   = 51;     ///< various controls; see manual
    int            iout                     = 2;      ///< amount of screen output
    int            idump                    = 0;      ///< dumping of matrices/operators
    double         eps                      = 1.0e-4; ///< <solution tolerance
    double         chktol                   = -1.0e0; ///< switch off internal sanity checks
    double         a_cmplx                  = 2.2e0;  ///< dimension estimate
    double         g_cmplx                  = 1.7e0;  ///< dimension estimate
    double         p_cmplx                  = 0.0e0;  ///< dimension estimate; irrelevant in the nsys=1 case
    double         w_avrge                  = 2.4e0;  ///< dimension estimate

};

// Convert Internal mat(csr-like) to CSR mat 
class ScalarSamgSolver : public SamgSolver
{
public:
    ScalarSamgSolver(const USI& blockDim, const USI& model) {
        blockdim = blockDim;
        nsys     = 1;
        ndiu     = 1;
        ndip     = 1;
    }

    /// Assemble coefficient matrix.
    void AssembleMat(const vector<vector<USI>>& colId,
        const vector<vector<OCP_DBL>>& val,
        const OCP_USI& dim,
        vector<OCP_DBL>& rhs,
        vector<OCP_DBL>& u) override;
};

// Convert Internal mat(bsr-like) to CSR mat 
class VectorSamgSolver : public SamgSolver
{
public:
    VectorSamgSolver(const USI& blockDim, const USI& model) {
        blockdim = blockDim;        
        if (model == ISOTHERMALMODEL) {
            nsys = 2; 
            nunknown_description.resize(nsys);
            nunknown_description[0] = 0;   // Pressure
            nunknown_description[1] = 2;   // Concentration(init)
            iu_tmp.resize(blockdim, 2);    // point to Concentration(init)
            iu_tmp[0] = 1;                 // point to Pressure
        }
        else if (model == THERMALMODEL) {
            nsys = 3;
            nunknown_description.resize(nsys);
            nunknown_description[0] = 0;     // Pressure
            nunknown_description[1] = 2;     // Concentration(init)
            nunknown_description[2] = 100;   // Temperature
            iu_tmp.resize(blockdim, 2);       // point to Concentration(init)
            iu_tmp[0] = 1;                   // point to Pressure
            iu_tmp.front() = 3;              // point to Temperature
        }
        else                            OCP_ABORT("Wrong Model for SAMG Solver!");              
    }
    /// Assemble coefficient matrix.
    void AssembleMat(const vector<vector<USI>>& colId,
        const vector<vector<OCP_DBL>>& val,
        const OCP_USI& dim,
        vector<OCP_DBL>& rhs,
        vector<OCP_DBL>& u) override;
};

#endif // WITH_SAMG

#endif


 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Apr/07/2023      Create file                          */
 /*----------------------------------------------------------------------------*/