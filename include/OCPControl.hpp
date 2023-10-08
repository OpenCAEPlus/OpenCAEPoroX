/*! \file    OCPControl.hpp
 *  \brief   OCPControl class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPCONTROL_HEADER__
#define __OCPCONTROL_HEADER__

// Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "ParamControl.hpp"
#include "Reservoir.hpp"

using namespace std;

/// continue simulating
const OCP_INT OCP_CONTINUE          = 0;  
/// Reset to last time step
const OCP_INT OCP_RESET             = -1; 
/// Reset with cut time(because of failed newton iterations)
const OCP_INT OCP_RESET_CUTTIME     = -2; 
/// Reset with cut time(because of out-ranged cfl number)
const OCP_INT OCP_RESET_CUTTIME_CFL = -3;   


/// shortcut instructions from the command line
class FastControl
{
public:
    FastControl(const USI& argc, const char* optset[]);

public:
    /// If use fastcontrol
    OCP_BOOL ifUse{ OCP_FALSE };
    /// IMPEC, FIM or AIM
    USI      method;
    /// length of the first time step beginning the next TSTEP
    OCP_DBL  timeInit;
    /// Maximum time step during running
    OCP_DBL  timeMax;
    /// Minimum time step during running
    OCP_DBL  timeMin;
    /// Print level
    USI      printLevel{ 0 };
};


/// Record iteration information
class ItersInfo
{
    friend class OCPControl;

public:
    /// Update all Iters
    void UpdateTotal();
    /// Reset all Iters
    void Reset();
    /// Update the number of Newton iterations.
    void UpdateNR() { NR++; }
    /// Update the number of linear iterations.
    void UpdateLS(const USI& num) { LS += num; }
    /// Get num of time steps
    auto GetTimeStep() const { return numTstep; }
    /// Return the number of Newton iterations in one time step.
    auto GetNR() const { return NR; }
    /// Return the total number of Newton iterations.
    auto GetNRt() const { return NRt; }
    /// Return the total number of wasted Newton iterations.
    auto GetNRwt() const { return NRwt; }
    /// Return the number of linear iterations in one time step.
    auto GetLS() const { return LS; }
    /// Return the total number of linear iterations.
    auto GetLSt() const { return LSt; }
    /// Return the total number of wasted linear iterations.
    auto GetLSwt() const { return LSwt; }

protected:
    /// number of time step
    USI numTstep{ 0 };
    /// number of Newton-Raphson iterations
    USI NR{ 0 };
    /// total number of Newton iterations
    USI NRt{ 0 };
    /// total number of wasted Newton iterations
    USI NRwt{ 0 };
    /// number of linear solver iterations
    USI LS{ 0 };
    /// total number of iterations of linear solver
    USI LSt{ 0 };
    /// totalnumber of wasted linear iterations
    USI LSwt{ 0 };
};



/// Params for choosing time stepsize in time marching.
/// Note: Most commonly used params are the first three
class ControlTimeParam
{
    friend class ControlTime;

public:
    /// Default constructor
    ControlTimeParam() = default;
    /// Set only control params
    ControlTimeParam(const TuningPair& src, const vector<OCP_DBL>& Tstep, const USI& i);
    /// Setup for FastControl
    void SetFastControl(const FastControl& fCtrl);

protected:

    /// length of the first time step beginning the next TSTEP
    OCP_DBL timeInit; 
    /// Maximum time step during running
    OCP_DBL timeMax; 
    /// Minimum time step during running
    OCP_DBL timeMin; 
    /// Maximum timestep increase factor
    OCP_DBL maxIncreFac; 
    /// Minimum timestep cutback factor
    OCP_DBL minChopFac;
    /// cutback factor after a convergence failure
    OCP_DBL cutFacNR;    

    // Params for time step prediction, i.e. Limits for changes at next time step
    /// Ideal max Pressure change
    OCP_DBL dPlim;
    /// Ideal max Temperature change
    OCP_DBL dTlim;
    /// Ideal max Saturation change
    OCP_DBL dSlim;
    /// Ideal max relative Ni (moles of components) change
    OCP_DBL dNlim;
    /// Ideal max relative Verr (pore - fluid) change
    OCP_DBL eVlim;

    /// num of TSTEP interval
    USI     numTstepI;
    /// total simulation time
    OCP_DBL total_time;
    /// Begin of TSTEP interval
    OCP_DBL begin_time;
    /// End of TSTEP interval
    OCP_DBL end_time;
};


/// OCP time controler
class ControlTime
{
public:
    /// Setup communicator
    void SetupComm(const Domain& domain) {
        myComm  = domain.myComm;
        numproc = domain.numproc;
        myrank  = domain.myrank;
    }
    void SetCtrlParam(const TuningPair& src, const vector<OCP_DBL>& Tstep, const USI& i) {
        ps.push_back(ControlTimeParam(src, Tstep, i));
    }
    void SetFastControl(const FastControl& fCtrl) {
        for (auto& p : ps)  p.SetFastControl(fCtrl);
    }
    void SetNextTSTEP(const USI& i, const AllWells& wells);
    /// Calculate next time step
    void CalNextTimeStep(const Reservoir& rs, const ItersInfo& iters, const initializer_list<string>& il);
    /// Get total simulation time
    auto GetTotalTime() const { return ps.front().total_time; }
    /// Get number of TSTEP interval
    auto GetNumTstepInterval() const { return ps.front().numTstepI; }

protected:
    MPI_Comm         myComm;
    OCP_INT          numproc, myrank;

public:
    /// cut time
    void CutDt(const OCP_DBL& fac = -1) {
        if (fac < 0) current_dt *= ps[w].cutFacNR;
        else         current_dt *= fac;
    }
    /// Return the current time.
    auto GetCurrentTime() const { return current_time; }
    /// Return current time step size.
    auto GetCurrentDt() const { return current_dt; }
    /// Return last time step size.
    auto GetLastDt() const { return last_dt; }
    /// Determine whether the critical time point has been reached.
    auto IfEnd() { return ((ps[w].end_time - current_time) < TINY); }

protected:
    /// control param set
    vector<ControlTimeParam> ps;
    /// work index
    USI                      w{0};
    /// from prediction for next TSTEP
    OCP_DBL                  predict_dt;
    /// current time step
    OCP_DBL                  current_dt;
    /// last time step
    OCP_DBL                  last_dt;
    /// current time
    OCP_DBL                  current_time{ 0 };
};


/// Params for Newton iterations controls
/// Note: Important for convergence of solution methods
class ControlNR
{
    friend class OCPControl;
public:
    ControlNR() = default;
    ControlNR(const vector<OCP_DBL>& src);
    auto Tol() const { return tol; }
    auto MaxIter() const { return maxIter; }
    auto DPmax() const { return dPmax; }
    auto DSmax() const { return dSmax; }
    auto DPmin() const { return dPmin; }
    auto DSmin() const { return dSmin; }

protected:
    /// Maximum number of Newton iterations in a time step
    USI     maxIter; 
    /// Maximum non-linear convergence error
    OCP_DBL tol;
    /// Maximum Pressure change in a Newton iteration
    OCP_DBL dPmax;
    /// Maximum Saturation change in a Newton iteration
    OCP_DBL dSmax;
    /// Minimum Pressure change in a Newton iteration
    OCP_DBL dPmin;
    /// Minimum Saturation change in a Newton iteration
    OCP_DBL dSmin;
    /// Maximum Verr (vol error b/w fluid and pore) in a Newton step
    OCP_DBL Verrmax;   
};


/// All parameters used for solution control
class OCPControl
{
public:
    /// Input parameters for control.
    void InputParam(const ParamControl& CtrlParam);
    /// Setup Comm
    void Setup(const Domain& domain);
    /// Apply control for time step i.
    void ApplyControl(const USI& i, const Reservoir& rs);
    // Check order is important
    OCP_BOOL Check(Reservoir& rs, const initializer_list<string>& il);
    // Calculate next time step
    void CalNextTimeStep(const Reservoir& rs, const initializer_list<string>& il) {
        time.CalNextTimeStep(rs, iters, il);
    }

public:
    MPI_Comm         myComm;
    OCP_INT          numproc, myrank;
    OCP_INT          workState;           ///< work state of global process
    OCP_INT          workState_loc;       ///< work state of current process

public:
    /// Print level
    USI              printLevel{0};
    /// num of time steps, nonlinear iterations and linear solver iterations
    ItersInfo        iters;
    /// Time control 
    ControlTime      time;
    /// NR control    
    ControlNR        NR;


public:  
    /// Get model
    auto GetModel() const { return model; }
    /// Get type of the solution method.
    auto GetMethod() const { return method; }
    /// Get work dir name.
    auto GetWorkDir() const { return workDir; }
    /// Get OCP file name.
    auto GetOCPFile() const { return ocpFile; }
    /// Get linear solver file name.
    auto GetLsFile() const { return lsFile; }
    /// Setup fast Control.
    void SetupFastControl(const USI& argc, const char* optset[]);

protected:
    /// model: isothermal, thermal
    OCPModel model{ OCPModel::none };
    /// Discrete method
    USI      method;
    /// Current work directory
    string   workDir;
    /// Current file name
    string   ocpFile;
    /// File name of linear Solver
    string   lsFile;
    /// NR control set   
    vector<ControlNR>   ctrlNRSet;
};

#endif /* end if __OCPControl_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/08/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/