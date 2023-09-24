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

#ifndef __OCP_CONTROL_HEADER__
#define __OCP_CONTROL_HEADER__

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


/// Params for choosing time stepsize in time marching.
/// Note: Most commonly used params are the first three
class ControlTime
{
public:
    ControlTime() = default;
    ControlTime(const vector<OCP_DBL>& src_t, const vector<OCP_DBL>& src_pt);

public:
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
};


/// Params for Newton iterations controls
/// Note: Important for convergence of solution methods
class ControlNR
{
public:
    ControlNR() = default;
    ControlNR(const vector<OCP_DBL>& src);

public:
    /// Maximum number of Newton iterations in a time step
    USI     maxNRiter; 
    /// Maximum non-linear convergence error
    OCP_DBL NRtol;
    /// Maximum Pressure change in a Newton iteration
    OCP_DBL NRdPmax;
    /// Maximum Saturation change in a Newton iteration
    OCP_DBL NRdSmax;
    /// Minimum Pressure change in a Newton iteration
    OCP_DBL NRdPmin;
    /// Minimum Saturation change in a Newton iteration
    OCP_DBL NRdSmin;
    /// Maximum Verr (vol error b/w fluid and pore) in a Newton step
    OCP_DBL Verrmax;   
};


/// shortcut instructions from the command line
class FastControl
{
public:
    void ReadParam(const USI& argc, const char* optset[]);

public:
    /// If use fastcontrol
    OCP_BOOL ifUse{OCP_FALSE};
    /// IMPEC, FIM or AIM
    USI      method; 
    /// length of the first time step beginning the next TSTEP
    OCP_DBL  timeInit;  
    /// Maximum time step during running
    OCP_DBL  timeMax;  
    /// Minimum time step during running
    OCP_DBL  timeMin;   
    /// Print level
    USI      printLevel{0};
};


/// Record iteration information
class ItersInfo
{
public:
    /// number of time step
    USI numTstep{ 0 };  
    /// number of linear solver iterations
    USI LS{ 0 };  
    /// total number of iterations of linear solver
    USI LSt{ 0 };
    /// totalnumber of wasted linear iterations
    USI LSwt{ 0 };
    /// number of Newton-Raphson iterations
    USI NR{ 0 }; 
    /// total number of Newton iterations
    USI NRt{ 0 }; 
    /// total number of wasted Newton iterations
    USI NRwt{ 0 };
};


/// All parameters used for solution control
class OCPControl
{
    friend class OpenCAEPoroX;
    friend class OCPOutput;
    friend class Out4RPT;

    friend class IsoT_FIM;
    friend class IsoT_IMPEC;
    friend class IsoT_AIMc;
    friend class T_FIM;
    // temp
    friend class Solver;

public:
    /// Input parameters for control.
    void InputParam(const ParamControl& CtrlParam);

    /// Setup Comm
    void Setup(const Domain& domain);

    /// Get model
    auto GetModel() const { return model; }

    /// Apply control for time step i.
    void ApplyControl(const USI& i, const Reservoir& rs);

    /// Setup fast Control.
    void SetupFastControl(const USI& argc, const char* optset[]);

    /// Return type of the solution method.
    USI GetMethod() const { return method; }

    /// Return number of TSTEPs.
    USI GetNumTSteps() const { return criticalTime.size(); }

    /// Return the current time.
    OCP_DBL GetCurTime() const { return current_time; }

    /// Return current time step size.
    OCP_DBL GetCurDt() const { return current_dt; }

    /// Return last time step size.
    OCP_DBL GetLastDt() const { return last_dt; }

    /// Return the number of linear iterations in one time step.
    USI GetLSiter() const { return iters.LS; }

    /// Return the total number of linear iterations.
    USI GetLSiterT() const { return iters.LSt; }

    /// Return the number of Newton iterations in one time step.
    USI GetNRiter() const { return iters.NR; }

    /// Return the total number of Newton iterations.
    USI GetNRiterT() const { return iters.NRt; }

    /// Update the number of iterations.
    void UpdateIters();

    /// Update the number of linear iterations.
    void UpdateIterLS(const USI& num) { iters.LS += num; }

    /// Update the number of Newton iterations.
    void UpdateIterNR() { iters.NR++; }

    /// Reset the number of iterations.
    void ResetIterNRLS();

    /// Determine whether the critical time point has been reached.
    OCP_BOOL IsCriticalTime(const USI& d)
    {
        return ((criticalTime[d] - current_time) < TINY);
    }

    /// Return work dir name.
    string GetWorkDir() const { return workDir; }

    /// Return linear solver file name.
    string GetLsFile() const { return lsFile; }

    // Check order is important
    OCP_BOOL Check(Reservoir& rs, initializer_list<string> il);

    // Calculate next time step
    void CalNextTimeStep(Reservoir& rs, initializer_list<string> il);

protected:
    MPI_Comm         myComm;
    OCP_INT          numproc, myrank;
    OCP_INT          workState;           ///< work state of global process
    OCP_INT          workState_loc;       ///< work state of current process

protected:
    /// model: ifThermal, isothermal
    OCPModel model{ OCPModel::none };
    /// Discrete method
    USI      method;  
    /// Current work directory
    string   workDir; 
    /// Current file name
    string   fileName;  
    /// File name of linear Solver
    string   lsFile; 

    vector<OCP_DBL> criticalTime; ///< Set of Critical time by user

    // Record time information
    OCP_DBL predict_dt;      ///< from prediction for next TSTEP
    OCP_DBL current_dt;      ///< Current time step
    OCP_DBL last_dt;         ///< last time step
    OCP_DBL current_time{0}; ///< Current time
    OCP_DBL end_time;        ///< Next Critical time

    // Print level
    USI printLevel{0};

    ItersInfo           iters;

    /// Time control
    ControlTime         ctrlTime;
    /// Time control set 
    vector<ControlTime> ctrlTimeSet;
    /// NR control       
    ControlNR           ctrlNR;
    /// NR control set   
    vector<ControlNR>   ctrlNRSet;

    // Receive directly from command lines, which will overwrite others
    FastControl         ctrlFast;
};

#endif /* end if __OCP_Control_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/08/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/