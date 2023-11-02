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
#include "OCPNRsuite.hpp"

using namespace std;


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
    /// Setup control param
    void SetCtrlParam(const TuningPair& src, const vector<OCP_DBL>& Tstep, const USI& i) {
        ps.push_back(ControlTimeParam(src, Tstep, i));
    }
    /// Setup fast control param
    void SetFastControl(const FastControl& fCtrl) {
        for (auto& p : ps)  p.SetFastControl(fCtrl);
    }
    /// cut time
    void CutDt(const OCP_DBL& fac = -1);
    void CutDt(const OCPNRsuite& NRs);
    /// Set param for next TSTEP
    void SetNextTSTEP(const USI& i, const AllWells& wells);
    /// Calculate next time step
    void CalNextTimeStep(const OCPNRsuite& NRs, const initializer_list<string>& il);
    /// Get total simulation time
    auto GetTotalTime() const { return ps.back().end_time; }
    /// Get number of TSTEP interval
    auto GetNumTstepInterval() const { return ps.size(); }

protected:
    MPI_Comm         myComm;
    OCP_INT          numproc, myrank;

public:

    /// Return the current time.
    auto GetCurrentTime() const { return current_time; }
    /// Return current time step size.
    auto GetCurrentDt() const { return current_dt; }
    /// Return last time step size.
    auto GetLastDt() const { return last_dt; }
    /// Determine whether the critical time point has been reached.
    auto IfEnd() { return ((wp->end_time - current_time) < TINY); }

protected:
    /// control param set
    vector<ControlTimeParam> ps;
    /// work param
    const ControlTimeParam*  wp;
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
class ControlNRParam
{
    friend class ControlNR;
public:
    ControlNRParam() = default;
    ControlNRParam(const vector<OCP_DBL>& src);

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


class ControlNR
{
public:
    /// Setup control param
    void SetCtrlParam(const vector<OCP_DBL>& src) {
        ps.push_back(ControlNRParam(src));
    }
    /// Setup communicator
    void SetupComm(const Domain& domain) {
        myComm = domain.myComm;
        numproc = domain.numproc;
        myrank = domain.myrank;
    }
    /// Set param for next TSTEP
    void SetNextTSTEP(const USI& i) { wp = &ps[i]; }
    /// Get dSmax
    auto DSmax() const { return wp->dSmax; }
    /// Get dPmax
    auto DPmax() const { return wp->dPmax; }
    /// If NR iterations converge
    OCPNRStateC CheckConverge(const OCPNRsuite& NRs, const initializer_list<string>& il) const;

protected:
    MPI_Comm         myComm;
    OCP_INT          numproc, myrank;

protected:
    /// control param set
    vector<ControlNRParam> ps;
    /// current param
    const ControlNRParam*  wp;
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
    /// Check if converge
    OCPNRStateC CheckConverge(const OCPNRsuite& NRs, const initializer_list<string>& il) {
        return NR.CheckConverge(NRs, il);
    }
    // Calculate next time step
    void CalNextTimeStep(const OCPNRsuite& NRs, const initializer_list<string>& il) {
        time.CalNextTimeStep(NRs, il);
    }

public:
    MPI_Comm         myComm;
    OCP_INT          numproc, myrank;

public:
    /// Print level
    USI         printLevel{0};
    /// Time control 
    ControlTime time;
    /// NR control    
    ControlNR   NR;
    /// Stop simulation
    OCP_BOOL    StopSim{ OCP_FALSE };
    /// Stop time
    OCP_DBL     MaxSimTime{ 1E20 };


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