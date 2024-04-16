/*! \file    OCPControlTime.hpp
 *  \brief   OCPControlTime class declaration
 *  \author  Shizhe Li
 *  \date    Dec/05/2023 
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPCONTROLTIME_HEADER__
#define __OCPCONTROLTIME_HEADER__

 // Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "ParamControl.hpp"
#include "OCPControlFast.hpp"
#include "OCPNRsuite.hpp"

using namespace std;


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
        myComm = domain.global_comm;
        numproc = domain.global_numproc;
        myrank = domain.global_rank;
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
    INT SetNextTSTEP();
    /// Calculate initial time step for next TSTEP
    void CalInitTimeStep4TSTEP(const OCP_BOOL& wellOptChange_loc);
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
    /// Set current time
    void SetStartTime(const OCP_DBL& startT) { current_time = startT; }
    /// Return the current time.
    auto GetCurrentTime() const { return current_time; }
    /// Return current time step size.
    auto GetCurrentDt() const { return current_dt; }
    /// Return last time step size.
    auto GetLastDt() const { return last_dt; }
    /// Determine whether the critical time point has been reached.
    auto IfEndTSTEP() { return ((wp->end_time - current_time) < TINY); }

protected:
    /// control param set
    vector<ControlTimeParam> ps;
    /// work param
    const ControlTimeParam* wp;
    /// from prediction for next TSTEP
    OCP_DBL                  predict_dt;
    /// current time step
    OCP_DBL                  current_dt;
    /// last time step
    OCP_DBL                  last_dt{ 0 };
    /// current time
    OCP_DBL                  current_time{ 0 };
};



#endif /* end if __OCPControlTime_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/05/2023      Create file                          */
/*----------------------------------------------------------------------------*/