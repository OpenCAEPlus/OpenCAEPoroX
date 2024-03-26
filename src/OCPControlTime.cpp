/*! \file    OCPControlTime.cpp
 *  \brief   OCPControlTime class definition
 *  \author  Shizhe Li
 *  \date    Dec/05/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPControlTime.hpp"


ControlTimeParam::ControlTimeParam(const TuningPair& src, const vector<OCP_DBL>& Tstep, const USI& i)
{
    const auto& src_t = src.Tuning[0];
    timeInit = src_t[0];
    timeMax = src_t[1];
    timeMin = src_t[2];
    maxIncreFac = src_t[3];
    minChopFac = src_t[4];
    cutFacNR = src_t[5];

    const auto& src_pt = src.Tuning[1];
    dPlim = src_pt[0];
    dSlim = src_pt[1];
    dNlim = src_pt[2];
    eVlim = src_pt[3];

    begin_time = Tstep[i];
    end_time = Tstep[i + 1];
}


void ControlTimeParam::SetFastControl(const FastControl& fCtrl)
{
    timeInit = fCtrl.timeInit;
    timeMax = fCtrl.timeMax;
    timeMin = fCtrl.timeMin;
}


void ControlTime::CutDt(const OCP_DBL& fac)
{
    const OCP_DBL ldt = current_dt;

    if (fac < 0) current_dt *= wp->cutFacNR;
    else         current_dt *= fac;

    if (CURRENT_RANK == MASTER_PROCESS) {
        cout << "### WARNING: Cut time step size: " << scientific
            << setprecision(3) << ldt << TIMEUNIT + " -> "
            << setprecision(3) << current_dt << TIMEUNIT + "!\n";
    }
}


void ControlTime::CutDt(const OCPNRsuite& NRs)
{
    const OCP_DBL ldt = current_dt;
    switch (NRs.GetWorkState())
    {
    case OCPNRStateP::reset:
        // do not cut
        return;

    case OCPNRStateP::resetCut:
        current_dt *= wp->cutFacNR;
        break;

    case OCPNRStateP::resetCutCFL:
        current_dt /= (NRs.GetMaxCFL() + 1);
        break;

    default:
        OCP_ABORT("WRONG work state!");
    }

    if (CURRENT_RANK == MASTER_PROCESS) {
        cout << "### WARNING: Cut time step size: " << scientific
            << setprecision(3) << ldt << TIMEUNIT + " -> "
            << setprecision(3) << current_dt << TIMEUNIT + "!\n";
    }
}


INT ControlTime::SetNextTSTEP() 
{ 
    USI d = 0;
    for (d = 0; d < ps.size(); d++) {
        if (current_time > ps[d].begin_time - TINY && current_time < ps[d].end_time) {
            wp = &ps[d];
            return d;
        }
    }

    return -1;
}


void ControlTime::CalInitTimeStep4TSTEP(const OCP_BOOL& wellOptChange_loc)
{
    /// Set initial time step for next TSTEP
    GetWallTime timer;
    timer.Start();

    OCP_BOOL       wellOptChange;
    MPI_Allreduce(&wellOptChange_loc, &wellOptChange, 1, OCPMPI_BOOL, MPI_LAND, myComm);

    timer.Stop();


    const OCP_DBL dt = wp->end_time - current_time;
    if (dt <= 0) OCP_ABORT("Non-positive time stepsize!");

    static OCP_BOOL firstflag = OCP_TRUE;
    if (wellOptChange || firstflag) {
        current_dt = min(dt, wp->timeInit);
        firstflag = OCP_FALSE;
    }
    else {
        current_dt = min(dt, predict_dt);
    }
}


void ControlTime::CalNextTimeStep(const OCPNRsuite& NRs, const initializer_list<string>& il)
{
    last_dt = current_dt;
    current_time += current_dt;

    OCP_DBL factor = wp->maxIncreFac;

    const OCP_DBL dPmax = fabs(NRs.DPmaxT());
    const OCP_DBL dTmax = fabs(NRs.DTmaxT());
    const OCP_DBL dNmax = fabs(NRs.DNmaxT());
    const OCP_DBL dSmax = fabs(NRs.DSmaxT());
    const OCP_DBL eVmax = fabs(NRs.EVmaxT());

    for (auto& s : il) {
        if (s == "dP") {
            if (dPmax > TINY) factor = min(factor, wp->dPlim / dPmax);
        }
        else if (s == "dT") {
            // no input now -- no value
            if (dTmax > TINY) factor = min(factor, wp->dTlim / dTmax);
        }
        else if (s == "dN") {
            if (dNmax > TINY) factor = min(factor, wp->dNlim / dNmax);
        }
        else if (s == "dS") {
            if (dSmax > TINY) factor = min(factor, wp->dSlim / dSmax);
        }
        else if (s == "eV") {
            if (eVmax > TINY) factor = min(factor, wp->eVlim / eVmax);
        }
        else if (s == "iter") {
            if (NRs.GetIterNR() < 5)
                factor = min(factor, static_cast<OCP_DBL>(2.0));
            else if (NRs.GetIterNR() > 10)
                factor = min(factor, static_cast<OCP_DBL>(0.5));
            else
                factor = min(factor, static_cast<OCP_DBL>(1.5));
        }
        else {
            OCP_ABORT("Iterm not recognized!");
        }
    }

    factor = max(wp->minChopFac, factor);

    OCP_DBL dt_loc = current_dt * factor;
    if (dt_loc > wp->timeMax) dt_loc = wp->timeMax;
    if (dt_loc < wp->timeMin) dt_loc = wp->timeMin;

    GetWallTime timer;
    timer.Start();

    MPI_Allreduce(&dt_loc, &current_dt, 1, OCPMPI_DBL, MPI_MIN, myComm);

    OCPTIME_COMM_COLLECTIVE += timer.Stop();

    predict_dt = current_dt;

    if (current_dt > (wp->end_time - current_time))
        current_dt = (wp->end_time - current_time);
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/05/2023      Create file                          */
/*----------------------------------------------------------------------------*/