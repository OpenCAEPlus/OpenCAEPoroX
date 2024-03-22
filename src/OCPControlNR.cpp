/*! \file    OCPControl.cpp
 *  \brief   OCPControl class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPControlNR.hpp"


ControlNRParam::ControlNRParam(const vector<OCP_DBL>& src)
{
    maxIter = src[0];
    tol = src[1];
    dPmax = src[2];
    dSmax = src[3];
    dPmin = src[4];
    dSmin = src[5];
    Verrmax = src[6];
}


OCPNRStateC ControlNR::CheckConverge(const OCPNRsuite& NRs, const initializer_list<string>& il) const
{
    OCPNRStateC conflag_loc = OCPNRStateC::not_converge;
    for (auto& s : il) {
        if (s == "res") {
            if ((NRs.res.maxRelRes_V <= NRs.res.maxRelRes0_V * wp->tol ||
                NRs.res.maxRelRes_V <= wp->tol ||
                NRs.res.maxRelRes_N <= wp->tol) &&
                NRs.res.maxWellRelRes_mol <= wp->tol) {
                conflag_loc = OCPNRStateC::converge;
            }
        }
        else if (s == "resT") {
            if (((NRs.res.maxRelRes_V <= NRs.res.maxRelRes0_V * wp->tol ||
                NRs.res.maxRelRes_V <= wp->tol ||
                NRs.res.maxRelRes_N <= wp->tol) &&
                NRs.res.maxWellRelRes_mol <= wp->tol)) {
                conflag_loc = OCPNRStateC::converge;
            }
        }
        else if (s == "d") {
            if (fabs(NRs.DPBmaxNRc()) <= wp->dPmin && fabs(NRs.DSmaxNRc()) <= wp->dSmin) {
                conflag_loc = OCPNRStateC::converge;
            }
        }
        else if (s == "dT") {
            if (fabs(NRs.DPBmaxNRc()) <= wp->dPmin && fabs(NRs.DSmaxNRc()) <= wp->dSmin) {
                conflag_loc = OCPNRStateC::converge;
            }
        }
        else {
            OCP_ABORT("Iterm not recognized!");
        }
    }

    GetWallTime timer;
    timer.Start();

    OCPNRStateC conflag;
    MPI_Allreduce(&conflag_loc, &conflag, 1, OCPMPI_ENUM, MPI_MAX, myComm);

    OCPTIME_COMM_COLLECTIVE += timer.Stop();

    if (conflag == OCPNRStateC::converge) {
        // converge
        return OCPNRStateC::converge;
    }
    else if (NRs.GetIterNR() >= wp->maxIter) {
        // not converge in specified numbers of iterations, reset
        if (CURRENT_RANK == MASTER_PROCESS) {
            cout << "### WARNING: NR not fully converged!\n";
        }
        return OCPNRStateC::not_converge;
    }
    else {
        // not converge and go on iterativing
        return OCPNRStateC::continueIter;
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/05/2023      Create file                          */
/*----------------------------------------------------------------------------*/