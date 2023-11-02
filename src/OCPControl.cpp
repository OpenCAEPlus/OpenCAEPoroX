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

#include "OCPControl.hpp"


FastControl::FastControl(const USI& argc, const char* optset[])
{
    ifUse = OCP_FALSE;
    timeInit = timeMax = timeMin = -1.0;

    std::stringstream buffer;
    string            tmp;
    string            key;
    string            value;
    for (USI n = 2; n < argc; n++) {
        buffer << optset[n];
        buffer >> tmp;

        string::size_type pos = tmp.find_last_of('=');
        if (pos == string::npos) OCP_ABORT("Unknown Usage! See -h");

        key = tmp.substr(0, pos);
        value = tmp.substr(pos + 1, tmp.size() - pos);

        switch (Map_Str2Int(&key[0], key.size())) {

        case Map_Str2Int("method", 6):
            if (value == "FIM") {
                method = FIM;
            }
            else if (value == "IMPEC") {
                method = IMPEC;
            }
            else if (value == "AIMc") {
                method = AIMc;
            }
            else {
                OCP_ABORT("Wrong method param in command line!");
            }
            ifUse = OCP_TRUE;
            if (method == FIM || method == AIMc) {
                if (timeInit <= 0) timeInit = 1;
                if (timeMax <= 0) timeMax = 10.0;
                if (timeMin <= 0) timeMin = 0.1;
            }
            else {
                if (timeInit <= 0) timeInit = 0.1;
                if (timeMax <= 0) timeMax = 1.0;
                if (timeMin <= 0) timeMin = 0.1;
            }
            break;

        case Map_Str2Int("dtInit", 6):
            timeInit = stod(value);
            break;

        case Map_Str2Int("dtMin", 5):
            timeMin = stod(value);
            break;

        case Map_Str2Int("dtMax", 5):
            timeMax = stod(value);
            break;

        case Map_Str2Int("verbose", 7):
            printLevel = OCP_MIN(OCP_MAX(stoi(value), PRINT_NONE), PRINT_ALL);
            break;

        default:
            OCP_ABORT("Unknown Options: " + key + "   See -h");
            break;
        }

        buffer.clear();
    }
}


void OCPControl::SetupFastControl(const USI& argc, const char* optset[])
{
    FastControl ctrlFast(argc, optset);
    if (ctrlFast.ifUse) {
        method = ctrlFast.method;
        switch (method) {
        case IMPEC:
            lsFile = "./csr.fasp";
            break;
        case AIMc:
        case FIM:
            lsFile = "./bsr.fasp";
            break;
        default:
            OCP_ABORT("Wrong method specified from command line!");
            break;
        }
        time.SetFastControl(ctrlFast);
    }
    printLevel = ctrlFast.printLevel;
}


ControlTimeParam::ControlTimeParam(const TuningPair& src, const vector<OCP_DBL>& Tstep, const USI& i)
{
    const auto& src_t = src.Tuning[0];
    timeInit    = src_t[0];
    timeMax     = src_t[1];
    timeMin     = src_t[2];
    maxIncreFac = src_t[3];
    minChopFac  = src_t[4];
    cutFacNR    = src_t[5];

    const auto& src_pt = src.Tuning[1];
    dPlim       = src_pt[0];
    dSlim       = src_pt[1];
    dNlim       = src_pt[2];
    eVlim       = src_pt[3];

    begin_time = Tstep[i];
    end_time   = Tstep[i + 1];
}


void ControlTimeParam::SetFastControl(const FastControl& fCtrl)
{
    timeInit = fCtrl.timeInit;
    timeMax  = fCtrl.timeMax;
    timeMin  = fCtrl.timeMin;
}


void ControlTime::CutDt(const OCP_DBL& fac) 
{
    const OCP_DBL ldt = current_dt;

    if (fac < 0) current_dt *= wp->cutFacNR;
    else         current_dt *= fac;

    if (CURRENT_RANK == MASTER_PROCESS) {
        cout << "### WARNING: Cut time step size: " << fixed
            << setprecision(3) << ldt << TIMEUNIT +  " -> "
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
        cout << "### WARNING: Cut time step size: " << fixed
            << setprecision(3) << ldt << TIMEUNIT + " -> "
            << setprecision(3) << current_dt << TIMEUNIT + "!\n";
    }
}


void ControlTime::SetNextTSTEP(const USI& i, const AllWells& wells)
{
    wp = &ps[i];

    /// Set initial time step for next TSTEP
    GetWallTime timer;
    timer.Start();

    OCP_BOOL       wellOptChange;
    const OCP_BOOL wellChange_loc = wells.GetWellOptChange();
    MPI_Allreduce(&wellChange_loc, &wellOptChange, 1, OCPMPI_BOOL, MPI_LAND, myComm);

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

	OCPTIME_COMM_COLLECTIVE += timer.Stop() / TIME_S2MS;

	predict_dt = current_dt;

	if (current_dt > (wp->end_time - current_time))
		current_dt = (wp->end_time - current_time);
}


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

    OCPTIME_COMM_COLLECTIVE += timer.Stop() / TIME_S2MS;

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


void OCPControl::InputParam(const ParamControl& CtrlParam)
{
    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Input Control Params -- begin");
    }


    model   = CtrlParam.model;
    workDir = CtrlParam.workDir;
    ocpFile = CtrlParam.fileName;

    if (CtrlParam.method == "IMPEC") {
        method = IMPEC;
    } else if (CtrlParam.method == "FIM") {
        method = FIM;
    } else if (CtrlParam.method == "AIMc") {
        method = AIMc;
    } else {
        OCP_ABORT("Wrong method specified!");
    }

    lsFile = CtrlParam.lsFile;
    
    MaxSimTime = CtrlParam.MaxSimTime;

    const USI   n = CtrlParam.tuning_T.size();
    vector<USI> ctrlCriticalTime(n + 1);
    for (USI i = 0; i < n; i++) {
        ctrlCriticalTime[i] = CtrlParam.tuning_T[i].d;
    }
    ctrlCriticalTime.back() = CtrlParam.criticalTime.size() - 1;
    for (USI i = 0; i < n; i++) {
        for (USI d = ctrlCriticalTime[i]; d < ctrlCriticalTime[i + 1]; d++) {
            time.SetCtrlParam(CtrlParam.tuning_T[i], CtrlParam.criticalTime, d);
            NR.SetCtrlParam(CtrlParam.tuning_T[i].Tuning[2]);
        }
    }


    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Input Control Params -- end");
    }
}


void OCPControl::Setup(const Domain& domain) 
{
    myComm  = domain.myComm;
    numproc = domain.numproc;
    myrank  = domain.myrank;

    time.SetupComm(domain);
    NR.SetupComm(domain);
}


void OCPControl::ApplyControl(const USI& i, const Reservoir& rs)
{
    /// Apply ith tuning for ith TSTEP
    time.SetNextTSTEP(i, rs.allWells);
    NR.SetNextTSTEP(i);
}




/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/