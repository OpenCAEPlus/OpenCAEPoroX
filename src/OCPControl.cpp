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

ControlTime::ControlTime(const vector<OCP_DBL>& src_t, const vector<OCP_DBL>& src_pt)
{
    timeInit    = src_t[0];
    timeMax     = src_t[1];
    timeMin     = src_t[2];
    maxIncreFac = src_t[3];
    minChopFac  = src_t[4];
    cutFacNR    = src_t[5];

    dPlim       = src_pt[0];
    dSlim       = src_pt[1];
    dNlim       = src_pt[2];
    eVlim       = src_pt[3];
}


ControlNR::ControlNR(const vector<OCP_DBL>& src)
{
    maxIter = src[0];
    NRtol     = src[1];
    dPmax   = src[2];
    dSmax   = src[3];
    dPmin   = src[4];
    dSmin   = src[5];
    Verrmax   = src[6];
}

void FastControl::ReadParam(const USI& argc, const char* optset[])
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

        key   = tmp.substr(0, pos);
        value = tmp.substr(pos + 1, tmp.size() - pos);

        switch (Map_Str2Int(&key[0], key.size())) {

            case Map_Str2Int("method", 6):
                if (value == "FIM") {
                    method = FIM;
                } else if (value == "IMPEC") {
                    method = IMPEC;
                } else if (value == "AIMc") {
                    method = AIMc;
                } else {
                    OCP_ABORT("Wrong method param in command line!");
                }
                ifUse = OCP_TRUE;
                if (method == FIM || method == AIMc) {
                    if (timeInit <= 0) timeInit = 1;
                    if (timeMax <= 0) timeMax = 10.0;
                    if (timeMin <= 0) timeMin = 0.1;
                } else {
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


void ItersInfo::UpdateTotal()
{
    numTstep += 1;
    NRt      += NR;
    LSt      += LS;
    NR        = 0;
    LS        = 0;
}


void ItersInfo::Reset()
{
    NRwt += NR;  
    LSwt += LS;
    NR    = 0;
    LS    = 0;
}


void OCPControl::InputParam(const ParamControl& CtrlParam)
{
    model    = CtrlParam.model;
    workDir  = CtrlParam.workDir;
    fileName = CtrlParam.fileName;

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
    criticalTime     = CtrlParam.criticalTime;

    USI t = CtrlParam.criticalTime.size();
    ctrlTimeSet.resize(t);
    ctrlNRSet.resize(t);

    USI         n = CtrlParam.tuning_T.size();
    vector<USI> ctrlCriticalTime(n + 1);
    for (USI i = 0; i < n; i++) {
        ctrlCriticalTime[i] = CtrlParam.tuning_T[i].d;
    }
    ctrlCriticalTime.back() = t;
    for (USI i = 0; i < n; i++) {
        for (USI d = ctrlCriticalTime[i]; d < ctrlCriticalTime[i + 1]; d++) {
            ctrlTimeSet[d] = ControlTime(CtrlParam.tuning_T[i].Tuning[0], CtrlParam.tuning_T[i].Tuning[1]);
            ctrlNRSet[d]   = ControlNR(CtrlParam.tuning_T[i].Tuning[2]);
        }
    }
}


void OCPControl::Setup(const Domain& domain) 
{
    myComm  = domain.myComm;
    numproc = domain.numproc;
    myrank  = domain.myrank;
}


void OCPControl::ApplyControl(const USI& i, const Reservoir& rs)
{
    /// Apply ith tuning for ith TSTEP
    ctrlTime    = ctrlTimeSet[i];
    ctrlNR      = ctrlNRSet[i];
    end_time    = criticalTime[i + 1];

    /// Set initial time step for next TSTEP
    GetWallTime timer;
    timer.Start();

    OCP_BOOL       wellOptChange;
    const OCP_BOOL wellChange_loc  = rs.allWells.GetWellOptChange();
    MPI_Allreduce(&wellChange_loc, &wellOptChange, 1, MPI_INT, MPI_LAND, rs.domain.myComm);

    timer.Stop();

    OCP_DBL dt = criticalTime[i + 1] - current_time;
    if (dt <= 0) OCP_ABORT("Non-positive time stepsize!");

    static OCP_BOOL firstflag = OCP_TRUE;
    if (wellOptChange || firstflag) {
        current_dt = min(dt, ctrlTime.timeInit);
        firstflag = OCP_FALSE;
    }
    else {
        current_dt = min(dt, predict_dt);
    }
}


void OCPControl::SetupFastControl(const USI& argc, const char* optset[])
{
    ctrlFast.ReadParam(argc, optset);
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
        USI n = ctrlTimeSet.size();
        for (USI i = 0; i < n; i++) {
            ctrlTimeSet[i].timeInit = ctrlFast.timeInit;
            ctrlTimeSet[i].timeMax  = ctrlFast.timeMax;
            ctrlTimeSet[i].timeMin  = ctrlFast.timeMin;
        }
    }
    printLevel = ctrlFast.printLevel;
}


OCP_BOOL OCPControl::Check(Reservoir& rs, initializer_list<string> il)
{
    workState_loc = OCP_CONTINUE;
    OCP_INT flag;
    for (auto& s : il) {

        if (s == "BulkP")
            flag = rs.bulk.CheckP();
        else if (s == "BulkT")
            flag = rs.bulk.CheckT();
        else if (s == "BulkNi")
            flag = rs.bulk.CheckNi();
        else if (s == "BulkVe")
            flag = rs.bulk.CheckVe(0.01);
        else if (s == "CFL")
            flag = rs.bulk.CheckCFL(1.0);
        else if (s == "WellP")
            flag = rs.allWells.CheckP(rs.bulk);
        else
            OCP_ABORT("Check iterm not recognized!");

        switch (flag) {
            // Bulk
            case BULK_SUCCESS:
                break;

            case BULK_NEGATIVE_PRESSURE:
            case BULK_NEGATIVE_TEMPERATURE:
            case BULK_NEGATIVE_COMPONENTS_MOLES:
            case BULK_OUTRANGED_VOLUME_ERROR:
                workState_loc = OCP_RESET_CUTTIME;
                break;

            case BULK_OUTRANGED_CFL:
                workState_loc = OCP_RESET_CUTTIME_CFL;
                break;

            // Well
            case WELL_SUCCESS:
                break;

            case WELL_NEGATIVE_PRESSURE:
                workState_loc = OCP_RESET_CUTTIME;
                break;

            case WELL_SWITCH_TO_BHPMODE:
            case WELL_CROSSFLOW:
                workState_loc = OCP_RESET;
                break;

            default:
                break;
        }

        if (workState_loc != OCP_CONTINUE)
            break;
    }

    GetWallTime timer;
    timer.Start();

    MPI_Allreduce(&workState_loc, &workState, 1, MPI_INT, MPI_MIN, myComm);

    OCPTIME_COMM_COLLECTIVE += timer.Stop() / 1000;
    OCPTIME_COMM_1ALLREDUCE += timer.Stop() / 1000;

    switch (workState)
    {
    case OCP_CONTINUE:
        return OCP_TRUE;

    case OCP_RESET:
        return OCP_FALSE;

    case OCP_RESET_CUTTIME:
        current_dt *= ctrlTime.cutFacNR;
        return OCP_FALSE;

    case OCP_RESET_CUTTIME_CFL:
        current_dt /= (rs.bulk.GetMaxCFL() + 1);
        return OCP_FALSE;

    default:
        OCP_ABORT("WRONG work state!");
    }


    return OCP_TRUE;
}

void OCPControl::CalNextTimeStep(Reservoir& rs, initializer_list<string> il)
{
    last_dt = current_dt;
    current_time += current_dt;

    OCP_DBL factor = ctrlTime.maxIncreFac;

    const OCP_DBL dPmax = max(rs.bulk.GetdPmax(), rs.allWells.GetdBHPmax());
    const OCP_DBL dTmax = rs.bulk.GetdTmax();
    const OCP_DBL dNmax = rs.bulk.GetdNmax();
    const OCP_DBL dSmax = rs.bulk.GetdSmax();
    const OCP_DBL eVmax = rs.bulk.GeteVmax();

    for (auto& s : il) {
        if (s == "dP") {
            if (dPmax > TINY) factor = min(factor, ctrlTime.dPlim / dPmax);
        } else if (s == "dT") {
            // no input now -- no value
            if (dTmax > TINY) factor = min(factor, ctrlTime.dTlim / dTmax);
        } else if (s == "dN") {
            if (dNmax > TINY) factor = min(factor, ctrlTime.dNlim / dNmax);
        } else if (s == "dS") {
            if (dSmax > TINY) factor = min(factor, ctrlTime.dSlim / dSmax);
        } else if (s == "eV") {
            if (eVmax > TINY) factor = min(factor, ctrlTime.eVlim / eVmax);
        } else if (s == "iter") {
            if (iters.NR < 5)
                factor = min(factor, 2.0);
            else if (iters.NR > 10)
                factor = min(factor, 0.5);
            else
                factor = min(factor, 1.5);
        }
    }

    factor = max(ctrlTime.minChopFac, factor);

    OCP_DBL dt_loc = current_dt * factor;
    if (dt_loc > ctrlTime.timeMax) dt_loc = ctrlTime.timeMax;
    if (dt_loc < ctrlTime.timeMin) dt_loc = ctrlTime.timeMin;

    GetWallTime timer;
    timer.Start();

    MPI_Allreduce(&dt_loc, &current_dt, 1, MPI_DOUBLE, MPI_MIN, myComm);

    OCPTIME_COMM_COLLECTIVE += timer.Stop() / 1000;

    predict_dt = current_dt;

    if (current_dt > (end_time - current_time)) current_dt = (end_time - current_time);
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/