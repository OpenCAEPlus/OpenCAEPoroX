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


void OCPControl::SetupFastControl(const USI& argc, const char* optset[])
{
    FastControl ctrlFast(argc, optset);
    if (ctrlFast.ifUse) {
        method.push_back(ctrlFast.method);
        switch (method[0]) {
        case OCPNLMethod::IMPEC:
            lsFile.push_back("./csr.fasp");
            break;
        case OCPNLMethod::AIMc:
        case OCPNLMethod::FIM:
            lsFile.push_back("./bsr.fasp");
            break;
        default:
            OCP_ABORT("Wrong method specified from command line!");
            break;
        }
        time.SetFastControl(ctrlFast);
    }
    printLevel = ctrlFast.printLevel;
}


void OCPControl::InputParam(const ParamControl& CtrlParam)
{
    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Input Control Params -- begin");
    }


    model   = CtrlParam.model;
    workDir = CtrlParam.workDir;
    ocpFile = CtrlParam.fileName;

    for (const auto& m : CtrlParam.method) {
        if (m == "IMPEC") {
            method.push_back(OCPNLMethod::IMPEC);
        }
        else if (m == "FIM") {
            method.push_back(OCPNLMethod::FIM);
        }
        else if (m == "AIMc") {
            method.push_back(OCPNLMethod::AIMc);
        }
        else if (m == "FIMddm") {
            method.push_back(OCPNLMethod::FIMddm);
        }
        else {
            OCP_ABORT("Wrong method specified!");
        }
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