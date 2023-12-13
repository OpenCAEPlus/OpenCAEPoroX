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
        time.SetFastControl(ctrlFast);
    }
    printLevel = ctrlFast.printLevel;
}



void OCPControl::InputParam(const ParamControl& CtrlParam)
{
    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Input Control Params -- begin");
    }

    workDir = CtrlParam.workDir;
    ocpFile = CtrlParam.fileName;    
    MaxSimTime = CtrlParam.MaxSimTime;

    SM.SetCtrlParam(CtrlParam);

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


void OCPControl::OutputModelMethodInfo() const
{
    SM.OutputInfo(printLevel);
}


void OCPControl::ApplyControl(const USI& i, const OCP_BOOL& wellOptChange_loc)
{
    /// Apply ith tuning for ith TSTEP
    time.SetNextTSTEP(i, wellOptChange_loc);
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