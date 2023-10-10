/*! \file    OCPNRsuite.cpp
 *  \brief   data structure used in non-linear iterations
 *  \author  Shizhe Li
 *  \date    Oct/07/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


 // Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "OCPNRsuite.hpp"


void OCPNRsuite::Setup(const OCP_BOOL& ifthermal, const BulkVarSet& bvs, const OCP_USI& nw, const Domain& domain)
{
    myComm  = domain.myComm;
    numproc = domain.numproc;
    myrank  = domain.myrank;

    ifThermal = ifthermal;
    nb        = bvs.nbI;
    np        = bvs.np;
    nc        = bvs.nc;

    if (ifThermal) {
        res.SetupT(nb, nw, nc);
    }
    else {
        res.SetupIsoT(nb, nw, nc);
    }


    lP.resize(nb);
    lT.resize(nb);
    lN.resize(nb * nc);
    lS.resize(nb * np);
    dP.resize(nb);
    dT.resize(nb);
    dN.resize(nb * nc);
    dS.resize(nb * np);
}


void OCPNRsuite::InitStep(const BulkVarSet& bvs)
{
    GetWallTime timer;
    timer.Start();
    
    OCP_DBL tmploc = res.maxRelRes_V;
    MPI_Allreduce(&tmploc, &res.maxRelRes0_V, 1, MPI_DOUBLE, MPI_MIN, myComm);
    
    OCPTIME_COMM_COLLECTIVE += timer.Stop() / 1000;


    lP = bvs.lP;
    lT = bvs.lT;
    lN = bvs.lNi;
    lS = bvs.lS;

    dPBmaxNR.clear();
    dTmaxNR.clear();
    dNmaxNR.clear();
    dSmaxNR.clear();
    dPWmaxNR.clear();
}


void OCPNRsuite::CalMaxChangeNR(const Reservoir& rs)
{
    OCP_DBL dPmaxTmp = 0;
    OCP_DBL dTmaxTmp = 0;
    OCP_DBL dNmaxTmp = 0;
    OCP_DBL dSmaxTmp = 0;

    // for bulk
    const BulkVarSet& bvs = rs.bulk.GetVarSet();

    for (OCP_USI n = 0; n < nb; n++) {
        dP[n] = bvs.P[n] - lP[n];
        if (fabs(dPmaxTmp) < fabs(dP[n]))  dPmaxTmp = dP[n];
        dT[n] = bvs.T[n] - lT[n];
        if (fabs(dTmaxTmp) < fabs(dT[n]))  dTmaxTmp = dT[n];

        OCP_DBL  Nt   = 0;
        for (USI i = 0; i < nc; i++) {
            Nt             += lN[n * nc + i];
            dN[n * nc + i] = bvs.Ni[n * nc + i] - lN[n * nc + i];
        }
        for (USI i = 0; i < nc; i++) {
            if (fabs(dNmaxTmp) < fabs(dN[n * nc + i] / Nt))  dNmaxTmp = dN[n * nc + i] / Nt;
        }

        for (USI j = 0; j < np; j++) {
            dS[n * np + j] = bvs.S[n * np + j] - lS[n * np + j];
            if (fabs(dSmaxTmp) < fabs(dS[n * np + j]))  dSmaxTmp = dS[n * np + j];
        }
    }

    lP = bvs.P;
    lT = bvs.T;
    lN = bvs.Ni;
    lS = bvs.S;

    dPBmaxNR.push_back(dPmaxTmp);
    dTmaxNR.push_back(dTmaxTmp);
    dNmaxNR.push_back(dNmaxTmp);
    dSmaxNR.push_back(dSmaxTmp);

    // for well   -- wrong now
    dPmaxTmp = 0;
    const AllWells& well = rs.allWells;
    for (const auto& w : well.wells) {
        const OCP_DBL dPw = w->CalMaxChangeNR();
        if (fabs(dPmaxTmp) < fabs(dPw)) {
            dPmaxTmp = dPw;
        }
    }
    dPWmaxNR.push_back(dPmaxTmp);
}


void OCPNRsuite::CalMaxChangeTime(const Reservoir& rs)
{
    OCP_FUNCNAME;

    dPmaxT  = 0;
    dPBmaxT = 0;
    dPWmaxT = 0;
    dTmaxT  = 0;
    dNmaxT  = 0;
    dSmaxT  = 0;
    eVmaxT  = 0;

    // for bulk
    const BulkVarSet& bvs = rs.bulk.GetVarSet();

    OCP_USI id;
    for (OCP_USI n = 0; n < nb; n++) {

        // dP
        if (fabs(dPBmaxT) < fabs(bvs.P[n] - bvs.lP[n])) {
            dPBmaxT = bvs.P[n] - bvs.lP[n];
        }

        // dT
        if (fabs(dTmaxT) < fabs(bvs.T[n] - bvs.lT[n])) {
            dTmaxT = bvs.T[n] - bvs.lT[n];
        }

        // dS
        for (USI j = 0; j < np; j++) {
            id = n * np + j;
            if (fabs(dSmaxT) < fabs(bvs.S[id] - bvs.lS[id])) {
                dSmaxT = bvs.S[id] - bvs.lS[id];
            }
        }

        // dN
        for (USI i = 0; i < nc; i++) {
            id = n * nc + i;
            const OCP_DBL tmp = fabs(max(bvs.Ni[id], bvs.lNi[id]));
            if (tmp > TINY) {
                if (fabs(dNmaxT) < fabs((bvs.Ni[id] - bvs.lNi[id]) / tmp)) {
                    dNmaxT = (bvs.Ni[id] - bvs.lNi[id]) / tmp;
                }
            }
        }

        // Ve
        if (fabs(eVmaxT) < fabs((bvs.vf[n] - bvs.rockVp[n]) / bvs.rockVp[n])) {
            eVmaxT = (bvs.vf[n] - bvs.rockVp[n]) / bvs.rockVp[n];
        }
    }

    // for well
    const AllWells& well = rs.allWells;
    for (const auto& w : well.wells) {
        OCP_DBL dPw = w->CalMaxChangeTime();
        if (fabs(dPWmaxT) < fabs(dPw)) {
            dPWmaxT = dPw;
        }
    }

    if (fabs(dPBmaxT) < fabs(dPWmaxT)) {
        dPmaxT = dPWmaxT;
    }
    else {
        dPmaxT = dPBmaxT;
    }
}



void OCPNRsuite::InitIter() {
    iterNR  = 0;
    iterLS  = 0;
    iterNRw = 0;
    iterLSw = 0;
    iterNRLS.clear();
}


void OCPNRsuite::UpdateIter(const USI& lsIter) 
{
    iterNR++;
    iterLS += lsIter;
    iterNRLS.push_back(lsIter);
}


void OCPNRsuite::ResetIter() 
{
    iterNRw += iterNR;
    iterLSw += iterLS;
    iterNR = 0;
    iterLS = 0;
    iterNRLS.clear();
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/07/2023      Create file                          */
/*----------------------------------------------------------------------------*/