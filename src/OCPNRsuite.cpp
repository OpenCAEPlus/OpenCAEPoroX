/*! \file    OCPNLsuite.cpp
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


void OCPRes::SetupIsoT(const OCP_USI& nb, const OCP_USI& nw, const OCP_USI& nc)
{
    OCP_USI reslen = (nb + nw) * (nc + 1);
    resAbs.resize(reslen);
    resRelV.resize(nb);
    resRelN.resize(nb);
}


void OCPRes::SetupT(const OCP_USI& nb, const OCP_USI& nw, const OCP_USI& nc)
{
    OCP_USI reslen = (nb + nw) * (nc + 2);
    resAbs.resize(reslen);
    resRelV.resize(nb);
    resRelN.resize(nb);
    resRelE.resize(nb);
}


void OCPRes::SetZero()
{
    fill(resAbs.begin(), resAbs.end(), 0);
    fill(resRelV.begin(), resRelV.end(), 0);
    fill(resRelN.begin(), resRelN.end(), 0);
    fill(resRelE.begin(), resRelE.end(), 0);
    maxRelRes_V       = 0;
    maxRelRes_N       = 0;
    maxRelRes_E       = 0;
    maxWellRelRes_mol = 0;
    maxId_V           = 0;
    maxId_N           = 0;
    maxId_E           = 0;
}


void OCPNLsuite::Setup(const OCP_BOOL& ifthermal, const BulkVarSet& bvs, const OCP_USI& nw, const Domain& domain)
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


void OCPNLsuite::InitStep(const BulkVarSet& bvs)
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

    dPmaxNR.clear();
    dTmaxNR.clear();
    dNmaxNR.clear();
    dSmaxNR.clear();
}


void OCPNLsuite::CaldMax(const BulkVarSet& bvs)
{
    OCP_DBL dPmaxTmp = 0;
    OCP_DBL dTmaxTmp = 0;
    OCP_DBL dNmaxTmp = 0;
    OCP_DBL dSmaxTmp = 0;

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

    dPmaxNR.push_back(dPmaxTmp);
    dTmaxNR.push_back(dTmaxTmp);
    dNmaxNR.push_back(dNmaxTmp);
    dSmaxNR.push_back(dSmaxTmp);
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/07/2023      Create file                          */
/*----------------------------------------------------------------------------*/