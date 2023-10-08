/*! \file    OCPNRsuite.cpp
 *  \brief   data structure used in NR iterations
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


void OCPRes::SetInitRes() 
{ 
    maxRelRes0_V = maxRelRes_V; 
}


void OCPNRsuite::SetupIsoT(const OCP_USI& nbin, const USI& npin, const USI& ncin) 
{
    nb = nbin;
    np = npin;
    nc = ncin;

    lP.resize(nb);
    lN.resize(nb * nc);
    lS.resize(nb * np);
    dPNR.resize(nb);
    dNNR.resize(nb * nc);
    dSNR.resize(nb * np);
}


void OCPNRsuite::SetupT(const OCP_USI& nbin, const USI& npin, const USI& ncin)
{
    nb = nbin;
    np = npin;
    nc = ncin;

    lP.resize(nb);
    lT.resize(nb);
    lN.resize(nb * nc);
    lS.resize(nb * np);
    dPNR.resize(nb);
    dTNR.resize(nb);
    dNNR.resize(nb * nc);
    dSNR.resize(nb * np);
}


void OCPNRsuite::Reset(const BulkVarSet& bvs)
{
    lP = bvs.lP;
    lT = bvs.lT;
    lN = bvs.lNi;
    lS = bvs.lS;
}


void OCPNRsuite::CaldMaxIsoT(const BulkVarSet& bvs)
{
    dPmax = 0;
    dNmax = 0;
    dSmax = 0;

    for (OCP_USI n = 0; n < nb; n++) {
        dPNR[n] = bvs.P[n] - lP[n];
        if (dPmax < fabs(dPNR[n]))  dPmax = dPNR[n];

        OCP_DBL  Nt = 0;
        for (USI i = 0; i < nc; i++) {
            Nt += lN[n * nc + i];
            dNNR[n * nc + i] = bvs.Ni[n * nc + i] - lN[n * nc + i];
        }
        for (USI i = 0; i < nc; i++) {
            if (fabs(dNmax) < fabs(dNNR[n * nc + i] / Nt))  dNmax = dNNR[n * nc + i] / Nt;
        }

        for (USI j = 0; j < np; j++) {
            if (fabs(dSmax) < fabs(dSNR[n]))  dSmax = dSNR[n];
        }
    }

    lP = bvs.P;
    lN = bvs.Ni;
    lS = bvs.S;
}


void OCPNRsuite::CaldMaxT(const BulkVarSet& bvs)
{
    dPmax = 0;
    dTmax = 0;
    dNmax = 0;
    dSmax = 0;

    for (OCP_USI n = 0; n < nb; n++) {
        dPNR[n] = bvs.P[n] - lP[n];
        if (fabs(dPmax) < fabs(dPNR[n]))  dPmax = dPNR[n];
        dTNR[n] = bvs.T[n] - lT[n];
        if (fabs(dTmax) < fabs(dTNR[n]))  dTmax = dTNR[n];

        OCP_DBL  Nt   = 0;
        for (USI i = 0; i < nc; i++) {
            Nt               += lN[n * nc + i];
            dNNR[n * nc + i] = bvs.Ni[n * nc + i] - lN[n * nc + i];
        }
        for (USI i = 0; i < nc; i++) {
            if (fabs(dNmax) < fabs(dNNR[n * nc + i] / Nt))  dNmax = dNNR[n * nc + i] / Nt;
        }

        for (USI j = 0; j < np; j++) {
            if (fabs(dSmax) < fabs(dSNR[n]))  dSmax = dSNR[n];
        }
    }

    lP = bvs.P;
    lT = bvs.T;
    lN = bvs.Ni;
    lS = bvs.S;
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/07/2023      Create file                          */
/*----------------------------------------------------------------------------*/