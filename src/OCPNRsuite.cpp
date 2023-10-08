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


void OCPNRsuite::SetupIsoT(const BulkVarSet& bvs, const OCP_USI& nw)
{
    nb = bvs.nbI;
    np = bvs.np;
    nc = bvs.nc;

    res.SetupIsoT(nb, nw, nc);

    lP.resize(nb);
    lN.resize(nb * nc);
    lS.resize(nb * np);
    dP.resize(nb);
    dN.resize(nb * nc);
    dS.resize(nb * np);

}


void OCPNRsuite::SetupT(const BulkVarSet& bvs, const OCP_USI& nw)
{
    nb = bvs.nbI;
    np = bvs.np;
    nc = bvs.nc;

    res.SetupT(nb, nw, nc);

    lP.resize(nb);
    lT.resize(nb);
    lN.resize(nb * nc);
    lS.resize(nb * np);
    dP.resize(nb);
    dT.resize(nb);
    dN.resize(nb * nc);
    dS.resize(nb * np);
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
        dP[n] = bvs.P[n] - lP[n];
        if (dPmax < fabs(dP[n]))  dPmax = dP[n];

        OCP_DBL  Nt = 0;
        for (USI i = 0; i < nc; i++) {
            Nt += lN[n * nc + i];
            dN[n * nc + i] = bvs.Ni[n * nc + i] - lN[n * nc + i];
        }
        for (USI i = 0; i < nc; i++) {
            if (fabs(dNmax) < fabs(dN[n * nc + i] / Nt))  dNmax = dN[n * nc + i] / Nt;
        }

        for (USI j = 0; j < np; j++) {
            dS[n * np + j] = bvs.S[n * np + j] - lS[n * np + j];
            if (fabs(dSmax) < fabs(dS[n * np + j]))  dSmax = dS[n * np + j];
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
        dP[n] = bvs.P[n] - lP[n];
        if (fabs(dPmax) < fabs(dP[n]))  dPmax = dP[n];
        dT[n] = bvs.T[n] - lT[n];
        if (fabs(dTmax) < fabs(dT[n]))  dTmax = dT[n];

        OCP_DBL  Nt   = 0;
        for (USI i = 0; i < nc; i++) {
            Nt               += lN[n * nc + i];
            dN[n * nc + i] = bvs.Ni[n * nc + i] - lN[n * nc + i];
        }
        for (USI i = 0; i < nc; i++) {
            if (fabs(dNmax) < fabs(dN[n * nc + i] / Nt))  dNmax = dN[n * nc + i] / Nt;
        }

        for (USI j = 0; j < np; j++) {
            dS[n * np + j] = bvs.S[n * np + j] - lS[n * np + j];
            if (fabs(dSmax) < fabs(dS[n * np + j]))  dSmax = dS[n * np + j];
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