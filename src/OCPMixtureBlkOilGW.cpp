/*! \file    OCPMixtureBlkOilGW.cpp
 *  \brief   OCPMixtureBlkOilGW class declaration
 *  \author  Shizhe Li
 *  \date    Sep/22/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPMixtureBlkOilGW.hpp"


 /////////////////////////////////////////////////////
 // OCPMixtureBlkOilGWMethod01
 /////////////////////////////////////////////////////


OCPMixtureBlkOilGWMethod01::OCPMixtureBlkOilGWMethod01(const ParamReservoir& rs_param, const USI& i, OCPMixtureVarSet& vs)
{
    vs.Init(2, 2, OCPMixtureType::BO_GW);


    PVTCO2.Setup(rs_param.PVTCO2.data[i]);
    PVTH2O.Setup(rs_param.PVTH2O.data[i]);

    garciaw.Setup(rs_param.GARCIAW);

    // calculate stdVg and stdVw
    OCP_DBL dummy;
    PVTCO2.CalRhoMuSol(rs_param.Psurf, rs_param.Tsurf, stdVg, dummy, dummy);
    stdVg = 1 / stdVg;
    OCP_DBL xGw;
    PVTH2O.CalRhoMuSol(rs_param.Psurf, rs_param.Tsurf, stdVw, dummy, xGw);
    garciaw.CalRho(rs_param.Tsurf, xGw, stdVw);
    stdVw = 1 / stdVw;

}


void OCPMixtureBlkOilGWMethod01::CalNi(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
    OCP_DBL dummy;
    OCP_DBL xWg, xGw;
    PVTCO2.CalRhoMuSol(vs.P, vs.T, vs.rho[0], dummy, xWg);
    PVTH2O.CalRhoMuSol(vs.P, vs.T, vs.rho[1], dummy, xGw);
    // correct the water phase density
    garciaw.CalRho(vs.T, xGw, vs.rho[1]);

    vs.Ni[0] = Vp * (vs.S[0] * vs.rho[0] * (1 - xWg) + vs.S[1] * vs.rho[1] * xGw);
    vs.Ni[1] = Vp * (vs.S[0] * vs.rho[0] * xWg + vs.S[1] * vs.rho[1] * (1 - xGw));
}


void OCPMixtureBlkOilGWMethod01::InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
    CalNi(Vp, vs);
    Flash(vs);
}


void OCPMixtureBlkOilGWMethod01::Flash(OCPMixtureVarSet& vs)
{  
    FlashDer(vs);
}


void OCPMixtureBlkOilGWMethod01::InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
    CalNi(Vp, vs);
    FlashDer(vs);
}


void OCPMixtureBlkOilGWMethod01::FlashDer(OCPMixtureVarSet& vs)
{

    fill(vs.rhoP.begin(), vs.rhoP.end(), 0.0);
    fill(vs.xiP.begin(), vs.xiP.end(), 0.0);
    fill(vs.muP.begin(), vs.muP.end(), 0.0);
    fill(vs.rhox.begin(), vs.rhox.end(), 0.0);
    fill(vs.xix.begin(), vs.xix.end(), 0.0);
    fill(vs.mux.begin(), vs.mux.end(), 0.0);
    fill(vs.dXsdXp.begin(), vs.dXsdXp.end(), 0.0);

    // xWg, xGw, d xWg / dP, d xGw / dP
    OCP_DBL xWg, xGw;
    OCP_DBL xWgP, xGwP;

    // Gas Properties
    PVTCO2.CalRhoMuSolDer(vs.P, vs.T, vs.rho[0], vs.mu[0], xWg, vs.rhoP[0], vs.muP[0], xWgP);

    // Water Properties
    PVTH2O.CalRhoMuSolDer(vs.P, vs.T, vs.rho[1], vs.mu[1], xGw, vs.rhoP[1], vs.muP[1], xGwP);


    vs.Nt = vs.Ni[0] + vs.Ni[1];
    OCP_DBL dummy;
    const OCP_DBL rgw = vs.Ni[0] / vs.Nt;
    if (rgw < xGw) {
        // water is unsaturated
        vs.phaseExist[0] = OCP_FALSE;
        vs.phaseExist[1] = OCP_TRUE;

        vs.x[1 * 2 + 0] = rgw;
        vs.x[1 * 2 + 1] = 1 - rgw;

        // correct the water phase density
        garciaw.CalRhoDer(vs.T, rgw, 0, vs.rho[1], vs.rhoP[1], vs.rhox[1 * 2 + 0]);

        vs.xi  = vs.rho;
        vs.xiP = vs.rhoP;
        vs.xix = vs.rhox;
    
        vs.vj[0] = 0;
        vs.vj[1] = vs.Nt / vs.rho[1];
        vs.Vf    = vs.vj[1];
        vs.S[0]  = 0.0;
        vs.S[1]  = 1.0;

        vs.vjP[1] = -vs.Ni[1] * vs.rhoP[1] / (vs.rho[1] * vs.rho[1]);
        vs.vfP    = vs.vjP[1];

        vs.vji[1][0] = (vs.rho[1] + vs.rhox[1 * 2 + 0] * vs.Ni[1] / vs.Nt) / (vs.rho[1] * vs.rho[1]);
        vs.vji[1][1] = (vs.rho[1] - vs.rhox[1 * 2 + 0] * vs.Ni[0] / vs.Nt) / (vs.rho[1] * vs.rho[1]);
        vs.vfi[0]    = vs.vji[1][0];
        vs.vfi[1]    = vs.vji[1][1];

        vs.dXsdXp[4 * 3 + 1] = vs.Ni[1] / pow(vs.Nt, 2);   // d XGw / d Ng
        vs.dXsdXp[4 * 3 + 2] = -vs.Ni[0] / pow(vs.Nt, 2);  // d XGw / d Nw
        vs.dXsdXp[5 * 3 + 2] = -vs.dXsdXp[4 * 3 + 1];      // d XWw / d Ng
        vs.dXsdXp[5 * 3 + 2] = -vs.dXsdXp[4 * 3 + 2];      // d XWw / d Nw

    }
    else {
        // water is saturated and gas is always saturated(now)
        vs.phaseExist[0] = OCP_TRUE;
        vs.phaseExist[1] = OCP_TRUE;

        vs.x[0 * 2 + 0] = 1 - xWg;
        vs.x[0 * 2 + 1] = xWg;
        vs.x[1 * 2 + 0] = xGw;
        vs.x[1 * 2 + 1] = 1 - xGw;

        // correct the water phase density
        garciaw.CalRhoDer(vs.T, xGw, xGwP, vs.rho[1], vs.rhoP[1], vs.rhox[1 * 2 + 0]);

        vs.xi  = vs.rho;
        vs.xiP = vs.rhoP;
        vs.xix = vs.rhox;

        // calculate phase mass first
        vs.nj[0] = (vs.Ni[1] - vs.Nt * (1 - xGw)) / (xWg - (1 - xGw));
        vs.nj[1] = vs.Nt - vs.nj[0];

        vs.vj[0]  = vs.nj[0] / vs.rho[0];
        vs.vj[1]  = vs.nj[1] / vs.rho[1];
        vs.Vf     = vs.vj[0] + vs.vj[1];
        vs.S[0]   = vs.vj[0] / vs.Vf;
        vs.S[1]   = vs.vj[1] / vs.Vf;

        const OCP_DBL n0P = (vs.Nt * xGwP - vs.nj[0] * (xWgP + xGwP)) / (xWg - (1 - xGw));
        const OCP_DBL n1P = -n0P;
        vs.vjP[0]         = (n0P * vs.rho[0] + vs.nj[0] * vs.rhoP[0]) / (vs.rho[0] * vs.rho[0]);
        vs.vjP[1]         = (n1P * vs.rho[1] + vs.nj[1] * vs.rhoP[1]) / (vs.rho[1] * vs.rho[1]);
        vs.vfP            = vs.vjP[0] + vs.vjP[1];

        const OCP_DBL n0N0 = -(1 - xGw) / (xWg - (1 - xGw));
        const OCP_DBL n0N1 = xGw / (xWg - (1 - xGw));
        const OCP_DBL n1N0 = 1 - n0N0;
        const OCP_DBL n1N1 = 1 - n0N1;
        vs.vji[0][0]       = n0N0 / vs.rho[0];
        vs.vji[0][1]       = n0N1 / vs.rho[0];
        vs.vji[1][0]       = n1N0 / vs.rho[1];
        vs.vji[1][1]       = n1N1 / vs.rho[1];

        vs.vfi[0] = vs.vji[0][0] + vs.vji[1][0];
        vs.vfi[1] = vs.vji[0][1] + vs.vji[1][1];

        vs.dXsdXp[0 * 3 + 0] = (vs.vjP[0] - vs.S[0] * vs.vfP) / vs.Vf;           // dSg / dP
        vs.dXsdXp[0 * 3 + 1] = (vs.vji[0][0] - vs.S[0] * vs.vfi[0]) / vs.Vf;     // dSg / dNg
        vs.dXsdXp[0 * 3 + 2] = (vs.vji[0][1] - vs.S[0] * vs.vfi[1]) / vs.Vf;     // dSg / dNw
                      
        vs.dXsdXp[1 * 3 + 0] = (vs.vjP[1] - vs.S[1] * vs.vfP) / vs.Vf;           // dSw / dP  
        vs.dXsdXp[1 * 3 + 1] = (vs.vji[1][0] - vs.S[1] * vs.vfi[0]) / vs.Vf;     // dSw / dNg
        vs.dXsdXp[1 * 3 + 2] = (vs.vji[1][1] - vs.S[1] * vs.vfi[1]) / vs.Vf;     // dSw / dNw
                      
        vs.dXsdXp[2 * 3 + 0] = -xWgP;                                            // dXGg / dP
        vs.dXsdXp[3 * 3 + 0] = xWgP;                                             // dXWg / dP
        vs.dXsdXp[4 * 3 + 0] = xGw;                                              // dXGw / dP
        vs.dXsdXp[5 * 3 + 0] = -xGw;                                             // dXWw / dP
    }


}


OCP_DBL OCPMixtureBlkOilGWMethod01::CalXi(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt)
{
    if (pt == PhaseType::gas)         return CalXiG(P, T);
    else if (pt == PhaseType::wat)    return CalXiW(P, T);
    else                              OCP_ABORT("WRONG PHASE TYPE");
}


OCP_DBL OCPMixtureBlkOilGWMethod01::CalRho(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt)
{
    if (pt == PhaseType::gas)         return CalRhoG(P, T);
    else if (pt == PhaseType::wat)    return CalRhoW(P, T);
    else                              OCP_ABORT("WRONG PHASE TYPE");
}


OCP_DBL OCPMixtureBlkOilGWMethod01::CalRhoW(const OCP_DBL& P, const OCP_DBL& T) const
{
    OCP_DBL rhow, xGw, dummy;
    PVTH2O.CalRhoMuSol(P, T, rhow, dummy, xGw);
    garciaw.CalRho(T, xGw, rhow);
    return rhow;
}


OCP_DBL OCPMixtureBlkOilGWMethod01::CalVmStd(const PhaseType& pt)
{
    if (pt == PhaseType::gas)       return stdVg;
    else if (pt == PhaseType::wat)  return stdVw;
    else                            OCP_ABORT("Wrong Phase Type!");
}


void OCPMixtureBlkOilGWMethod01::CalVStd(OCPMixtureVarSet& vs)
{
    vs.vj[0] = vs.Ni[0] * stdVg;
    vs.vj[1] = vs.Ni[1] * stdVw;
}


/////////////////////////////////////////////////////
// OCPMixtureBlkOilGW 
/////////////////////////////////////////////////////

void OCPMixtureBlkOilGW::Setup(const ParamReservoir& rs_param, const USI& i)
{
    if (rs_param.PVTCO2.data.size() > 0 && rs_param.PVTH2O.data.size() > 0) {
        pmMethod = new OCPMixtureBlkOilGWMethod01(rs_param, i, vs);
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/22/2023      Create file                          */
/*----------------------------------------------------------------------------*/