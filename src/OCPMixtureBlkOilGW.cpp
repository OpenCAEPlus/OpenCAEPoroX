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
    PVTCO2.Setup(rs_param.PVTCO2.data[i]);
    PVTH2O.Setup(rs_param.PVTH2O.data[i]);

    MWCO2 = rs_param.comsParam.MW.data[0][0];
    MWH2O = rs_param.comsParam.MW.data[0][1];
}


void OCPMixtureBlkOilGWMethod01::InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
    //vs.Ni[0] = Vp * vs.S[0] * PVDO->CalXiO(vs.P);
    //vs.Ni[1] = Vp * vs.S[1] * PVTW.CalXiW(vs.P);

    Flash(vs);
}


void OCPMixtureBlkOilGWMethod01::Flash(OCPMixtureVarSet& vs)
{

    //// Oil Properties
    //PVDO->CalRhoXiMuDer(vs.P, vs.rho[0], vs.xi[0], vs.mu[0], vs.rhoP[0], vs.xiP[0], vs.muP[0]);

    //// Water Properties
    //PVTW.CalRhoXiMuDer(vs.P, vs.rho[1], vs.xi[1], vs.mu[1], vs.rhoP[1], vs.xiP[1], vs.muP[1]);

    vs.Nt = vs.Ni[0] + vs.Ni[1];
    vs.vj[0] = vs.Ni[0] / vs.xi[0];
    vs.vj[1] = vs.Ni[1] / vs.xi[1];
    vs.Vf = vs.vj[0] + vs.vj[1];
    vs.S[0] = vs.vj[0] / vs.Vf;
    vs.S[1] = vs.vj[1] / vs.Vf;

    vs.vji[0][0] = 1 / vs.xi[0];
    vs.vji[1][1] = 1 / vs.xi[1];
    vs.vjP[0] = -vs.Ni[0] * vs.xiP[0] / (vs.xi[0] * vs.xi[0]);
    vs.vjP[1] = -vs.Ni[1] * vs.xiP[1] / (vs.xi[1] * vs.xi[1]);

    vs.vfi[0] = vs.vji[0][0];
    vs.vfi[1] = vs.vji[1][1];
    vs.vfP = vs.vjP[0] + vs.vjP[1];
}


void OCPMixtureBlkOilGWMethod01::InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
    //vs.Ni[0] = Vp * vs.S[0] * PVDO->CalXiO(vs.P);
    //vs.Ni[1] = Vp * vs.S[1] * PVTW.CalXiW(vs.P);

    FlashDer(vs);
}


void OCPMixtureBlkOilGWMethod01::FlashDer(OCPMixtureVarSet& vs)
{
    Flash(vs);

    vs.dXsdXp[0] = (vs.vjP[0] - vs.S[0] * vs.vfP) / vs.Vf;        // dSo / dP
    vs.dXsdXp[1] = (vs.vji[0][0] - vs.S[0] * vs.vfi[0]) / vs.Vf;  // dSo / dNo
    vs.dXsdXp[2] = -vs.S[0] * vs.vfi[1] / vs.Vf;                  // dSo / dNw

    vs.dXsdXp[3] = -vs.dXsdXp[0];  // dSw / dP
    vs.dXsdXp[4] = -vs.dXsdXp[1];  // dSw / dNo
    vs.dXsdXp[5] = -vs.dXsdXp[2];  // dSw / dNw
}


OCP_DBL OCPMixtureBlkOilGWMethod01::CalXi(const OCP_DBL& P, const PhaseType& pt)
{
    if (pt == PhaseType::oil)         return CalXiO(P);
    else if (pt == PhaseType::wat)    return CalXiW(P);
    else                              OCP_ABORT("WRONG PHASE TYPE");
}


OCP_DBL OCPMixtureBlkOilGWMethod01::CalRho(const OCP_DBL& P, const PhaseType& pt)
{
    if (pt == PhaseType::oil)         return CalRhoO(P);
    else if (pt == PhaseType::wat)    return CalRhoW(P);
    else                              OCP_ABORT("WRONG PHASE TYPE");
}


OCP_DBL OCPMixtureBlkOilGWMethod01::CalVmStd(const PhaseType& pt)
{
    //if (pt == PhaseType::oil)  return (stdVo * CONV1);
    //else if (pt == PhaseType::wat)  return (stdVw * CONV1);
    //else    OCP_ABORT("Wrong Phase Type!");
}


void OCPMixtureBlkOilGWMethod01::CalVStd(OCPMixtureVarSet& vs)
{
    //vs.vj[0] = vs.Ni[0] * stdVo * CONV1;
    //vs.vj[1] = vs.Ni[1] * stdVw * CONV1;
}


/////////////////////////////////////////////////////
// OCPMixtureBlkOilGW 
/////////////////////////////////////////////////////

void OCPMixtureBlkOilGW::Setup(const ParamReservoir& rs_param, const USI& i)
{
    vs.Init(2, 2, mixtureType);
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