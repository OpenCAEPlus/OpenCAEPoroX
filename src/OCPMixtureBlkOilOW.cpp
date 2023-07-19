/*! \file    OCPMixtureBlkOilOW.cpp
 *  \brief   OCPMixtureBlkOilOW class declaration
 *  \author  Shizhe Li
 *  \date    Jul/13/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPMixtureBlkOilOW.hpp"


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOWMethod01
/////////////////////////////////////////////////////


OCPMixtureBlkOilOWMethod01::OCPMixtureBlkOilOWMethod01(const vector<vector<OCP_DBL>>& PVDOin,
    const vector<vector<OCP_DBL>>& PVTWin,
    const OCP_DBL& stdRhoO,
    const OCP_DBL& stdRhoW,
    OCPMixtureVarSet& vs) 
{
    PVDO.Setup(PVDOin, stdRhoO);
    PVTW.Setup(PVTWin, stdRhoW);

    vs.phaseExist[0]  = OCP_TRUE;
    vs.phaseExist[1]  = OCP_TRUE;

    vs.xij[0 * 2 + 0] = 1;
    vs.xij[0 * 2 + 1] = 0;
    vs.xij[1 * 2 + 0] = 0;
    vs.xij[1 * 2 + 1] = 1;
}


void OCPMixtureBlkOilOWMethod01::InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
    vs.Ni[0] = Vp * (1 - vs.S[1]) * PVDO.CalXiO(vs.P);
    vs.Ni[1] = Vp * vs.S[1] * PVTW.CalXiW(vs.P);

    Flash(vs);
}


void OCPMixtureBlkOilOWMethod01::Flash(OCPMixtureVarSet& vs)
{

    // Oil Properties
    PVDO.CalRhoXiMuDer(vs.P, vs.rho[0], vs.xi[0], vs.mu[0], vs.rhoP[0], vs.xiP[0], vs.muP[0]);

    // Water Properties
    PVTW.CalRhoXiMuDer(vs.P, vs.rho[1], vs.xi[1], vs.mu[1], vs.rhoP[1], vs.xiP[1], vs.muP[1]); 

    vs.Nt             = vs.Ni[0] + vs.Ni[1];
    vs.vj[0]          = vs.Ni[0] / vs.xi[0];
    vs.vj[1]          = vs.Ni[1] / vs.xi[1];
    vs.vf             = vs.vj[0] + vs.vj[1];
    vs.S[0]           = vs.vj[0] / vs.vf;
    vs.S[1]           = vs.vj[1] / vs.vf;

    vs.vji[0 * 2 + 0] = 1 / vs.xi[0];
    vs.vji[1 * 2 + 1] = 1 / vs.xi[1];
    vs.vjP[0]         = -vs.Ni[0] * vs.xiP[0] / (vs.xi[0] * vs.xi[0]);
    vs.vjP[1]         = -vs.Ni[1] * vs.xiP[1] / (vs.xi[1] * vs.xi[1]);

    vs.vfi[0]         = vs.vji[0 * 2 + 0];
    vs.vfi[1]         = vs.vji[1 * 2 + 1];
    vs.vfP            = vs.vjP[0] + vs.vjP[1];
}


void OCPMixtureBlkOilOWMethod01::InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
    vs.Ni[0] = Vp * (1 - vs.S[1]) * PVDO.CalXiO(vs.P);
    vs.Ni[1] = Vp * vs.S[1] * PVTW.CalXiW(vs.P);
    
    FlashDer(vs);
}


void OCPMixtureBlkOilOWMethod01::FlashDer(OCPMixtureVarSet& vs)
{
    Flash(vs);

    vs.dXsdXp[0]      = (vs.vjP[0] - vs.S[0] * vs.vfP) / vs.vf;             // dSo / dP
    vs.dXsdXp[1]      = (vs.vji[0 * 2 + 0] - vs.S[0] * vs.vfi[0]) / vs.vf;  // dSo / dNo
    vs.dXsdXp[2]      = -vs.S[0] * vs.vfi[1] / vs.vf;                       // dSo / dNw
                      
    vs.dXsdXp[3]      = (vs.vjP[1] - vs.S[1] * vs.vfP) / vs.vf;             // dSw / dP
    vs.dXsdXp[4]      = -vs.S[1] * vs.vfi[0] / vs.vf;                       // dSw / dNo
    vs.dXsdXp[5]      = (vs.vji[1 * 2 + 1] - vs.S[1] * vs.vfi[1]) / vs.vf;  // dSw / dNw
}


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOW 
/////////////////////////////////////////////////////

void OCPMixtureBlkOilOW::Setup(const ParamReservoir& rs_param, const USI& i)
{  
    vs.Init(2, 2, OCP_FALSE);
    GetStdRhoOW(rs_param);
	if (rs_param.PVTW_T.data.size() > 0 && rs_param.PVDO_T.data.size() > 0) {
		pmMethod = new OCPMixtureBlkOilOWMethod01(rs_param.PVDO_T.data[i], 
                                                  rs_param.PVTW_T.data[i], 
                                                  stdRhoO, stdRhoW, vs);
	}
}


void OCPMixtureBlkOilOW::GetStdRhoOW(const ParamReservoir& rs_param)
{
    if (rs_param.density.activity) {
        stdRhoO = rs_param.density.data[0];
        stdRhoW = rs_param.density.data[1];
    }
    else {
        stdRhoO = (141.5 * RHOW_STD) / (rs_param.gravity.data[0] + 131.5);
        stdRhoW = RHOW_STD * rs_param.gravity.data[1];
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/13/2023      Create file                          */
/*----------------------------------------------------------------------------*/