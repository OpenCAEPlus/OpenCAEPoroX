/*! \file    OCPMixtureBlkOilOGW.cpp
 *  \brief   OCPMixtureBlkOilOGW class declaration
 *  \author  Shizhe Li
 *  \date    Jul/19/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPMixtureBlkOilOGW.hpp"

 /////////////////////////////////////////////////////
 // OCPMixtureBlkOilOGWMethod01
 /////////////////////////////////////////////////////


OCPMixtureBlkOilOGWMethod01::OCPMixtureBlkOilOGWMethod01(const vector<vector<OCP_DBL>>& PVCOin,
	const vector<vector<OCP_DBL>>& PVDGin,
	const vector<vector<OCP_DBL>>& PVTWin,
	const OCP_DBL& stdRhoO,
	const OCP_DBL& stdRhoG,
	const OCP_DBL& stdRhoW,
	OCPMixtureVarSet& vs)
{
	PVCO.Setup(PVCOin, stdRhoO, stdRhoG);
	PVDG.Setup(PVDGin, stdRhoG);
	PVTW.Setup(PVTWin, stdRhoW);

	vs.phaseExist[2] = OCP_TRUE;

	vs.x[1 * 3 + 0] = 0;
	vs.x[1 * 3 + 1] = 1;
	vs.x[1 * 3 + 2] = 0;
	vs.x[2 * 3 + 0] = 0;
	vs.x[2 * 3 + 1] = 0;
	vs.x[2 * 3 + 2] = 1;
}


void OCPMixtureBlkOilOGWMethod01::CalNi(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
	// water always exists
	if (vs.S[0] < TINY) {
		if (vs.S[1] < TINY) {
			// water
			vs.Ni[0] = 0;
			vs.Ni[1] = 0;
			vs.Ni[2] = Vp * vs.S[2] * PVTW.CalXiW(vs.P);
		}
		else {
			// gas ,water
			vs.Ni[0] = 0;
			vs.Ni[1] = Vp * vs.S[1] * PVDG.CalXiG(vs.P);
			vs.Ni[2] = Vp * vs.S[2] * PVTW.CalXiW(vs.P);
		}
	}
	else if (vs.S[1] < TINY) {
		// oil, water
		const OCP_DBL rs = PVCO.CalRs(vs.Pb);
		vs.Ni[0] = Vp * vs.S[0] * PVCO.CalXiO(vs.P, vs.Pb) / (1 + rs);
		vs.Ni[1] = vs.Ni[0] * rs;
		vs.Ni[2] = Vp * vs.S[2] * PVTW.CalXiW(vs.P);
	}
	else {
		// oil, gas, water
		const OCP_DBL rs = PVCO.CalRs(vs.P);
		vs.Ni[0] = Vp * vs.S[0] * PVCO.CalXiO(vs.P, vs.P) / (1 + rs);
		vs.Ni[1] = Vp * vs.S[1] * PVDG.CalXiG(vs.P) + vs.Ni[0] * rs;
		vs.Ni[2] = Vp * vs.S[2] * PVTW.CalXiW(vs.P);
	}
}


void OCPMixtureBlkOilOGWMethod01::InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
	CalNi(Vp, vs);
	Flash(vs);
}


void OCPMixtureBlkOilOGWMethod01::Flash(OCPMixtureVarSet& vs)
{
	FlashDer(vs);
}


void OCPMixtureBlkOilOGWMethod01::InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
	CalNi(Vp, vs);
	FlashDer(vs);
}


void OCPMixtureBlkOilOGWMethod01::FlashDer(OCPMixtureVarSet& vs)
{

	fill(vs.rhoP.begin(), vs.rhoP.end(), 0.0);
	fill(vs.xiP.begin(), vs.xiP.end(), 0.0);
	fill(vs.muP.begin(), vs.muP.end(), 0.0);
	fill(vs.rhox.begin(), vs.rhox.end(), 0.0);
	fill(vs.xix.begin(), vs.xix.end(), 0.0);
	fill(vs.mux.begin(), vs.mux.end(), 0.0);
	fill(vs.dXsdXp.begin(), vs.dXsdXp.end(), 0.0);

	vs.Nt = vs.Ni[0] + vs.Ni[1] + vs.Ni[2];

	// some derivatives are ignored, which lead to a more roubust result

	const OCP_DBL Rs_sat = PVCO.CalRs(vs.P);
	if (vs.Ni[0] < vs.Nt * TINY) 
	{
		if (vs.Ni[1] <= vs.Ni[0] * Rs_sat) 
		{
			// only water
			vs.phaseExist[0]  = OCP_FALSE;
			vs.phaseExist[1]  = OCP_FALSE;
			vs.phaseExist[2]  = OCP_TRUE;
			vs.S[0]           = 0;
			vs.S[1]           = 0;
			vs.S[2]           = 1;

			// water property
			PVTW.CalRhoXiMuDer(vs.P, vs.rho[2], vs.xi[2], vs.mu[2], vs.rhoP[2], vs.xiP[2], vs.muP[2]);

			// total
			vs.vj[0]          = 0;
			vs.vj[1]          = 0;
			vs.vj[2]          = vs.Ni[2] / vs.xi[2];
			vs.Vf             = vs.vj[2];

			vs.vjP[2]         = -vs.Ni[2] * vs.xiP[2] / (vs.xi[2] * vs.xi[2]);
			vs.vfP            = vs.vjP[2];

			vs.vji[2][2]      = 1 / vs.xi[2];			
			vs.vfi[0]         = 0;   
			vs.vfi[1]         = 0;   
			vs.vfi[2]         = vs.vji[2][2];

			vs.dXsdXp[1 * 4 + 2] = 1 / PVDG.CalXiG(vs.P) / vs.Vf;                  // dSg / dNg
																			       
			vs.dXsdXp[2 * 4 + 0] = (vs.vjP[2] - vs.S[2] * vs.vfP) / vs.Vf;         // dSw / dP
			vs.dXsdXp[2 * 4 + 1] = -vs.S[2] * vs.vfi[0] / vs.Vf;                   // dSw / dNo
			vs.dXsdXp[2 * 4 + 2] = -vs.S[2] * vs.vfi[1] / vs.Vf;                   // dSw / dNg
			vs.dXsdXp[2 * 4 + 3] = (vs.vji[2][2] - vs.S[2] * vs.vfi[2]) / vs.Vf;   // dSw / dNw
		}
		else {
			// dry gas and water
			vs.phaseExist[0] = OCP_FALSE;
			vs.phaseExist[1] = OCP_TRUE;
			vs.phaseExist[2] = OCP_TRUE;

			// gas property
			PVDG.CalRhoXiMuDer(vs.P, vs.rho[1], vs.xi[1], vs.mu[1], vs.rhoP[1], vs.xiP[1], vs.muP[1]);

			// water property
			PVTW.CalRhoXiMuDer(vs.P, vs.rho[2], vs.xi[2], vs.mu[2], vs.rhoP[2], vs.xiP[2], vs.muP[2]);

			// hypothetical Oil property
			OCP_DBL rs, rsP;
			PVCO.CalRhoXiMuRsDer(vs.P, vs.rho[0], vs.xi[0], vs.mu[0], rs, vs.rhoP[0], vs.xiP[0], vs.muP[0], rsP);

			// total
			vs.vj[0]  = 0;
			vs.vj[1]  = (vs.Ni[1] - rs * vs.Ni[0]) / vs.xi[1];
			vs.vj[2]  = vs.Ni[2] / vs.xi[2];

			vs.Vf     = vs.vj[1] + vs.vj[2];
			vs.S[0]   = 0;
			vs.S[1]   = vs.vj[1] / vs.Vf;
			vs.S[2]   = vs.vj[2] / vs.Vf;

			vs.vjP[0] = vs.Ni[0] * (rsP * vs.xi[0] - (1 + rs) * vs.xiP[0]) / (vs.xi[0] * vs.xi[0]);
			vs.vjP[1] = (-rsP * vs.Ni[0] * vs.xi[1] - (vs.Ni[1] - rs * vs.Ni[0]) * vs.xiP[1]) / (vs.xi[1] * vs.xi[1]);
			vs.vjP[2] = -vs.Ni[2] * vs.xiP[2] / (vs.xi[2] * vs.xi[2]);
			vs.vfP    = vs.vjP[0] + vs.vjP[1] + vs.vjP[2];

			vs.vji[0][0] = (1 + rs) / vs.xi[0];
			vs.vji[1][0] = -rs / vs.xi[1];
			vs.vji[1][1] = 1 / vs.xi[1];
			vs.vji[2][2] = 1 / vs.xi[2];
			vs.vfi[0]    = vs.vji[0][0] + vs.vji[1][0];
			vs.vfi[1]    = vs.vji[1][1];
			vs.vfi[2]    = vs.vji[2][2];

			vs.dXsdXp[1] = vs.vji[0][0] / vs.Vf;   // dSo / dNo

			vs.dXsdXp[1 * 4 + 0] = (vs.vjP[1] - vs.S[1] * vs.vfP) / vs.Vf;        // dSg / dP
			vs.dXsdXp[1 * 4 + 1] = (vs.vji[1][0] - vs.S[1] * vs.vfi[0]) / vs.Vf;  // dSg / dNo
			vs.dXsdXp[1 * 4 + 2] = (vs.vji[1][1] - vs.S[1] * vs.vfi[1]) / vs.Vf;  // dSg / dNg
			vs.dXsdXp[1 * 4 + 3] = -vs.S[1] * vs.vfi[2] / vs.Vf;                  // dSg / dNw

			vs.dXsdXp[2 * 4 + 0] = (vs.vjP[2] - vs.S[2] * vs.vfP) / vs.Vf;        // dSw / dP
			vs.dXsdXp[2 * 4 + 1] = -vs.S[2] * vs.vfi[0] / vs.Vf;                  // dSw / dNo
			vs.dXsdXp[2 * 4 + 2] = -vs.S[2] * vs.vfi[1] / vs.Vf;                  // dSw / dNg
			vs.dXsdXp[2 * 4 + 3] = (vs.vji[2][2] - vs.S[2] * vs.vfi[2]) / vs.Vf;  // dSw / dNw

			vs.dXsdXp[3 * 4 + 0]    = -rsP / ((1 + rs) * (1 + rs));               // d Xoo / dP
			vs.dXsdXp[4 * 4 + 0]    = -vs.dXsdXp[3 * 4 + 0];                      // d Xgo / dP
		}
	}
	else if (vs.Ni[1] <= vs.Ni[0] * Rs_sat) {
		// unsaturated oil and water
		vs.phaseExist[0]  = OCP_TRUE;
		vs.phaseExist[1]  = OCP_FALSE;
		vs.phaseExist[2]  = OCP_TRUE;
		vs.x[0 * 3 + 0] = vs.Ni[0] / (vs.Ni[0] + vs.Ni[1]);
		vs.x[0 * 3 + 1] = 1 - vs.x[0 * 3 + 0];

		// oil property
		const OCP_DBL rs = vs.Ni[1] / vs.Ni[0];
		OCP_DBL rhooRs, xioRs, muoRs;
		PVCO.CalRhoXiMuDer(rs, vs.P, vs.rho[0], vs.xi[0], vs.mu[0], vs.rhoP[0], vs.xiP[0], vs.muP[0], rhooRs, xioRs, muoRs);

		// water property
		PVTW.CalRhoXiMuDer(vs.P, vs.rho[2], vs.xi[2], vs.mu[2], vs.rhoP[2], vs.xiP[2], vs.muP[2]);


		// total
		vs.vj[0]  = vs.Ni[0] * (1 + rs) / vs.xi[0];
		vs.vj[1]  = 0;
		vs.vj[2]  = vs.Ni[2] / vs.xi[2];
		vs.Vf     = vs.vj[0] + vs.vj[2];
		vs.S[0]   = vs.vj[0] / vs.Vf;
		vs.S[1]   = 0;
		vs.S[2]   = vs.vj[2] / vs.Vf;

		vs.vjP[0] = -vs.Ni[0] * (1 + rs) * vs.xiP[0] / (vs.xi[0] * vs.xi[0]);
		vs.vjP[2] = -vs.Ni[2] * vs.xiP[2] / (vs.xi[2] * vs.xi[2]);
		vs.vfP    = vs.vjP[0] + vs.vjP[2];

		vs.vji[0][0] = (1 + rs) / vs.xi[0] - (vs.xi[0] - (1 + rs) * xioRs) / (vs.xi[0] * vs.xi[0]) * rs;
		vs.vji[0][1] = (vs.xi[0] - (1 + rs) * xioRs) / (vs.xi[0] * vs.xi[0]);
		vs.vji[2][2] = 1 / vs.xi[2];
		vs.vfi[0]    = vs.vji[0][0];
		vs.vfi[1]    = vs.vji[0][1];
		vs.vfi[2]    = vs.vji[2][2];

		vs.dXsdXp[0]         = (vs.vjP[0] - vs.S[0] * vs.vfP) / vs.Vf;         // dSo / dP
		vs.dXsdXp[1]         = (vs.vji[0][0] - vs.S[0] * vs.vfi[0]) / vs.Vf;   // dSo / dNo
		vs.dXsdXp[2]         = (vs.vji[0][1] - vs.S[0] * vs.vfi[1]) / vs.Vf;   // dSo / dNg
		vs.dXsdXp[3]         = -vs.S[0] / vs.Vf * vs.vfi[2];                   // dSo / dNw

		vs.dXsdXp[2 * 4 + 0] = (vs.vjP[2] - vs.S[2] * vs.vfP) / vs.Vf;         // dSw / dP
		vs.dXsdXp[2 * 4 + 1] = -vs.S[2] * vs.vfi[0] / vs.Vf;                   // dSw / dNo
		vs.dXsdXp[2 * 4 + 2] = -vs.S[2] * vs.vfi[1] / vs.Vf;                   // dSw / dNg
		vs.dXsdXp[2 * 4 + 3] = (vs.vji[2][2] - vs.S[2] * vs.vfi[2]) / vs.Vf;   // dSw / dNw

		vs.dXsdXp[3 * 4 + 1] = vs.Ni[1] / pow((vs.Ni[0] + vs.Ni[1]), 2);       // d Xoo / d No
		vs.dXsdXp[3 * 4 + 2] = -vs.Ni[0] / pow((vs.Ni[0] + vs.Ni[1]), 2);      // d Xoo / d Ng
		vs.dXsdXp[4 * 4 + 1] = -vs.dXsdXp[3 * 4 + 1];                          // d Xgo / d No
		vs.dXsdXp[4 * 4 + 2] = -vs.dXsdXp[3 * 4 + 2];                          // d Xgo / d Ng

		OCP_DBL tmp_new = (1 + rs) * (1 + rs);
		vs.mux[0]       = -muoRs * tmp_new;         // dMuo / dXoo
		vs.mux[1]       = muoRs * tmp_new;          // dMuo / dXgo				        
		vs.xix[0]       = -xioRs * tmp_new;         // dXio / dXoo
		vs.xix[1]       = xioRs * tmp_new;          // dXio / dXgo
		vs.rhox[0]      = -rhooRs * tmp_new;        // dRhoo / dXoo
		vs.rhox[1]      = rhooRs * tmp_new;         // dRhoo / dXgo
	}
	else {
		// saturated oil, gas and water
		vs.phaseExist[0] = OCP_TRUE;
		vs.phaseExist[1] = OCP_TRUE;
		vs.phaseExist[2] = OCP_TRUE;

		// oil property
		OCP_DBL rs, rsP;
		PVCO.CalRhoXiMuRsDer(vs.P, vs.rho[0], vs.xi[0], vs.mu[0], rs, vs.rhoP[0], vs.xiP[0], vs.muP[0], rsP);

		// gas property
		PVDG.CalRhoXiMuDer(vs.P, vs.rho[1], vs.xi[1], vs.mu[1], vs.rhoP[1], vs.xiP[1], vs.muP[1]);

		// water property
		PVTW.CalRhoXiMuDer(vs.P, vs.rho[2], vs.xi[2], vs.mu[2], vs.rhoP[2], vs.xiP[2], vs.muP[2]);
	
		vs.x[0 * 3 + 0] = 1 / (1 + rs);
		vs.x[0 * 3 + 1] = 1 - vs.x[0 * 3 + 0];

		// total
		vs.vj[0]  = vs.Ni[0] * (1 + rs) / vs.xi[0];
		vs.vj[1]  = (vs.Ni[1] - rs * vs.Ni[0]) / vs.xi[1];
		vs.vj[2]  = vs.Ni[2] / vs.xi[2];
		vs.Vf     = vs.vj[0] + vs.vj[1] + vs.vj[2];
		vs.S[0]   = vs.vj[0] / vs.Vf;
		vs.S[1]   = vs.vj[1] / vs.Vf;
		vs.S[2]   = vs.vj[2] / vs.Vf;

		vs.vjP[0] = vs.Ni[0] * (rsP * vs.xi[0] - (1 + rs) * vs.xiP[0]) / (vs.xi[0] * vs.xi[0]);
		vs.vjP[1] = (-rsP * vs.Ni[0] * vs.xi[1] - (vs.Ni[1] - rs * vs.Ni[0]) * vs.xiP[1]) / (vs.xi[1] * vs.xi[1]);
		vs.vjP[2] = -vs.Ni[2] * vs.xiP[2] / (vs.xi[2] * vs.xi[2]);
		vs.vfP    = vs.vjP[0] + vs.vjP[1] + vs.vjP[2];

		vs.vji[0][0] = (1 + rs) / vs.xi[0];
		vs.vji[1][0] = -rs / vs.xi[1];
		vs.vji[1][1] = 1 / vs.xi[1];
		vs.vji[2][2] = 1 / vs.xi[2];
		vs.vfi[0]    = vs.vji[0][0] + vs.vji[1][0];
		vs.vfi[1]    = vs.vji[1][1];
		vs.vfi[2]    = vs.vji[2][2];

		vs.dXsdXp[0 * 4 + 0] = (vs.vjP[0] - vs.S[0] * vs.vfP) / vs.Vf;           // dSo / dP
		vs.dXsdXp[0 * 4 + 1] = (vs.vji[0][0] - vs.S[0] * vs.vfi[0]) / vs.Vf;     // dSo / dNo
		vs.dXsdXp[0 * 4 + 2] = -vs.S[0] * vs.vfi[1] / vs.Vf;                     // dSo / dNg
		vs.dXsdXp[0 * 4 + 3] = -vs.S[0] * vs.vfi[2] / vs.Vf;                     // dSo / dNw

		vs.dXsdXp[1 * 4 + 0] = (vs.vjP[1] - vs.S[1] * vs.vfP) / vs.Vf;           // dSg / dP  
		vs.dXsdXp[1 * 4 + 1] = (vs.vji[1][0] - vs.S[1] * vs.vfi[0]) / vs.Vf;     // dSg / dNo
		vs.dXsdXp[1 * 4 + 2] = (vs.vji[1][1] - vs.S[1] * vs.vfi[1]) / vs.Vf;     // dSg / dNg
		vs.dXsdXp[1 * 4 + 3] = -vs.S[1] * vs.vfi[2] / vs.Vf;                     // dSg / dNw

		vs.dXsdXp[2 * 4 + 0] = (vs.vjP[2] - vs.S[2] * vs.vfP) / vs.Vf;           // dSw / dP
		vs.dXsdXp[2 * 4 + 1] = -vs.S[2] * vs.vfi[0] / vs.Vf;                     // dSw / dNo
		vs.dXsdXp[2 * 4 + 2] = -vs.S[2] * vs.vfi[1] / vs.Vf;                     // dSw / dNg
		vs.dXsdXp[2 * 4 + 3] = (vs.vji[2][2] - vs.S[2] * vs.vfi[2]) / vs.Vf;     // dSw / dNw

		vs.dXsdXp[3 * 4 + 0] = -rsP / ((1 + rs) * (1 + rs));                     // d Xoo / dP
		vs.dXsdXp[4 * 4 + 0] = -vs.dXsdXp[3 * 4 + 0];                            // d Xgo / dP
	}
}


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOGW 
/////////////////////////////////////////////////////


void OCPMixtureBlkOilOGW::Setup(const ParamReservoir& rs_param, const USI& i)
{
	vs.Init(3, 3, OCP_FALSE);
	GetStdRhoOGW(rs_param);
	if (rs_param.PVCO_T.data.size() > 0 &&
		rs_param.PVDG_T.data.size() > 0 &&
		rs_param.PVTW_T.data.size() > 0) {
		pmMethod = new OCPMixtureBlkOilOGWMethod01(
			rs_param.PVCO_T.data[i],
			rs_param.PVDG_T.data[i],
			rs_param.PVTW_T.data[i],
			stdRhoO, stdRhoG, stdRhoW, vs);
	}
}


void OCPMixtureBlkOilOGW::GetStdRhoOGW(const ParamReservoir& rs_param)
{
	if (rs_param.density.activity) {
		stdRhoO = rs_param.density.data[0];
		stdRhoW = rs_param.density.data[1];
		stdRhoG = rs_param.density.data[2];
	}
	else {
		stdRhoO = (141.5 * RHOW_STD) / (rs_param.gravity.data[0] + 131.5);
		stdRhoW = RHOW_STD * rs_param.gravity.data[1];
		stdRhoG = RHOAIR_STD * rs_param.gravity.data[2];
	}
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/19/2023      Create file                          */
/*----------------------------------------------------------------------------*/