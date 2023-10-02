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


OCPMixtureBlkOilOGWMethod01::OCPMixtureBlkOilOGWMethod01(const ParamReservoir& rs_param, const USI& i,
	OCPMixtureVarSet& vs)
{
	vs.Init(3, 3, OCPMixtureType::BO_OGW);


	OCP_DBL stdRhoO, stdRhoW, stdRhoG;
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

	PVCO.Setup(rs_param.PVCO_T.data[i], stdRhoO, stdRhoG, stdVo, stdVg);
	PVDG.Setup(rs_param.PVDG_T.data[i], stdRhoG, stdVg);
	PVTW.Setup(rs_param.PVTW_T.data[i], stdRhoW, stdVw);

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
		x  = PVCO.CalRs(vs.Pb) * stdVo / stdVg;
		vs.Ni[0] = Vp * vs.S[0] * PVCO.CalXiO(vs.P, vs.Pb) / (1 + x);
		vs.Ni[1] = vs.Ni[0] * x;
		vs.Ni[2] = Vp * vs.S[2] * PVTW.CalXiW(vs.P);
	}
	else {
		// oil, gas, water
		x = PVCO.CalRs(vs.P) * stdVo / stdVg;
		vs.Ni[0] = Vp * vs.S[0] * PVCO.CalXiO(vs.P, vs.P) / (1 + x);
		vs.Ni[1] = Vp * vs.S[1] * PVDG.CalXiG(vs.P) + vs.Ni[0] * x;
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
	x = PVCO.CalRs(vs.P) * stdVo / stdVg;

	if (vs.Ni[0] < vs.Nt * TINY) 
	{
		if (vs.Ni[1] <= vs.Ni[0] * x)
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

			x = rs * (stdVo / stdVg);
			const OCP_DBL xP = rsP * (stdVo / stdVg);

			// total
			vs.vj[0]  = 0;
			vs.vj[1]  = (vs.Ni[1] - x * vs.Ni[0]) / vs.xi[1];
			vs.vj[2]  = vs.Ni[2] / vs.xi[2];

			vs.Vf     = vs.vj[1] + vs.vj[2];
			vs.S[0]   = 0;
			vs.S[1]   = vs.vj[1] / vs.Vf;
			vs.S[2]   = vs.vj[2] / vs.Vf;

			vs.vjP[0] = vs.Ni[0] * (xP * vs.xi[0] - (1 + x) * vs.xiP[0]) / (vs.xi[0] * vs.xi[0]);
			vs.vjP[1] = (-xP * vs.Ni[0] * vs.xi[1] - (vs.Ni[1] - x * vs.Ni[0]) * vs.xiP[1]) / (vs.xi[1] * vs.xi[1]);
			vs.vjP[2] = -vs.Ni[2] * vs.xiP[2] / (vs.xi[2] * vs.xi[2]);
			vs.vfP    = vs.vjP[0] + vs.vjP[1] + vs.vjP[2];

			vs.vji[0][0] = (1 + x) / vs.xi[0];
			vs.vji[1][0] = -x / vs.xi[1];
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

			vs.dXsdXp[3 * 4 + 0]    = -xP / ((1 + x) * (1 + x));               // d Xoo / dP
			vs.dXsdXp[4 * 4 + 0]    = -vs.dXsdXp[3 * 4 + 0];                      // d Xgo / dP
		}
	}
	else if (vs.Ni[1] <= vs.Ni[0] * x) {
		// unsaturated oil and water
		x = vs.Ni[1] / vs.Ni[0];

		vs.phaseExist[0]  = OCP_TRUE;
		vs.phaseExist[1]  = OCP_FALSE;
		vs.phaseExist[2]  = OCP_TRUE;
		vs.x[0 * 3 + 0]   = 1 / (1 + x);
		vs.x[0 * 3 + 1]   = 1 - vs.x[0 * 3 + 0];

		// oil property
		const OCP_DBL rs = x * (stdVg / stdVo);
		OCP_DBL rhooRs, xioRs, muoRs;
		PVCO.CalRhoXiMuDer(rs, vs.P, vs.rho[0], vs.xi[0], vs.mu[0], vs.rhoP[0], vs.xiP[0], vs.muP[0], rhooRs, xioRs, muoRs);
		const OCP_DBL xiox  = xioRs * (stdVg / stdVo);
		const OCP_DBL rhoox = rhooRs * (stdVg / stdVo);
		const OCP_DBL muox  = muoRs * (stdVg / stdVo);


		// water property
		PVTW.CalRhoXiMuDer(vs.P, vs.rho[2], vs.xi[2], vs.mu[2], vs.rhoP[2], vs.xiP[2], vs.muP[2]);

		// total
		vs.vj[0]  = (vs.Ni[0] + vs.Ni[1]) / vs.xi[0];
		vs.vj[1]  = 0;
		vs.vj[2]  = vs.Ni[2] / vs.xi[2];
		vs.Vf     = vs.vj[0] + vs.vj[2];
		vs.S[0]   = vs.vj[0] / vs.Vf;
		vs.S[1]   = 0;
		vs.S[2]   = vs.vj[2] / vs.Vf;

		vs.vjP[0] = -(vs.Ni[0] + vs.Ni[1]) * vs.xiP[0] / (vs.xi[0] * vs.xi[0]);
		vs.vjP[2] = -vs.Ni[2] * vs.xiP[2] / (vs.xi[2] * vs.xi[2]);
		vs.vfP    = vs.vjP[0] + vs.vjP[2];

		vs.vji[0][0] = (vs.xi[0] + x * (1 + x) * xiox) / (vs.xi[0] * vs.xi[0]);
		vs.vji[0][1] = (vs.xi[0] - (1 + x) * xiox) / (vs.xi[0] * vs.xi[0]);
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

		const OCP_DBL tmp_new = (1 + x) * (1 + x);
		vs.mux[0]       = -muox * tmp_new;         // dMuo / dXoo
		vs.mux[1]       = muox * tmp_new;          // dMuo / dXgo				        
		vs.xix[0]       = -xiox * tmp_new;         // dXio / dXoo
		vs.xix[1]       = xiox * tmp_new;          // dXio / dXgo
		vs.rhox[0]      = -rhoox * tmp_new;        // dRhoo / dXoo
		vs.rhox[1]      = rhoox * tmp_new;         // dRhoo / dXgo
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
		
		x = rs * (stdVo / stdVg);
		const OCP_DBL xP = rsP * (stdVo / stdVg);

		vs.x[0 * 3 + 0] = 1 / (1 + x);
		vs.x[0 * 3 + 1] = 1 - vs.x[0 * 3 + 0];

		// total
		vs.vj[0]  = vs.Ni[0] * (1 + x) / vs.xi[0];
		vs.vj[1]  = (vs.Ni[1] - x * vs.Ni[0]) / vs.xi[1];
		vs.vj[2]  = vs.Ni[2] / vs.xi[2];
		vs.Vf     = vs.vj[0] + vs.vj[1] + vs.vj[2];
		vs.S[0]   = vs.vj[0] / vs.Vf;
		vs.S[1]   = vs.vj[1] / vs.Vf;
		vs.S[2]   = vs.vj[2] / vs.Vf;

		vs.vjP[0] = vs.Ni[0] * (xP * vs.xi[0] - (1 + x) * vs.xiP[0]) / (vs.xi[0] * vs.xi[0]);
		vs.vjP[1] = (-xP * vs.Ni[0] * vs.xi[1] - (vs.Ni[1] - x * vs.Ni[0]) * vs.xiP[1]) / (vs.xi[1] * vs.xi[1]);
		vs.vjP[2] = -vs.Ni[2] * vs.xiP[2] / (vs.xi[2] * vs.xi[2]);
		vs.vfP    = vs.vjP[0] + vs.vjP[1] + vs.vjP[2];

		vs.vji[0][0] = (1 + x) / vs.xi[0];
		vs.vji[1][0] = -x / vs.xi[1];
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

		vs.dXsdXp[3 * 4 + 0] = -xP / ((1 + x) * (1 + x));                        // d Xoo / dP
		vs.dXsdXp[4 * 4 + 0] = -vs.dXsdXp[3 * 4 + 0];                            // d Xgo / dP
	}
}


OCP_DBL OCPMixtureBlkOilOGWMethod01::CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const PhaseType& pt)
{
	if (pt == PhaseType::oil)         return CalXiO(P, Pb);
	else if (pt == PhaseType::gas)    return CalXiG(P);
	else if (pt == PhaseType::wat)  return CalXiW(P);
	else                              OCP_ABORT("WRONG PHASE TYPE");
}


OCP_DBL OCPMixtureBlkOilOGWMethod01::CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const PhaseType& pt)
{
	if (pt == PhaseType::oil)         return CalRhoO(P, Pb);
	else if (pt == PhaseType::gas)    return CalRhoG(P);
	else if (pt == PhaseType::wat)  return CalRhoW(P);
	else                              OCP_ABORT("WRONG PHASE TYPE");
}


OCP_DBL OCPMixtureBlkOilOGWMethod01::CalVmStd(const PhaseType& pt)
{
	if      (pt == PhaseType::oil)  return (stdVo * CONV1);
    else if (pt == PhaseType::gas)  return (stdVg * 1000);
    else if (pt == PhaseType::wat)  return (stdVw * CONV1);
    else    OCP_ABORT("Wrong Phase Type!");
}


void OCPMixtureBlkOilOGWMethod01::CalVStd(OCPMixtureVarSet& vs)
{
	vs.vj[0] = vs.Ni[0] * stdVo * CONV1;
	vs.vj[1] = vs.Ni[1] * stdVg * 1000;
	vs.vj[2] = vs.Ni[2] * stdVw * CONV1;
}


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOGW 
/////////////////////////////////////////////////////


void OCPMixtureBlkOilOGW::Setup(const ParamReservoir& rs_param, const USI& i)
{
	if (rs_param.PVCO_T.data.size() > 0 &&
		rs_param.PVDG_T.data.size() > 0 &&
		rs_param.PVTW_T.data.size() > 0) {
		pmMethod = new OCPMixtureBlkOilOGWMethod01(rs_param, i, vs);
	}
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/19/2023      Create file                          */
/*----------------------------------------------------------------------------*/