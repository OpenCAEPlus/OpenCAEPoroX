/*! \file    OCPMixtureKMethod.cpp
 *  \brief   OCPMixtureKMethod class declaration
 *  \author  Shizhe Li
 *  \date    Oct/02/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


#include "OCPMixtureMethodK.hpp"


/////////////////////////////////////////////////////
// OCPMixtureMethodK_OW01
/////////////////////////////////////////////////////


OCPMixtureMethodK_OW01::OCPMixtureMethodK_OW01(const ParamReservoir& rs_param, const USI& i, OCPMixtureVarSet& vs) 
{

    vs.Init(2, 2, OCPMixtureType::BO_OW);

    OCP_DBL stdRhoO, stdRhoW;
    if (rs_param.density.activity) {
        stdRhoO = rs_param.density.data[0];
        stdRhoW = rs_param.density.data[1];
    }
    else {
        stdRhoO = (141.5 * RHOW_STD) / (rs_param.gravity.data[0] + 131.5);
        stdRhoW = RHOW_STD * rs_param.gravity.data[1];
    }
    if (rs_param.PVDO_T.data.size() > 0) {
        PVDO = new OCP_PVDO();
        PVDO->Setup(rs_param.PVDO_T.data[i], stdRhoO, stdVo);
    }
    else if (rs_param.PVCDO_T.data.size() > 0) {
        PVDO = new OCP_PVCDO();
        PVDO->Setup(rs_param.PVCDO_T.data[i], stdRhoO, stdVo);
    }
    else {
        OCP_ABORT("PARAMS are not enough!");
    }
    
    PVTW.Setup(rs_param.PVTW_T.data[i], stdRhoW, stdVw);

    vs.phaseExist[0]  = OCP_TRUE;
    vs.phaseExist[1]  = OCP_TRUE;

    vs.x[0 * 2 + 0] = 1;
    vs.x[0 * 2 + 1] = 0;
    vs.x[1 * 2 + 0] = 0;
    vs.x[1 * 2 + 1] = 1;
}


void OCPMixtureMethodK_OW01::SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const
{
	mvs.P = bvs.P[bId];
	mvs.T = bvs.T[bId];
	copy(&bvs.Ni[bId * bvs.nc], &bvs.Ni[bId * bvs.nc] + bvs.nc, mvs.Ni.begin());
	copy(&bvs.S[bId * bvs.np], &bvs.S[bId * bvs.np] + bvs.np, mvs.S.begin());
}


void OCPMixtureMethodK_OW01::InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
    vs.Ni[0] = Vp * vs.S[0] * PVDO->CalXiO(vs.P);
    vs.Ni[1] = Vp * vs.S[1] * PVTW.CalXiW(vs.P);

    Flash(vs);
}


void OCPMixtureMethodK_OW01::Flash(OCPMixtureVarSet& vs)
{

    // Oil Properties
    PVDO->CalRhoXiMuDer(vs.P, vs.rho[0], vs.xi[0], vs.mu[0], vs.rhoP[0], vs.xiP[0], vs.muP[0]);

    // Water Properties
    PVTW.CalRhoXiMuDer(vs.P, vs.rho[1], vs.xi[1], vs.mu[1], vs.rhoP[1], vs.xiP[1], vs.muP[1]); 

    vs.Nt             = vs.Ni[0] + vs.Ni[1];
    vs.vj[0]          = vs.Ni[0] / vs.xi[0];
    vs.vj[1]          = vs.Ni[1] / vs.xi[1];
    vs.Vf             = vs.vj[0] + vs.vj[1];
    vs.S[0]           = vs.vj[0] / vs.Vf;
    vs.S[1]           = vs.vj[1] / vs.Vf;

    vs.vji[0][0]      = 1 / vs.xi[0];
    vs.vji[1][1]      = 1 / vs.xi[1];
    vs.vjP[0]         = -vs.Ni[0] * vs.xiP[0] / (vs.xi[0] * vs.xi[0]);
    vs.vjP[1]         = -vs.Ni[1] * vs.xiP[1] / (vs.xi[1] * vs.xi[1]);

    vs.vfi[0]         = vs.vji[0][0];
    vs.vfi[1]         = vs.vji[1][1];
    vs.vfP            = vs.vjP[0] + vs.vjP[1];
}


void OCPMixtureMethodK_OW01::InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
    vs.Ni[0] = Vp * vs.S[0] * PVDO->CalXiO(vs.P);
    vs.Ni[1] = Vp * vs.S[1] * PVTW.CalXiW(vs.P);
    
    FlashDer(vs);
}


void OCPMixtureMethodK_OW01::FlashDer(OCPMixtureVarSet& vs)
{
    Flash(vs);

    vs.dXsdXp[0] = (vs.vjP[0] - vs.S[0] * vs.vfP) / vs.Vf;        // dSo / dP
    vs.dXsdXp[1] = (vs.vji[0][0] - vs.S[0] * vs.vfi[0]) / vs.Vf;  // dSo / dNo
    vs.dXsdXp[2] = -vs.S[0] * vs.vfi[1] / vs.Vf;                  // dSo / dNw
                 
    vs.dXsdXp[3] = -vs.dXsdXp[0];  // dSw / dP
    vs.dXsdXp[4] = -vs.dXsdXp[1];  // dSw / dNo
    vs.dXsdXp[5] = -vs.dXsdXp[2];  // dSw / dNw
}


OCP_DBL OCPMixtureMethodK_OW01::CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt)
{
    if (pt == PhaseType::oil)         return CalXiO(P);
    else if (pt == PhaseType::wat)    return CalXiW(P);
    else                              OCP_ABORT("WRONG PHASE TYPE");
}


OCP_DBL OCPMixtureMethodK_OW01::CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt)
{
    if (pt == PhaseType::oil)         return CalRhoO(P);
    else if (pt == PhaseType::wat)    return CalRhoW(P);
    else                              OCP_ABORT("WRONG PHASE TYPE");
}


OCP_DBL OCPMixtureMethodK_OW01::CalVmStd(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt)
{
    if      (pt == PhaseType::oil)  return (stdVo * CONV1);
    else if (pt == PhaseType::wat)  return (stdVw * CONV1);
    else    OCP_ABORT("Wrong Phase Type!");
}


void OCPMixtureMethodK_OW01::CalVStd(OCPMixtureVarSet& vs)
{
    vs.vj[0] = vs.Ni[0] * stdVo * CONV1;
    vs.vj[1] = vs.Ni[1] * stdVw * CONV1;
}



/////////////////////////////////////////////////////
// OCPMixtureMethodK_OW01T
/////////////////////////////////////////////////////


OCPMixtureMethodK_OW01T::OCPMixtureMethodK_OW01T(const ComponentParam& param, const USI& tarId, OCPMixtureVarSet& vs)
{

	vs.Init(2, 2, OCPMixtureType::THERMALK_OW);

	if (param.molden.activity)
		xi_ref = param.molden.data[tarId];
	else
		OCP_ABORT("ACF hasn't been input!");

	if (param.cp.activity)
		cp = param.cp.data[tarId];
	else
		cp.resize(2, 0);

	if (param.ct1.activity)
		ct1 = param.ct1.data[tarId];
	else
		ct1.resize(2, 0);

	if (param.ct2.activity)
		ct2 = param.ct2.data[tarId];
	else
		ct2.resize(2, 0);

	if (param.cpt.activity)
		cpt = param.cpt.data[tarId];
	else
		cpt.resize(2, 0);

	if (param.MW.activity)
		MWc = param.MW.data[tarId];
	else
		OCP_ABORT("MW hasn't been input!");

	Tref = param.Tref[tarId] + CONV5;
	Pref = param.Pref[tarId];

	vC.Setup(param, tarId);
	eC.Setup(param, tarId);

	MWp = MWc;

	// Init
	vs.phaseExist[0] = true;
	vs.phaseExist[1] = true;

	vs.x[0 * 2 + 0] = 1;
	vs.x[0 * 2 + 1] = 0;
	vs.x[1 * 2 + 0] = 0;
	vs.x[1 * 2 + 1] = 1;

	// d mu / dP
	fill(vs.rhox.begin(), vs.rhox.end(), 0.0);
	fill(vs.xix.begin(), vs.xix.end(), 0.0);
	fill(vs.muP.begin(), vs.muP.end(), 0.0);
	fill(vs.mux.begin(), vs.mux.end(), 0.0);
}


void OCPMixtureMethodK_OW01T::SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const
{
	mvs.P = bvs.P[bId];
	mvs.T = bvs.T[bId] + CONV5;
	copy(&bvs.Ni[bId * bvs.nc], &bvs.Ni[bId * bvs.nc] + bvs.nc, mvs.Ni.begin());
	copy(&bvs.S[bId * bvs.np], &bvs.S[bId * bvs.np] + bvs.np, mvs.S.begin());
}



OCP_DBL OCPMixtureMethodK_OW01T::CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt)
{
	if (pt == PhaseType::oil)       return CalRhoO(P, T);
	else if (pt == PhaseType::wat)  return CalRhoW(P, T);
	else                            OCP_ABORT("WRONG TarPhase");
}


OCP_DBL OCPMixtureMethodK_OW01T::CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt)
{
	if (pt == PhaseType::oil)       return CalXiO(P, T);
	else if (pt == PhaseType::wat)  return CalXiW(P, T);
	else                            OCP_ABORT("WRONG TarPhase");
}


OCP_DBL OCPMixtureMethodK_OW01T::CalXiO(const OCP_DBL& P, const OCP_DBL& T)
{
	const OCP_DBL dP = P - Pref;
	const OCP_DBL dT = T - Tref;
	return xi_ref[0] * exp(cp[0] * dP - ct1[0] * dT - ct2[0] * (pow(T, 2) - pow(Tref, 2)) / 2 + cpt[0] * dP * dT);
}


OCP_DBL OCPMixtureMethodK_OW01T::CalXiW(const OCP_DBL& P, const OCP_DBL& T)
{
	const OCP_DBL dP = P - Pref;
	const OCP_DBL dT = T - Tref;
	return xi_ref[1] * exp(cp[1] * dP - ct1[1] * dT - ct2[1] * (pow(T, 2) - pow(Tref, 2)) / 2 + cpt[1] * dP * dT);
}


OCP_DBL OCPMixtureMethodK_OW01T::CalRhoO(const OCP_DBL& P, const OCP_DBL& T)
{
	return MWp[0] * CalXiO(P, T);
}


OCP_DBL OCPMixtureMethodK_OW01T::CalRhoW(const OCP_DBL& P, const OCP_DBL& T)
{
	return MWp[1] * CalXiW(P, T);
}


void OCPMixtureMethodK_OW01T::InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
	vs.Ni[0] = Vp * vs.S[0] * CalXiO(vs.P, vs.T);
	vs.Ni[1] = Vp * vs.S[1] * CalXiW(vs.P, vs.T);
	Flash(vs);
}


void OCPMixtureMethodK_OW01T::Flash(OCPMixtureVarSet& vs)
{
	FlashDer(vs);
}


void OCPMixtureMethodK_OW01T::InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
	vs.Ni[0] = Vp * vs.S[0] * CalXiO(vs.P, vs.T);
	vs.Ni[1] = Vp * vs.S[1] * CalXiW(vs.P, vs.T);
	FlashDer(vs);
}


void OCPMixtureMethodK_OW01T::FlashDer(OCPMixtureVarSet& vs)
{
	// Assign value
	const OCP_DBL dP = vs.P - Pref;
	const OCP_DBL dT = vs.T - Tref;

	vs.Nt = vs.Ni[0] + vs.Ni[1];

	// phase viscosity
	vs.mu[0] = vC.CalViscosity(ViscosityParams(&vs.P, &vs.T, &vs.x[0 * 2]), vs.muP[0], vs.muT[0], &vs.mux[0 * 2]);
	vs.mu[1] = vC.CalViscosity(ViscosityParams(&vs.P, &vs.T, &vs.x[1 * 2]), vs.muP[1], vs.muT[1], &vs.mux[1 * 2]);

	// phase molar density
	vs.xi[0] = xi_ref[0] * exp(cp[0] * dP - ct1[0] * dT - ct2[0] * (pow(vs.T, 2) - pow(Tref, 2)) / 2 + cpt[0] * dP * dT);
	vs.xi[1] = xi_ref[1] * exp(cp[1] * dP - ct1[1] * dT - ct2[1] * (pow(vs.T, 2) - pow(Tref, 2)) / 2 + cpt[1] * dP * dT);
	// d xi / dP
	vs.xiP[0] = vs.xi[0] * (cp[0] + cpt[0] * dT);
	vs.xiP[1] = vs.xi[1] * (cp[1] + cpt[1] * dT);
	// d xi / dT
	vs.xiT[0] = vs.xi[0] * (-ct1[0] - ct2[0] * vs.T + cpt[0] * dP);
	vs.xiT[1] = vs.xi[1] * (-ct1[1] - ct2[1] * vs.T + cpt[1] * dP);

	// phase mass density
	vs.rho[0] = MWp[0] * vs.xi[0];
	vs.rho[1] = MWp[1] * vs.xi[1];
	// d rho / d P
	vs.rhoP[0] = MWp[0] * vs.xiP[0];
	vs.rhoP[1] = MWp[1] * vs.xiP[1];
	// d rho / d T
	vs.rhoT[0] = MWp[0] * vs.xiT[0];
	vs.rhoT[1] = MWp[1] * vs.xiT[1];

	// phase volume
	vs.vj[0] = vs.Ni[0] / vs.xi[0];
	vs.vj[1] = vs.Ni[1] / vs.xi[1];
	// total volume
	vs.Vf = vs.vj[0] + vs.vj[1];

	// phase saturation
	vs.S[0] = vs.vj[0] / vs.Vf;
	vs.S[1] = vs.vj[1] / vs.Vf;

	// d vf/ d Ni
	vs.vfi[0] = 1 / vs.xi[0];
	vs.vfi[1] = 1 / vs.xi[1];
	// d vf / d P
	vs.vfP = -(vs.vj[0] * vs.xiP[0] / vs.xi[0] + vs.vj[1] * vs.xiP[1] / vs.xi[1]);
	// d vf / d T
	vs.vfT = -(vs.vj[0] * vs.xiT[0] / vs.xi[0] + vs.vj[1] * vs.xiT[1] / vs.xi[1]);

	// Derivative of secondary vars with respect to primary vars
	vs.dXsdXp[0] = (-vs.vj[0] * vs.xiP[0] / vs.xi[0] - vs.S[0] * vs.vfP) / vs.Vf; // dSo / dP
	vs.dXsdXp[1] = (1 / vs.xi[0] - vs.S[0] * vs.vfi[0]) / vs.Vf;            // dSo / dNo
	vs.dXsdXp[2] = -vs.S[0] * vs.vfi[1] / vs.Vf;                         // dSo / dNw
	vs.dXsdXp[3] = (-vs.vj[0] * vs.xiT[0] / vs.xi[0] - vs.S[0] * vs.vfT) / vs.Vf; // dSo / dT

	vs.dXsdXp[4] = -vs.dXsdXp[0]; // dSw / dP
	vs.dXsdXp[5] = -vs.dXsdXp[1]; // dSw / dNo
	vs.dXsdXp[6] = -vs.dXsdXp[2]; // dSw / dNw
	vs.dXsdXp[7] = -vs.dXsdXp[3]; // dSw / dT

	vs.H[0] = eC.CalEnthalpy(vs.T, &vs.x[0 * 2], vs.HT[0], &vs.Hx[0 * 2]);
	vs.H[1] = eC.CalEnthalpy(vs.T, &vs.x[1 * 2], vs.HT[1], &vs.Hx[1 * 2]);

	// Internal energy per unit volume of fluid

	// Uf, d Uf / d T, d Uf / d P
	vs.Uf  = -vs.P / (GRAVITY_FACTOR * CONV6);
	vs.UfP = -1 / (GRAVITY_FACTOR * CONV6);
	vs.UfT = 0;

	for (USI j = 0; j < 2; j++) {
		// Uf
		vs.Uf += vs.S[j] * vs.xi[j] * vs.H[j];
		// dUf / dP
		vs.UfP += -(vs.vj[j] * vs.xiP[j] / vs.xi[j] + vs.S[j] * vs.vfP) / vs.Vf * vs.xi[j] * vs.H[j];
		vs.UfP += vs.xiP[j] * vs.S[j] * vs.H[j];
		// dUf / dT
		vs.UfT += -(vs.vj[j] * vs.xiT[j] / vs.xi[j] + vs.S[j] * vs.vfT) / vs.Vf * vs.xi[j] * vs.H[j];
		vs.UfT += (vs.xiT[j] * vs.H[j] + vs.HT[j] * vs.xi[j]) * vs.S[j];
	}

	// d Uf / d Ni
	vs.Ufi[0] = vs.dXsdXp[1] * vs.xi[0] * vs.H[0] + vs.dXsdXp[5] * vs.xi[1] * vs.H[1];
	vs.Ufi[1] = vs.dXsdXp[2] * vs.xi[0] * vs.H[0] + vs.dXsdXp[6] * vs.xi[1] * vs.H[1];
}


void OCPMixtureMethodK_OW01T::CalVStd(OCPMixtureVarSet& vs)
{
	vs.vj[0] = vs.Ni[0] / CalXiO(vs.P, vs.T);
	vs.vj[1] = vs.Ni[1] / CalXiW(vs.P, vs.T);
}



/////////////////////////////////////////////////////
// OCPMixtureMethodK_OGW01
/////////////////////////////////////////////////////


OCPMixtureMethodK_OGW01::OCPMixtureMethodK_OGW01(const ParamReservoir& rs_param, const USI& i,
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


void OCPMixtureMethodK_OGW01::SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const
{
	mvs.Pb = bvs.Pb[bId];
	mvs.P  = bvs.P[bId];
	mvs.T  = bvs.T[bId];
	copy(&bvs.Ni[bId * bvs.nc], &bvs.Ni[bId * bvs.nc] + bvs.nc, mvs.Ni.begin());
	copy(&bvs.S[bId * bvs.np], &bvs.S[bId * bvs.np] + bvs.np, mvs.S.begin());
}


void OCPMixtureMethodK_OGW01::CalNi(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
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


void OCPMixtureMethodK_OGW01::InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
	CalNi(Vp, vs);
	Flash(vs);
}


void OCPMixtureMethodK_OGW01::Flash(OCPMixtureVarSet& vs)
{
	FlashDer(vs);
}


void OCPMixtureMethodK_OGW01::InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
	CalNi(Vp, vs);
	FlashDer(vs);
}


void OCPMixtureMethodK_OGW01::FlashDer(OCPMixtureVarSet& vs)
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


OCP_DBL OCPMixtureMethodK_OGW01::CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt)
{
	if (pt == PhaseType::oil)       return CalXiO(P, Pb);
	else if (pt == PhaseType::gas)  return CalXiG(P);
	else if (pt == PhaseType::wat)  return CalXiW(P);
	else                            OCP_ABORT("WRONG PHASE TYPE");
}


OCP_DBL OCPMixtureMethodK_OGW01::CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt)
{
	if (pt == PhaseType::oil)       return CalRhoO(P, Pb);
	else if (pt == PhaseType::gas)  return CalRhoG(P);
	else if (pt == PhaseType::wat)  return CalRhoW(P);
	else                            OCP_ABORT("WRONG PHASE TYPE");
}


OCP_DBL OCPMixtureMethodK_OGW01::CalVmStd(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt)
{
	if      (pt == PhaseType::oil)  return (stdVo * CONV1);
    else if (pt == PhaseType::gas)  return (stdVg * 1000);
    else if (pt == PhaseType::wat)  return (stdVw * CONV1);
    else                            OCP_ABORT("Wrong Phase Type!");
}


void OCPMixtureMethodK_OGW01::CalVStd(OCPMixtureVarSet& vs)
{
	vs.vj[0] = vs.Ni[0] * stdVo * CONV1;
	vs.vj[1] = vs.Ni[1] * stdVg * 1000;
	vs.vj[2] = vs.Ni[2] * stdVw * CONV1;
}


/////////////////////////////////////////////////////
// OCPMixtureMethodK_GW01
/////////////////////////////////////////////////////


OCPMixtureMethodK_GW01::OCPMixtureMethodK_GW01(const ParamReservoir& rs_param, const USI& i, OCPMixtureVarSet& vs)
{
	vs.Init(2, 2, OCPMixtureType::BO_GW);


	PVTCO2.Setup(rs_param.PVTCO2.data[i]);
	PVTH2O.Setup(rs_param.PVTH2O.data[i]);

	garciaw.Setup(rs_param.GARCIAW);
}


void OCPMixtureMethodK_GW01::SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const
{
	mvs.P = bvs.P[bId];
	mvs.T = bvs.T[bId];
	copy(&bvs.Ni[bId * bvs.nc], &bvs.Ni[bId * bvs.nc] + bvs.nc, mvs.Ni.begin());
	copy(&bvs.S[bId * bvs.np], &bvs.S[bId * bvs.np] + bvs.np, mvs.S.begin());
}


void OCPMixtureMethodK_GW01::CalNi(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
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


void OCPMixtureMethodK_GW01::InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
	CalNi(Vp, vs);
	Flash(vs);
}


void OCPMixtureMethodK_GW01::Flash(OCPMixtureVarSet& vs)
{
	FlashDer(vs);
}


void OCPMixtureMethodK_GW01::InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
	CalNi(Vp, vs);
	FlashDer(vs);
}


void OCPMixtureMethodK_GW01::FlashDer(OCPMixtureVarSet& vs)
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

		vs.xi = vs.rho;
		vs.xiP = vs.rhoP;
		vs.xix = vs.rhox;

		vs.vj[0] = 0;
		vs.vj[1] = vs.Nt / vs.rho[1];
		vs.Vf = vs.vj[1];
		vs.S[0] = 0.0;
		vs.S[1] = 1.0;

		vs.vjP[1] = -vs.Ni[1] * vs.rhoP[1] / (vs.rho[1] * vs.rho[1]);
		vs.vfP = vs.vjP[1];

		vs.vji[1][0] = (vs.rho[1] + vs.rhox[1 * 2 + 0] * vs.Ni[1] / vs.Nt) / (vs.rho[1] * vs.rho[1]);
		vs.vji[1][1] = (vs.rho[1] - vs.rhox[1 * 2 + 0] * vs.Ni[0] / vs.Nt) / (vs.rho[1] * vs.rho[1]);
		vs.vfi[0] = vs.vji[1][0];
		vs.vfi[1] = vs.vji[1][1];

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

		vs.xi = vs.rho;
		vs.xiP = vs.rhoP;
		vs.xix = vs.rhox;

		// calculate phase mass first
		vs.nj[0] = (vs.Ni[1] - vs.Nt * (1 - xGw)) / (xWg - (1 - xGw));
		vs.nj[1] = vs.Nt - vs.nj[0];

		vs.vj[0] = vs.nj[0] / vs.rho[0];
		vs.vj[1] = vs.nj[1] / vs.rho[1];
		vs.Vf = vs.vj[0] + vs.vj[1];
		vs.S[0] = vs.vj[0] / vs.Vf;
		vs.S[1] = vs.vj[1] / vs.Vf;

		const OCP_DBL n0P = (vs.Nt * xGwP - vs.nj[0] * (xWgP + xGwP)) / (xWg - (1 - xGw));
		const OCP_DBL n1P = -n0P;
		vs.vjP[0] = (n0P * vs.rho[0] + vs.nj[0] * vs.rhoP[0]) / (vs.rho[0] * vs.rho[0]);
		vs.vjP[1] = (n1P * vs.rho[1] + vs.nj[1] * vs.rhoP[1]) / (vs.rho[1] * vs.rho[1]);
		vs.vfP = vs.vjP[0] + vs.vjP[1];

		const OCP_DBL n0N0 = -(1 - xGw) / (xWg - (1 - xGw));
		const OCP_DBL n0N1 = xGw / (xWg - (1 - xGw));
		const OCP_DBL n1N0 = 1 - n0N0;
		const OCP_DBL n1N1 = 1 - n0N1;
		vs.vji[0][0] = n0N0 / vs.rho[0];
		vs.vji[0][1] = n0N1 / vs.rho[0];
		vs.vji[1][0] = n1N0 / vs.rho[1];
		vs.vji[1][1] = n1N1 / vs.rho[1];

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


OCP_DBL OCPMixtureMethodK_GW01::CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt)
{
	if (pt == PhaseType::gas)         return CalXiG(P, T);
	else if (pt == PhaseType::wat)    return CalXiW(P, T);
	else                              OCP_ABORT("WRONG PHASE TYPE");
}


OCP_DBL OCPMixtureMethodK_GW01::CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt)
{
	if (pt == PhaseType::gas)         return CalRhoG(P, T);
	else if (pt == PhaseType::wat)    return CalRhoW(P, T);
	else                              OCP_ABORT("WRONG PHASE TYPE");
}


OCP_DBL OCPMixtureMethodK_GW01::CalRhoW(const OCP_DBL& P, const OCP_DBL& T) const
{
	OCP_DBL rhow, xGw, dummy;
	PVTH2O.CalRhoMuSol(P, T, rhow, dummy, xGw);
	garciaw.CalRho(T, xGw, rhow);
	return rhow;
}


OCP_DBL OCPMixtureMethodK_GW01::CalVmStd(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt)
{
}


void OCPMixtureMethodK_GW01::CalVStd(OCPMixtureVarSet& vs)
{
}






/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/02/2023      Create file                          */
/*----------------------------------------------------------------------------*/