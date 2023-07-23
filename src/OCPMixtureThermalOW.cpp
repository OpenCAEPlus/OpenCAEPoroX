/*! \file    OCPMixtureBlkOilOW.cpp
 *  \brief   OCPMixtureBlkOilOW class declaration
 *  \author  Shizhe Li
 *  \date    Jul/20/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPMixtureThermalOW.hpp"


/////////////////////////////////////////////////////
// OCPMixtureUnitThermalOWMethod01
/////////////////////////////////////////////////////


OCPMixtureUnitThermalOWMethod01::OCPMixtureUnitThermalOWMethod01(const ComponentParam& param, const USI& tarId, OCPMixtureVarSet& vs)
{
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

    vs.xij[0 * 2 + 0] = 1;
    vs.xij[0 * 2 + 1] = 0;
    vs.xij[1 * 2 + 0] = 0;
    vs.xij[1 * 2 + 1] = 1;

    // d mu / dP
    fill(vs.rhox.begin(), vs.rhox.end(), 0.0);
    fill(vs.xix.begin(), vs.xix.end(), 0.0);
    fill(vs.muP.begin(), vs.muP.end(), 0.0);
    fill(vs.mux.begin(), vs.mux.end(), 0.0);
}


OCP_DBL OCPMixtureUnitThermalOWMethod01::CalRhoO(const OCP_DBL& P, const OCP_DBL& T)
{
    return MWp[0] * CalXiO(P, T);
}


OCP_DBL OCPMixtureUnitThermalOWMethod01::CalXiO(const OCP_DBL& P, const OCP_DBL& T)
{
    const OCP_DBL dP = P - Pref;
    const OCP_DBL dT = T - Tref + CONV5;
    return xi_ref[0] * exp(cp[0] * dP - ct1[0] * dT - ct2[0] * (pow(T, 2) - pow(Tref, 2)) / 2 + cpt[0] * dP * dT);
}


OCP_DBL OCPMixtureUnitThermalOWMethod01::CalRhoW(const OCP_DBL& P, const OCP_DBL& T)
{
    return MWp[1] * CalXiW(P, T);
}


OCP_DBL OCPMixtureUnitThermalOWMethod01::CalXiW(const OCP_DBL& P, const OCP_DBL& T)
{
    const OCP_DBL dP = P - Pref;
    const OCP_DBL dT = T - Tref + CONV5;
	return xi_ref[1] * exp(cp[1] * dP - ct1[1] * dT - ct2[1] * (pow(T, 2) - pow(Tref, 2)) / 2 + cpt[1] * dP * dT);
}


void OCPMixtureUnitThermalOWMethod01::InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
    vs.Ni[0] = Vp * vs.S[0] * CalXiO(vs.P, vs.T - CONV5);
    vs.Ni[1] = Vp * vs.S[1] * CalXiW(vs.P, vs.T - CONV5);
    Flash(vs);
}


void OCPMixtureUnitThermalOWMethod01::Flash(OCPMixtureVarSet& vs)
{
    FlashDer(vs);
}


void OCPMixtureUnitThermalOWMethod01::InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
    vs.Ni[0] = Vp * vs.S[0] * CalXiO(vs.P, vs.T - CONV5);
    vs.Ni[1] = Vp * vs.S[1] * CalXiW(vs.P, vs.T - CONV5);
    FlashDer(vs);
}


void OCPMixtureUnitThermalOWMethod01::FlashDer(OCPMixtureVarSet& vs)
{
    // Assign value
    const OCP_DBL dP = vs.P - Pref;
    const OCP_DBL dT = vs.T - Tref;

    vs.Nt = vs.Ni[0] + vs.Ni[1];

    // phase viscosity
    vs.mu[0] = vC.CalViscosity(vs.P, vs.T, &vs.xij[0 * 2], vs.muP[0], vs.muT[0], &vs.mux[0 * 2]);
    vs.mu[1] = vC.CalViscosity(vs.P, vs.T, &vs.xij[1 * 2], vs.muP[1], vs.muT[1], &vs.mux[1 * 2]);

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
    vs.vf    = vs.vj[0] + vs.vj[1];

    // phase saturation
    vs.S[0] = vs.vj[0] / vs.vf;
    vs.S[1] = vs.vj[1] / vs.vf;

    // d vf/ d Ni
    vs.vfi[0] = 1 / vs.xi[0];
    vs.vfi[1] = 1 / vs.xi[1];
    // d vf / d P
    vs.vfP = -(vs.vj[0] * vs.xiP[0] / vs.xi[0] + vs.vj[1] * vs.xiP[1] / vs.xi[1]);
    // d vf / d T
    vs.vfT = -(vs.vj[0] * vs.xiT[0] / vs.xi[0] + vs.vj[1] * vs.xiT[1] / vs.xi[1]);

    // Derivative of secondary vars with respect to primary vars
    vs.dXsdXp[0] = (-vs.vj[0] * vs.xiP[0] / vs.xi[0] - vs.S[0] * vs.vfP) / vs.vf; // dSo / dP
    vs.dXsdXp[1] = (1 / vs.xi[0] - vs.S[0] * vs.vfi[0]) / vs.vf;            // dSo / dNo
    vs.dXsdXp[2] = -vs.S[0] * vs.vfi[1] / vs.vf;                         // dSo / dNw
    vs.dXsdXp[3] = (-vs.vj[0] * vs.xiT[0] / vs.xi[0] - vs.S[0] * vs.vfT) / vs.vf; // dSo / dT

    vs.dXsdXp[4] = -vs.dXsdXp[0]; // dSw / dP
    vs.dXsdXp[5] = -vs.dXsdXp[1]; // dSw / dNo
    vs.dXsdXp[6] = -vs.dXsdXp[2]; // dSw / dNw
    vs.dXsdXp[7] = -vs.dXsdXp[3]; // dSw / dT

    vs.H[0] = eC.CalEnthalpy(vs.T, &vs.xij[0 * 2], vs.HT[0], &vs.Hx[0 * 2]);
    vs.H[1] = eC.CalEnthalpy(vs.T, &vs.xij[1 * 2], vs.HT[1], &vs.Hx[1 * 2]);

    // Internal energy per unit volume of fluid

    // Uf, d Uf / d T, d Uf / d P
    vs.Uf  = -vs.P / (GRAVITY_FACTOR * CONV6);
    vs.UfP = -1 / (GRAVITY_FACTOR * CONV6);
    vs.UfT = 0;

    for (USI j = 0; j < 2; j++) {
        // Uf
        vs.Uf += vs.S[j] * vs.xi[j] * vs.H[j];
        // dUf / dP
        vs.UfP += -(vs.vj[j] * vs.xiP[j] / vs.xi[j] + vs.S[j] * vs.vfP) / vs.vf * vs.xi[j] * vs.H[j];
        vs.UfP += vs.xiP[j] * vs.S[j] * vs.H[j];
        // dUf / dT
        vs.UfT += -(vs.vj[j] * vs.xiT[j] / vs.xi[j] + vs.S[j] * vs.vfT) / vs.vf * vs.xi[j] * vs.H[j];
        vs.UfT += (vs.xiT[j] * vs.H[j] + vs.HT[j] * vs.xi[j]) * vs.S[j];
    }

    // d Uf / d Ni
    vs.Ufi[0] = vs.dXsdXp[1] * vs.xi[0] * vs.H[0] + vs.dXsdXp[5] * vs.xi[1] * vs.H[1];
    vs.Ufi[1] = vs.dXsdXp[2] * vs.xi[0] * vs.H[0] + vs.dXsdXp[6] * vs.xi[1] * vs.H[1];
}


/////////////////////////////////////////////////////
// OCPMixtureUnitThermalOW 
/////////////////////////////////////////////////////


void OCPMixtureUnitThermalOW::Setup(const ParamReservoir& rs_param, const USI& i)
{
    vs.Init(2, 2, OCP_TRUE);
    pmMethod = new OCPMixtureUnitThermalOWMethod01(rs_param.comsParam, i, vs);
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/20/2023      Create file                          */
/*----------------------------------------------------------------------------*/