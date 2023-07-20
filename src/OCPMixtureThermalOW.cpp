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
// OCPMixtureThermalOWMethod01
/////////////////////////////////////////////////////


OCPMixtureThermalOWMethod01::OCPMixtureThermalOWMethod01(const ParamReservoir& param, const USI& tarId, OCPMixtureVarSet& vs)
{
    if (param.comsParam.molden.activity)
        xi_ref = param.comsParam.molden.data[tarId];
    else
        OCP_ABORT("ACF hasn't been input!");
    if (param.comsParam.avisc.activity || param.comsParam.bvisc.activity) {
        if (param.comsParam.avisc.activity)
            avisc = param.comsParam.avisc.data[tarId];
        else
            avisc.resize(2, 0);
        if (param.comsParam.bvisc.activity)
            bvisc = param.comsParam.bvisc.data[tarId];
        else
            bvisc.resize(2, 0);
        useViscTab = OCP_FALSE;
    }
    else {
        if (param.comsParam.viscTab.data.size() <= tarId) {
            OCP_ABORT("VISCTAB hasn't been input for " + to_string(tarId + 1) +
                " th Region!");
        }
        useViscTab = OCP_TRUE;
        visc.Setup(param.comsParam.viscTab.data[tarId]);
        // unit convert: F -> R
        for (auto& v : visc.GetCol(0)) {
            v += CONV5;
        }
    }

    if (param.comsParam.cp.activity)
        cp = param.comsParam.cp.data[tarId];
    else
        cp.resize(2, 0);

    if (param.comsParam.ct1.activity)
        ct1 = param.comsParam.ct1.data[tarId];
    else
        ct1.resize(2, 0);

    if (param.comsParam.ct2.activity)
        ct2 = param.comsParam.ct2.data[tarId];
    else
        ct2.resize(2, 0);

    if (param.comsParam.cpt.activity)
        cpt = param.comsParam.cpt.data[tarId];
    else
        cpt.resize(2, 0);

    if (param.comsParam.MW.activity)
        MWc = param.comsParam.MW.data[tarId];
    else
        OCP_ABORT("MW hasn't been input!");

    Tref = param.comsParam.Tref[tarId] + CONV5;
    Pref = param.comsParam.Pref[tarId];
  
    data.resize(3, 0);
    cdata.resize(3, 0);

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


OCP_DBL OCPMixtureThermalOWMethod01::CalRhoO(const OCP_DBL& P, const OCP_DBL& T)
{
    return MWp[0] * CalXiO(P, T);
}


OCP_DBL OCPMixtureThermalOWMethod01::CalXiO(const OCP_DBL& P, const OCP_DBL& T)
{
    const OCP_DBL dP = P - Pref;
    const OCP_DBL dT = T - Tref + CONV5;
    return xi_ref[0] * exp(cp[0] * dP - ct1[0] * dT - ct2[0] * pow(dT, 2) / 2 + cpt[0] * dP * dT);
}


OCP_DBL OCPMixtureThermalOWMethod01::CalRhoW(const OCP_DBL& P, const OCP_DBL& T)
{
    return MWp[1] * CalXiW(P, T);
}


OCP_DBL OCPMixtureThermalOWMethod01::CalXiW(const OCP_DBL& P, const OCP_DBL& T)
{
    const OCP_DBL dP = P - Pref;
    const OCP_DBL dT = T - Tref + CONV5;
	return xi_ref[1] * exp(cp[1] * dP - ct1[1] * dT - ct2[1] * pow(dT, 2) / 2 + cpt[1] * dP * dT);
}


void OCPMixtureThermalOWMethod01::InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{

}


void OCPMixtureThermalOWMethod01::Flash(OCPMixtureVarSet& vs)
{

}


void OCPMixtureThermalOWMethod01::InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{

}


void OCPMixtureThermalOWMethod01::FlashDer(OCPMixtureVarSet& vs)
{
    // Assign value
    const OCP_DBL dP = vs.P - Pref;
    const OCP_DBL dT = vs.T - Tref;

    vs.Nt = vs.Ni[0] + vs.Ni[1];

    // phase viscosity
    if (useViscTab) {
        visc.Eval_All(0, vs.T, data, cdata);
        vs.mu[0] = data[1];
        vs.mu[1] = data[2];
        // d mu / dT
        vs.muT[0] = cdata[1];
        vs.muT[1] = cdata[2];

    } else {
        vs.mu[0] = avisc[0] * exp(bvisc[0] / vs.T);
        vs.mu[1] = avisc[1] * exp(bvisc[1] / vs.T);
        // d mu / dT
        vs.muT[0] = vs.mu[0] * (-bvisc[0] / (vs.T * vs.T));
        vs.muT[1] = vs.mu[1] * (-bvisc[1] / (vs.T * vs.T));
    }

    // phase molar density
    vs.xi[0] = xi_ref[0] *
            exp(cp[0] * dP - ct1[0] * dT - ct2[0] * pow(dT, 2) / 2 + cpt[0] * dP * dT);
    vs.xi[1] = xi_ref[1] *
            exp(cp[1] * dP - ct1[1] * dT - ct2[1] * pow(dT, 2) / 2 + cpt[1] * dP * dT);
    // d xi / dP
    vs.xiP[0] = vs.xi[0] * (cp[0] + cpt[0] * dT);
    vs.xiP[1] = vs.xi[1] * (cp[1] + cpt[1] * dT);
    // d xi / dT
    vs.xiT[0] = vs.xi[0] * (-ct1[0] - ct2[0] * dT + cpt[0] * dP);
    vs.xiT[1] = vs.xi[1] * (-ct1[1] - ct2[1] * dT + cpt[1] * dP);

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
    vs.vf = vs.vj[0] + vs.vj[1];

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


    // Calculate Enthalpy
    //fill(vs.H.begin(), vs.H.end(), 0.0);
    //fill(vs.HT.begin(), vs.HT.end(), 0.0);

    //const OCP_DBL dT = vs.T - Tref;

    //if (liquid_based || simple_hvap) {
    //    for (USI j = 0; j < 2; j++) {
    //        for (USI i = 0; i < 2; i++) {

    //            Hx[j * 2 + i] =
    //                (cpl1[i] * dT + 1.0 / 2 * cpl2[i] * pow(dT, 2) +
    //                 1.0 / 3 * cpl3[i] * pow(dT, 3) + 1.0 / 4 * cpl4[i] * pow(dT, 4));

    //            H[j] += xij[j * 2 + i] * Hx[j * 2 + i];

    //            HT[j] +=
    //                xij[j * 2 + i] * (cpl1[i] + cpl2[i] * dT +
    //                                       cpl3[i] * pow(dT, 2) + cpl4[i] * pow(dT, 3));
    //        }
    //    }
    //} else if (gas_based) {
    //    for (USI j = 0; j < 2; j++) {
    //        for (USI i = 0; i < 2; i++) {
    //            Hx[j * 2 + i] = cpg1[i] * dT + 1.0 / 2 * cpg2[i] * pow(dT, 2) +
    //                                 1.0 / 3 * cpg3[i] * pow(dT, 3) +
    //                                 1.0 / 4 * cpg4[i] * pow(dT, 4);
    //            H[j] += xij[j * 2 + i] * Hx[j * 2 + i];
    //            HT[j] +=
    //                xij[j * 2 + i] * (cpg1[i] + cpg2[i] * dT +
    //                                       cpg3[i] * pow(dT, 2) + cpg4[i] * pow(dT, 3));

    //            if (T < Tcrit[i]) {
    //                Hx[j * 2 + i] -= hvr[i] * pow((Tcrit[i] - T), ev[i]);

    //                H[j] -= xij[j * 2 + i] * hvr[i] * pow((Tcrit[i] - T), ev[i]);

    //                HT[j] += xij[j * 2 + i] * hvr[i] * ev[i] *
    //                         pow((Tcrit[i] - T), ev[i] - 1);
    //            }
    //        }
    //    }
    //} else {
    //    OCP_ABORT("WRONG Type !");
    //}

    //// Internal energy per unit volume of fluid

    //// Uf, d Uf / d T, d Uf / d P
    //Uf  = -P / (GRAVITY_FACTOR * CONV6);
    //UfP = -1 / (GRAVITY_FACTOR * CONV6);
    //UfT = 0;

    //for (USI j = 0; j < numPhase; j++) {
    //    // Uf
    //    Uf += S[j] * xi[j] * H[j];
    //    // dUf / dP
    //    UfP += -(vj[j] * xiP[j] / xi[j] + S[j] * vfP) / vf * xi[j] * H[j];
    //    UfP += xiP[j] * S[j] * H[j];
    //    // dUf / dT
    //    UfT += -(vj[j] * xiT[j] / xi[j] + S[j] * vfT) / vf * xi[j] * H[j];
    //    UfT += (xiT[j] * H[j] + HT[j] * xi[j]) * S[j];
    //}

    //// d Uf / d Ni
    //Ufi[0] = dXsdXp[1] * xi[0] * H[0] + dXsdXp[5] * xi[1] * H[1];
    //Ufi[1] = dXsdXp[2] * xi[0] * H[0] + dXsdXp[6] * xi[1] * H[1];
}


/////////////////////////////////////////////////////
// OCPMixtureThermalOW 
/////////////////////////////////////////////////////



void OCPMixtureThermalOW::Setup(const ParamReservoir& rs_param, const USI& i)
{
    vs.Init(2, 2, OCP_TRUE);
    pmMethod = new OCPMixtureThermalOWMethod01(rs_param, i, vs);
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/20/2023      Create file                          */
/*----------------------------------------------------------------------------*/