/*! \file    MixtureBO3_ODGW.cpp
 *  \brief   Used for the condition where oil, gas, disolve gas, water exist.
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "MixtureBO.hpp"

///////////////////////////////////////////////
// BOMixture_ODGW
///////////////////////////////////////////////

BOMixture_ODGW::BOMixture_ODGW(const ParamReservoir& rs_param, const USI& i)
{
    mixtureType = BLKOIL_ODGW;
    numPhase    = 3;
    numCom      = 3;
    BOMixtureInit(rs_param);

    PVTW.Setup(rs_param.PVTW_T.data[i], std_RhoW);  
    PVCO.Setup(rs_param.PVCO_T.data[i], std_RhoO, std_RhoG);
    PVDG.Setup(rs_param.PVDG_T.data[i], std_RhoG);
}

void BOMixture_ODGW::Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{
    FlashIMPEC(Pin, Tin, Niin, 0, 0, 0);
}


void BOMixture_ODGW::CalNi(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL* Sjin, const OCP_DBL& Vpore)
{
        // Calculate Ni only frist, then call FlashDer()
	for (USI i = 0; i < 3; i++) Ni[i] = 0;

	P    = Pin;
	S[1] = Sjin[1]; // gas saturation
	S[2] = Sjin[2]; // water saturation

	// Water Property
	Ni[2] = Vpore * S[2] * PVTW.CalXiW(P);

	USI phasecae;
	if (1 - S[1] - S[2] < TINY) {
		if (S[1] < TINY)
			phasecae = PHASE_W; // case 1 : water, no oil, no gas
		else
			phasecae = PHASE_GW; // case 2 : water, gas, no oil
	}
	else if (S[1] < TINY)
		phasecae = PHASE_OW; // case 3 : water, oil, no gas
	else
		phasecae = PHASE_ODGW; // case 4 : water, oil, gas

	switch (phasecae) {
	    case PHASE_W:
	    {
	    	// water
	    	break;
	    }
	    case PHASE_GW:
	    {
	    	// water, gas
	    	Ni[1] = Vpore * S[1] * PVDG.CalXiG(P);
	    	break;
	    }
	    case PHASE_OW:
	    {
	    	// water, oil, unsaturated
	    	// Ni[0] and Ni[1]
	    	const OCP_DBL Pbb = Pbbin;
	    	const OCP_DBL rs = PVCO.CalRs(Pbb);
	    	Ni[0] = Vpore * (1 - S[1] - S[2]) * PVCO.CalXiO(P, Pbb) / (1 + rs);
	    	Ni[1] = Ni[0] * rs;

	    	break;
	    }
	    case PHASE_ODGW:
	    {
		    const OCP_DBL rs = PVCO.CalRs(P);
		    Ni[0] = Vpore * (1 - S[1] - S[2]) * PVCO.CalXiO(P, P) / (1 + rs);
		    Ni[1] = Vpore * S[1] * PVDG.CalXiG(P) + Ni[0] * rs;
		    break;
	    }
	}
}


void BOMixture_ODGW::InitFlashIMPEC(const OCP_DBL& Pin,
                                    const OCP_DBL& Pbbin,
                                    const OCP_DBL& Tin,
                                    const OCP_DBL* Sjin,
                                    const OCP_DBL& Vpore,
                                    const OCP_DBL* Ziin,
                                    const OCP_USI& bId)
{
    CalNi(Pin, Pbbin, Sjin, Vpore);
    FlashFIM(Pin, Tin, &Ni[0], 0, 0, 0, 0);
}

void BOMixture_ODGW::InitFlashFIM(const OCP_DBL& Pin,
                                  const OCP_DBL& Pbbin,
                                  const OCP_DBL& Tin,
                                  const OCP_DBL* Sjin,
                                  const OCP_DBL& Vpore,
                                  const OCP_DBL* Ziin,
                                  const OCP_USI& bId)
{
    CalNi(Pin, Pbbin, Sjin, Vpore);
    FlashFIM(Pin, Tin, &Ni[0], 0, 0, 0, 0);
}

void BOMixture_ODGW::FlashIMPEC(const OCP_DBL& Pin,
                                const OCP_DBL& Tin,
                                const OCP_DBL* Niin,
                                const USI&     lastNP,
                                const OCP_DBL* xijin,
                                const OCP_USI& bId)
{
    for (USI j = 0; j < 3; j++) phaseExist[j] = OCP_FALSE;
    fill(xij.begin(), xij.end(), 0.0);

    P  = Pin;
    Nt = 0;
    for (USI i = 0; i < numCom; i++) {
        Ni[i] = Niin[i];
        Nt += Ni[i];
    }

    // Water property
    PVTW.CalRhoXiMuDer(P, rho[2], xi[2], mu[2], rhoP[2], xiP[2], muP[2]);

    USI     phasecase;
    const OCP_DBL Rs_sat = PVCO.CalRs(P);

    if (Ni[0] < Nt * TINY) {
        if (Ni[1] <= Ni[0] * Rs_sat)
            phasecase = PHASE_W; // water, no oil, no gas
        else
            phasecase = PHASE_GW; // water, gas, no oil
    } else if (Ni[1] <= Ni[0] * Rs_sat)
        phasecase = PHASE_OW; // water, oil, no gas
    else
        phasecase = PHASE_ODGW; // water, oil ,gas

    switch (phasecase) {
        case PHASE_W:
            {
                // water
                phaseExist[2]  = OCP_TRUE;
                S[0]           = 0;
                S[1]           = 0;
                S[2]           = 1;
                xij[2 * 3 + 2] = 1;

                // total
                vj[0]  = 0;
                vj[1]  = 0;
                vj[2]  = Ni[2] / xi[2];
                vf     = vj[2];
                vfP    = -Ni[2] * xiP[2] / (xi[2] * xi[2]);
                vfi[0] = 0;
                vfi[1] = 0;
                vfi[2] = 1 / xi[2];

                break;
            }
        case PHASE_GW:
            {
                // water, gas
                phaseExist[1]  = OCP_TRUE;
                phaseExist[2]  = OCP_TRUE;
                xij[1 * 3 + 1] = 1;
                xij[2 * 3 + 2] = 1;

                // hypothetical Oil property
                OCP_DBL rs, rsP;
                PVCO.CalRhoXiMuRsDer(P, rho[0], xi[0], mu[0], rs, rhoP[0], xiP[0], muP[0], rsP);

                // gas property
                PVDG.CalRhoXiMuDer(P, rho[1], xi[1], mu[1], rhoP[1], xiP[1], muP[1]);

                // total
                vj[0] = 0;
                vj[1] = (Ni[1] - rs * Ni[0]) / xi[1];
                vj[2] = Ni[2] / xi[2];

                vf     = vj[1] + vj[2];
                S[0]   = 0;
                S[1]   = vj[1] / vf;
                S[2]   = vj[2] / vf;

                vfP = Ni[0] * (rsP * xi[0] - (1 + rs) * xiP[0]) / (xi[0] * xi[0]) +
                    (-rsP * Ni[0] * xi[1] - (Ni[1] - rs * Ni[0]) * xiP[1]) / (xi[1] * xi[1])
                    - Ni[2] * xiP[2] / (xi[2] * xi[2]);
                vfi[0] = (1 + rs) / xi[0] - rs / xi[1];
                vfi[1] = 1 / xi[1];
                vfi[2] = 1 / xi[2];

                break;
            }
        case PHASE_OW:
            {
                // water, oil
                phaseExist[0]  = OCP_TRUE;
                phaseExist[2]  = OCP_TRUE;
                xij[0 * 3 + 0] = Ni[0] / (Ni[0] + Ni[1]);
                xij[0 * 3 + 1] = 1 - xij[0 * 3 + 0];
                xij[1 * 3 + 1] = 1;
                xij[2 * 3 + 2] = 1;

                // oil property
                const OCP_DBL rs = Ni[1] / Ni[0];
                OCP_DBL rhooRs, xioRs, muoRs;
                PVCO.CalRhoXiMuDer(rs, P, rho[0], xi[0], mu[0], rhoP[0], xiP[0], muP[0], rhooRs, xioRs, muoRs);

                // total
                vj[0]  = Ni[0] * (1 + rs) / xi[0];
                vj[1]  = 0;
                vj[2]  = Ni[2] / xi[2];
                vf     = vj[0] + vj[2];
                S[0]   = vj[0] / vf;
                S[1]   = 0;
                S[2]   = vj[2] / vf;

                vfP    = -Ni[0] * (1 + rs) * xiP[0] / (xi[0] * xi[0]) - Ni[2] * xiP[2] / (xi[2] * xi[2]);
                vfi[0] = (1 + rs) / xi[0] - (xi[0] - (1 + rs) * xioRs) / (xi[0] * xi[0]) * rs;
                vfi[1] = (xi[0] - (1 + rs) * xioRs) / (xi[0] * xi[0]);
                vfi[2] = 1 / xi[2];

                break;
            }
        case PHASE_ODGW:
            {
                phaseExist[0] = OCP_TRUE;
                phaseExist[1] = OCP_TRUE;
                phaseExist[2] = OCP_TRUE;

                // oil property
                OCP_DBL rs, rsP;
                PVCO.CalRhoXiMuRsDer(P, rho[0], xi[0], mu[0], rs, rhoP[0], xiP[0], muP[0], rsP);

                // gas property
                PVDG.CalRhoXiMuDer(P, rho[1], xi[1], mu[1], rhoP[1], xiP[1], muP[1]);

                // total
                xij[0 * 3 + 0] = 1 / (1 + rs);
                xij[0 * 3 + 1] = 1 - xij[0 * 3 + 0];
                xij[1 * 3 + 1] = 1;
                xij[2 * 3 + 2] = 1;

                vj[0] = Ni[0] * (1 + rs) / xi[0];
                vj[1] = (Ni[1] - rs * Ni[0]) / xi[1];
                vj[2] = Ni[2] / xi[2];
                vf    = vj[0] + vj[1] + vj[2];
                S[0]  = vj[0] / vf;
                S[1]  = vj[1] / vf;
                S[2]  = vj[2] / vf;

                vfP = Ni[0] * (rsP * xi[0] - (1 + rs) * xiP[0]) / (xi[0] * xi[0]) +
                    (-rsP * Ni[0] * xi[1] - (Ni[1] - rs * Ni[0]) * xiP[1]) / (xi[1] * xi[1]) +
                    -Ni[2] * xiP[2] / (xi[2] * xi[2]);
                vfi[0] = (1 + rs) / xi[0] - rs / xi[1];
                vfi[1] = 1 / xi[1];
                vfi[2] = 1 / xi[2];

                break;
            }
    }
}

void BOMixture_ODGW::FlashFIM(const OCP_DBL& Pin,
                              const OCP_DBL& Tin,
                              const OCP_DBL* Niin,
                              const OCP_DBL* Sjin,
                              const USI&     lastNP,
                              const OCP_DBL* xijin,
                              const OCP_USI& bId)
{
    for (USI j = 0; j < 3; j++) {
        phaseExist[j] = OCP_FALSE;
        rhoP[j]       = 0;
        xiP[j]        = 0;
        muP[j]        = 0;
    }
    fill(xij.begin(), xij.end(), 0.0);
    fill(rhox.begin(), rhox.end(), 0.0);
    fill(xix.begin(), xix.end(), 0.0);
    fill(mux.begin(), mux.end(), 0.0);
    fill(dXsdXp.begin(), dXsdXp.end(), 0.0);

    P  = Pin;
    Nt = 0;
    for (USI i = 0; i < numCom; i++) {
        Ni[i] = Niin[i];
        Nt += Ni[i];
    }

    // Water property
    PVTW.CalRhoXiMuDer(P, rho[2], xi[2], mu[2], rhoP[2], xiP[2], muP[2]);

    USI     phasecase;
    const OCP_DBL Rs_sat = PVCO.CalRs(P);

    if (Ni[0] < Nt * TINY) {
        if (Ni[1] <= Ni[0] * Rs_sat)
            phasecase = PHASE_W; // water, no oil, no gas
        else
            phasecase = PHASE_GW; // water, gas, no oil
    }
    // Beacareful when NO = NG = 0, if it's possible.
    else if (Ni[1] <= Ni[0] * Rs_sat)
        phasecase = PHASE_OW; // water, oil, no gas
    else
        phasecase = PHASE_ODGW; // water, oil ,gas

    switch (phasecase) {
        case PHASE_W:
            {
                // water
                phaseExist[2]  = OCP_TRUE;
                S[0]           = 0;
                S[1]           = 0;
                S[2]           = 1;
                xij[2 * 3 + 2] = 1;

                // unsaturated oil, water
                xi[1] = PVDG.CalXiG(P);
                // total
                vj[0] = 0;
                vj[1] = 0;
                vj[2] = Ni[2] / xi[2];
                vf    = vj[2];

                vfP = -Ni[2] * xiP[2] / (xi[2] * xi[2]);
                vfi[0] = 0;
                vfi[1] = 0;
                vfi[2] = 1 / xi[2];

                // dXsdXp[0] = 0;      // dSo / dP
                // dXsdXp[1] = 0;      // dSo / dNo
                // dXsdXp[2] = 0;      // dSo / dNg
                // dXsdXp[3] = 0;      // dSo / dNw

                // dXsdXp[1 * 4 + 0] = 0;             // dSg / dP
                // dXsdXp[1 * 4 + 1] = 0;             // dSg / dNo
                dXsdXp[1 * 4 + 2] = 1 / xi[1] / vf;   // dSg / dNg
                // dXsdXp[1 * 4 + 3] = 0;             // dSg / dNw

                dXsdXp[2 * 4 + 0] = (-Ni[2] * xiP[2] / (xi[2] * xi[2]) - S[2] * vfP) / vf; // dSw / dP
                dXsdXp[2 * 4 + 1] = -S[2] * vfi[0] / vf;               // dSw / dNo
                dXsdXp[2 * 4 + 2] = -S[2] * vfi[1] / vf;               // dSw / dNg
                dXsdXp[2 * 4 + 3] = (1 / xi[2] - S[2] * vfi[2]) / vf;  // dSw / dNw

                break;
            }
        case PHASE_GW:
            {
                // water, gas
                phaseExist[1]  = OCP_TRUE;
                phaseExist[2]  = OCP_TRUE;
                xij[1 * 3 + 1] = 1;
                xij[2 * 3 + 2] = 1;

                // hypothetical Oil property
                OCP_DBL rs, rsP;
                PVCO.CalRhoXiMuRsDer(P, rho[0], xi[0], mu[0], rs, rhoP[0], xiP[0], muP[0], rsP);

                // gas property
                PVDG.CalRhoXiMuDer(P, rho[1], xi[1], mu[1], rhoP[1], xiP[1], muP[1]);

                // total
                vj[0] = 0;
                vj[1] = (Ni[1] - rs * Ni[0]) / xi[1];
                vj[2] = Ni[2] / xi[2];

                vf   = vj[1] + vj[2];
                S[0] = 0;
                S[1] = vj[1] / vf;
                S[2] = vj[2] / vf;

                vfP    = Ni[0] * (rsP * xi[0] - (1 + rs) * xiP[0]) / (xi[0] * xi[0]) +
                    (-rsP * Ni[0] * xi[1] - (Ni[1] - rs * Ni[0]) * xiP[1]) / (xi[1] * xi[1])
                    -Ni[2] * xiP[2] / (xi[2] * xi[2]);
                vfi[0] = (1 + rs) / xi[0] - rs / xi[1];
                vfi[1] = 1 / xi[1];
                vfi[2] = 1 / xi[2];

                // discard small value
                // dXsdXp[0] = (Ni[0] * (rsP * xi[0] - (1 + rs) * xiP[0]) / (xi[0] * xi[0])) / vf; // dSo / dP
                dXsdXp[1] = ((1 + rs) / xi[0]) / vf;                  // dSo / dNo
                dXsdXp[2] = 0;                                        // dSo / dNg
                dXsdXp[3] = 0;                                        // dSo / dNw
               
                dXsdXp[1 * 4 + 0] = ((-rsP * Ni[0] * xi[1] - (Ni[1] - rs * Ni[0]) * xiP[1]) / (xi[1] * xi[1]) - S[1] * vfP) / vf; // dSg / dP
                dXsdXp[1 * 4 + 1] = (-rs / xi[1] - S[1] * vfi[0]) / vf;    // dSg / dNo
                dXsdXp[1 * 4 + 2] = (1 / xi[1] - S[1] * vfi[1]) / vf;      // dSg / dNg
                dXsdXp[1 * 4 + 3] = -S[1] * vfi[2] / vf;                   // dSg / dNw

                dXsdXp[2 * 4 + 0] = (-Ni[2] * xiP[2] / (xi[2] * xi[2]) - S[2] * vfP) / vf; // dSw / dP
                dXsdXp[2 * 4 + 1] = -S[2] * vfi[0] / vf;               // dSw / dNo
                dXsdXp[2 * 4 + 2] = -S[2] * vfi[1] / vf;               // dSw / dNg
                dXsdXp[2 * 4 + 3] = (1 / xi[2] - S[2] * vfi[2]) / vf;  // dSw / dNw

                //dXsdXp[3 * 4 + 0] = -rsP / ((1 + rs) * (1 + rs)); // d Xoo / dP
                //dXsdXp[4 * 4 + 0] = -dXsdXp[3 * 4 + 0];           // d Xgo / dP
                break;
            }
        case PHASE_OW:
            {
                // unsaturated oil, water
                phaseExist[0]  = OCP_TRUE;
                phaseExist[2]  = OCP_TRUE;
                xij[0 * 3 + 0] = Ni[0] / (Ni[0] + Ni[1]);
                xij[0 * 3 + 1] = 1 - xij[0 * 3 + 0];
                xij[1 * 3 + 1] = 1;
                xij[2 * 3 + 2] = 1;

                // oil property
                const OCP_DBL rs = Ni[1] / Ni[0];
                OCP_DBL rhooRs, xioRs, muoRs;
                PVCO.CalRhoXiMuDer(rs, P, rho[0], xi[0], mu[0], rhoP[0], xiP[0], muP[0], rhooRs, xioRs, muoRs);

                // total
                vj[0]  = Ni[0] * (1 + rs) / xi[0];
                vj[1]  = 0;
                vj[2]  = Ni[2] / xi[2];
                vf     = vj[0] + vj[2];
                S[0]   = vj[0] / vf;
                S[1]   = 0;
                S[2]   = vj[2] / vf;
                
                vfP = -Ni[0] * (1 + rs) * xiP[0] / (xi[0] * xi[0]) - Ni[2] * xiP[2] / (xi[2] * xi[2]);

                vfi[0] = (1 + rs) / xi[0] - (xi[0] - (1 + rs) * xioRs) / (xi[0] * xi[0]) * rs;
                vfi[1] = (xi[0] - (1 + rs) * xioRs) / (xi[0] * xi[0]);
                vfi[2] = 1 / xi[2];

                dXsdXp[0] =
                    (-Ni[0] * (1 + rs) * xiP[0] / (xi[0] * xi[0]) - S[0] * vfP) / vf; // dSo / dP
                dXsdXp[1] =
                    ((1 + rs) / xi[0] - (xi[0] - (1 + rs) * xioRs) / (xi[0] * xi[0]) * rs - S[0] * vfi[0]) /
                    vf;                                             // dSo / dNo
                dXsdXp[2] = ((xi[0] - (1 + rs) * xioRs) / (xi[0] * xi[0]) - S[0] * vfi[1]) / vf; // dSo / dNg
                dXsdXp[3] = -S[0] / vf * vfi[2];                    // dSo / dNw

                dXsdXp[2 * 4 + 0] = (-Ni[2] * xiP[2] / (xi[2] * xi[2]) - S[2] * vfP) / vf; // dSw / dP
                dXsdXp[2 * 4 + 1] = -S[2] * vfi[0] / vf;               // dSw / dNo
                dXsdXp[2 * 4 + 2] = -S[2] * vfi[1] / vf;               // dSw / dNg
                dXsdXp[2 * 4 + 3] = (1 / xi[2] - S[2] * vfi[2]) / vf; // dSw / dNw

                dXsdXp[3 * 4 + 1] = Ni[1] / pow((Ni[0] + Ni[1]), 2);  // d Xoo / d No
                dXsdXp[3 * 4 + 2] = -Ni[0] / pow((Ni[0] + Ni[1]), 2); // d Xoo / d Ng
                dXsdXp[4 * 4 + 1] = -dXsdXp[3 * 4 + 1];               // d Xgo / d No
                dXsdXp[4 * 4 + 2] = -dXsdXp[3 * 4 + 2];               // d Xgo / d Ng

                 OCP_DBL tmp_new = (1 + rs) * (1 + rs);

                 mux[0] = -muoRs * tmp_new;         // dMuo / dXoo
                 mux[1] = muoRs * tmp_new;          // dMuo / dXgo

                 xix[0] = -xioRs * tmp_new;         // dXio / dXoo
                 xix[1] = xioRs * tmp_new;          // dXio / dXgo

                 rhox[0] = -rhooRs * tmp_new;       // dRhoo / dXoo
                 rhox[1] = rhooRs * tmp_new;        // dRhoo / dXgo

                break;
            }
        case PHASE_ODGW:
            {
            // saturated oil, gas, water
            
                phaseExist[0] = OCP_TRUE;
                phaseExist[1] = OCP_TRUE;
                phaseExist[2] = OCP_TRUE;

                // oil property
                OCP_DBL rs, rsP;
                PVCO.CalRhoXiMuRsDer(P, rho[0], xi[0], mu[0], rs, rhoP[0], xiP[0], muP[0], rsP);

                // gas property
                PVDG.CalRhoXiMuDer(P, rho[1], xi[1], mu[1], rhoP[1], xiP[1], muP[1]);

                // total
                xij[0 * 3 + 0] = 1 / (1 + rs);
                xij[0 * 3 + 1] = 1 - xij[0 * 3 + 0];
                xij[1 * 3 + 1] = 1;
                xij[2 * 3 + 2] = 1;

                vj[0] = Ni[0] * (1 + rs) / xi[0];
                vj[1] = (Ni[1] - rs * Ni[0]) / xi[1];
                vj[2] = Ni[2] / xi[2];
                vf    = vj[0] + vj[1] + vj[2];
                S[0]  = vj[0] / vf;
                S[1]  = vj[1] / vf;
                S[2]  = vj[2] / vf;

                vfP = Ni[0] * (rsP * xi[0] - (1 + rs) * xiP[0]) / (xi[0] * xi[0]) +
                    (-rsP * Ni[0] * xi[1] - (Ni[1] - rs * Ni[0]) * xiP[1]) / (xi[1] * xi[1]) +
                    -Ni[2] * xiP[2] / (xi[2] * xi[2]);

                vfi[0] = (1 + rs) / xi[0] - rs / xi[1];
                vfi[1] = 1 / xi[1];
                vfi[2] = 1 / xi[2];

                dXsdXp[0] = (Ni[0] * (rsP * xi[0] - (1 + rs) * xiP[0]) / (xi[0] * xi[0]) - S[0] * vfP) / vf; // dSo / dP
                dXsdXp[1] = ((1 + rs) / xi[0] - S[0] * vfi[0]) / vf;    // dSo / dNo
                dXsdXp[2] = -S[0] * vfi[1] / vf;                        // dSo / dNg
                dXsdXp[3] = -S[0] * vfi[2] / vf;                        // dSo / dNw

                dXsdXp[1 * 4 + 0] = ((-rsP * Ni[0] * xi[1] - (Ni[1] - rs * Ni[0]) * xiP[1]) / (xi[1] * xi[1]) - S[1] * vfP) / vf;  // dSg / dP                          // dSg / dP
                dXsdXp[1 * 4 + 1] = (-rs / xi[1] - S[1] * vfi[0]) / vf; // dSg / dNo
                dXsdXp[1 * 4 + 2] = (1 / xi[1] - S[1] * vfi[1]) / vf;   // dSg / dNg
                dXsdXp[1 * 4 + 3] = -S[1] * vfi[2] / vf;                // dSg / dNw

                dXsdXp[2 * 4 + 0] = (-Ni[2] * xiP[2] / (xi[2] * xi[2]) - S[2] * vfP) / vf; // dSw / dP
                dXsdXp[2 * 4 + 1] = -S[2] * vfi[0] / vf;               // dSw / dNo
                dXsdXp[2 * 4 + 2] = -S[2] * vfi[1] / vf;               // dSw / dNg
                dXsdXp[2 * 4 + 3] = (1 / xi[2] - S[2] * vfi[2]) / vf;  // dSw / dNw

                dXsdXp[3 * 4 + 0] = -rsP / ((1 + rs) * (1 + rs)); // d Xoo / dP
                dXsdXp[4 * 4 + 0] = -dXsdXp[3 * 4 + 0];           // d Xgo / dP

                break;
            }
    }

#ifdef Debug
    for (auto vj : dXsdXp) {
        if (!isfinite(vj)) {
            OCP_ABORT("Nan in dXsdXp!");
        }
    }
#endif
}

OCP_DBL
BOMixture_ODGW::XiPhase(const OCP_DBL& Pin,
                        const OCP_DBL& Tin,
                        const OCP_DBL* Ziin,
                        const USI&     tarPhase)
{
    if (tarPhase == GAS) {
        return PVDG.CalXiG(Pin);
    } else if (tarPhase == WATER) {
        return PVTW.CalXiW(Pin);
    } else {
        OCP_ABORT("Wrong tarPhase!");
    }
}

OCP_DBL
BOMixture_ODGW::RhoPhase(const OCP_DBL& Pin,
                         const OCP_DBL& Pbbin,
                         const OCP_DBL& Tin,
                         const OCP_DBL* Ziin,
                         const USI&     tarPhase)
{

    if (tarPhase == OIL) {
        return PVCO.CalRhoO(Pin, Pbbin);
    } else if (tarPhase == GAS) {
        return PVDG.CalRhoG(Pin);
    } else if (tarPhase == WATER) {
        return PVTW.CalRhoW(Pin);
    } else {
        OCP_ABORT("WRONG tarPhase!");
    }
}

void BOMixture_ODGW::SetupWellOpt(WellOpt&                  opt,
                                  const vector<SolventINJ>& sols,
                                  const OCP_DBL&            Psurf,
                                  const OCP_DBL&            Tsurf)
{
    const USI wellType = opt.WellType();
    if (wellType == INJ) {
        const string fluidName = opt.InjFluidType();
        opt.SetInjFactor(1.0);

        if (fluidName == "WAT") {
            vector<OCP_DBL> tmpZi({0, 0, 1});
            opt.SetInjZi(tmpZi);
            opt.SetInjProdPhase(WATER);
        } else if (fluidName == "GAS") {
            vector<OCP_DBL> tmpZi({0, 1, 0});
            opt.SetInjZi(tmpZi);
            opt.SetInjProdPhase(GAS);
        } else {
            OCP_ABORT("WRONG Injecting Fluid!");
        }

    } else if (wellType == PROD) {
        vector<OCP_DBL> tmpWght(3, 0);
        switch (opt.OptMode()) {
            case BHP_MODE:
                break;
            case ORATE_MODE:
                tmpWght[0] = 1;
                break;
            case GRATE_MODE:
                tmpWght[1] = 1;
                break;
            case WRATE_MODE:
                tmpWght[2] = 1;
                break;
            case LRATE_MODE:
                tmpWght[0] = tmpWght[2] = 1;
                break;
            default:
                OCP_ABORT("WRONG Opt Mode!");
                break;
        }
        opt.SetProdPhaseWeight(tmpWght);
    } else {
        OCP_ABORT("Wrong Well Type!");
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/