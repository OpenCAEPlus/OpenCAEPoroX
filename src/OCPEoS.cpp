/*! \file    OCPEoS.cpp
 *  \brief   OCPEoS class declaration
 *  \author  Shizhe Li
 *  \date    Jul/23/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPEoS.hpp"

/////////////////////////////////////////////////////
// EoS_PR
/////////////////////////////////////////////////////


EoS_PR::EoS_PR(const ComponentParam& param, const USI& tarId)
{
    nc = param.numCom;
    if (param.Tc.activity)      Tc = param.Tc.data[tarId];
    else                        OCP_ABORT("TCRIT is Missing!");
    if (param.Pc.activity)      Pc = param.Pc.data[tarId];
    else                        OCP_ABORT("PCRIT is Missing!");
    if (param.OmegaA.activity)  OmegaA = param.OmegaA.data[tarId];
    else                        OmegaA.resize(nc, 0.457235529);
    if (param.OmegaB.activity)  OmegaB = param.OmegaB.data[tarId];
    else                        OmegaB.resize(nc, 0.077796074);
    if (param.Acf.activity)     acf = param.Acf.data[tarId];
    else                        OCP_ABORT("ACF is Missing!");

    BIC.resize(nc * nc, 0);
    if (param.BIC[tarId].size() != nc * nc) {
        USI iter = 0;
        for (USI i = 1; i < nc; i++) {
            for (USI j = 0; j < i; j++) {
                BIC[i * nc + j] = param.BIC[tarId][iter];
                BIC[j * nc + i] = BIC[i * nc + j];
                iter++;
            }
        }
    }
    else {
        BIC = param.BIC[tarId];
    }

    if (param.Vshift.activity) {
        Vshift = param.Vshift.data[tarId];
        for (USI i = 0; i < nc; i++)
            Vshift[i] *= (GAS_CONSTANT * OmegaB[i] * Tc[i] / Pc[i]);
    }
    else
        Vshift.resize(nc, 0);

    Ai.resize(nc);
    Bi.resize(nc);
    Z.resize(3);
    Ax.resize(nc);
    Bx.resize(nc);
    Zx.resize(nc);
    An.resize(nc);
    Bn.resize(nc);
    Zn.resize(nc);
    x.resize(nc);
}


void EoS_PR::CalAiBi(const OCP_DBL& P, const OCP_DBL& T)
{
    OCP_DBL mwi, Pri, Tri;
    for (USI i = 0; i < nc; i++) {
        if (acf[i] <= 0.49) {
            mwi = 0.37464 + 1.54226 * acf[i] - 0.26992 * pow(acf[i], 2);
        }
        else {
            mwi = 0.379642 + 1.48503 * acf[i] - 0.164423 * pow(acf[i], 2) +
                0.016667 * pow(acf[i], 3);
        }
        Pri = P / Pc[i];
        Tri = T / Tc[i];
        Ai[i] = OmegaA[i] * Pri / pow(Tri, 2) * pow((1 + mwi * (1 - sqrt(Tri))), 2);
        Bi[i] = OmegaB[i] * Pri / Tri;
    }
}


void EoS_PR::CalAjBj(const vector<OCP_DBL>& x)
{
    Aj = 0;
    Bj = 0;
    for (USI i1 = 0; i1 < nc; i1++) {
        Bj += Bi[i1] * x[i1];
        Aj += x[i1] * x[i1] * Ai[i1] * (1 - BIC[i1 * nc + i1]);

        for (USI i2 = 0; i2 < i1; i2++) {
            Aj += 2 * x[i1] * x[i2] * sqrt(Ai[i1] * Ai[i2]) * (1 - BIC[i1 * nc + i2]);
        }
    }
}


void EoS_PR::CalZj(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x)
{
    const OCP_DBL a = (delta1 + delta2 - 1) * Bj - 1;
    const OCP_DBL b = (Aj + delta1 * delta2 * Bj * Bj - (delta1 + delta2) * Bj * (Bj + 1));
    const OCP_DBL c = -(Aj * Bj + delta1 * delta2 * Bj * Bj * (Bj + 1));

    const USI flag = CubicRoot(a, b, c, OCP_TRUE, Z); // True with NT
    if (flag == 1) {
        Zj = Z[0];
    }
    else {
        OCP_DBL zj1 = Z[0];
        OCP_DBL zj2 = Z[2];
        OCP_DBL dG = (zj2 - zj1) + log((zj1 - Bj) / (zj2 - Bj)) -
            Aj / (Bj * (delta2 - delta1)) *
            log((zj1 + delta1 * Bj) * (zj2 + delta2 * Bj) /
                ((zj1 + delta2 * Bj) * (zj2 + delta1 * Bj)));
        if (dG > 0)
            Zj = zj1;
        else
            Zj = zj2;
    }
}


void EoS_PR::CalAjBjZj(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x)
{
    CalAiBi(P, T);
    CalAjBj(x);
    CalZj(P, T, x);
}


void EoS_PR::CalAxBxZx(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x)
{
    Bx = Bi;
    for (USI i = 0; i < nc; i++) {
        OCP_DBL tmp = 0;
        for (USI k = 0; k < nc; k++) {
            tmp += x[k] * (1 - BIC[i * nc + k]) * sqrt(Ai[i] * Ai[k]);
        }
        Ax[i] = 2 * tmp;
        Zx[i] = ((Bj - Zj) * Ax[i] + ((Aj + (delta1 * delta2) * (3 * Bj * Bj + 2 * Bj))
              + ((delta1 + delta2) * (2 * Bj + 1) - 2 * (delta1 * delta2) * Bj) * Zj
              - ((delta1 + delta2) - 1) * Zj * Zj) * Bx[i])
              / (3 * Zj * Zj + 2 * (((delta1 + delta2) - 1) * Bj - 1) * Zj
              + (Aj + (delta1 * delta2) * Bj * Bj - (delta1 + delta2) * Bj * (Bj + 1)));
    }
}


void EoS_PR::CalAnBnZn(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& n, const OCP_DBL& nt)
{
    for (USI i = 0; i < nc; i++) {
        OCP_DBL tmp = 0;
        for (USI m = 0; m < nc; m++) {
            tmp += (1 - BIC[i * nc + m]) * sqrt(Ai[i] * Ai[m]) * n[m];
        }
        An[i] = 2 / nt * (tmp - Aj);
        Bn[i] = 1 / nt * (Bi[i] - Bj);
        Zn[i] = ((Bj - Zj) * An[i] + ((Aj + delta1 * delta2 * (3 * Bj * Bj + 2 * Bj))
              + ((delta1 + delta2) * (2 * Bj + 1) - 2 * delta1 * delta2 * Bj) * Zj
              - (delta1 + delta2 - 1) * Zj * Zj) * Bn[i])
              / (3 * Zj * Zj + 2 * ((delta1 + delta2 - 1) * Bj - 1) * Zj
              + (Aj + (delta1 * delta2) * Bj * Bj - (delta1 + delta2) * Bj * (Bj + 1)));
    }
}


void EoS_PR::CalApBpZp(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x)
{
    Ap = Aj / P;
    Bp = Bj / P;
    Zp = ((Bj - Zj) * Ap + ((Aj + delta1 * delta2 * (3 * Bj * Bj + 2 * Bj))
        + ((delta1 + delta2) * (2 * Bj + 1) - 2 * delta1 * delta2 * Bj) * Zj
        - (delta1 + delta2 - 1) * Zj * Zj) * Bp)
        / (3 * Zj * Zj + 2 * ((delta1 + delta2 - 1) * Bj - 1) * Zj
            + (Aj + delta1 * delta2 * Bj * Bj - (delta1 + delta2) * Bj * (Bj + 1)));
}


void EoS_PR::CalFug(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x,
    vector<OCP_DBL>& fug)
{
    CalAjBjZj(P, T, x);

    for (USI i = 0; i < nc; i++) {
        OCP_DBL tmp = 0;
        for (USI k = 0; k < nc; k++) {
            tmp += 2 * (1 - BIC[i * nc + k]) * sqrt(Ai[i] * Ai[k]) * x[k];
        }
        fug[i] = exp(Bi[i] / Bj * (Zj - 1) - log(Zj - Bj) -
            Aj / (delta1 - delta2) / Bj * (tmp / Aj - Bi[i] / Bj) *
            log((Zj + delta1 * Bj) / (Zj + delta2 * Bj))) * x[i] * P;
    }
}


void EoS_PR::CalFugPhi(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x, 
                       vector<OCP_DBL>& fug, vector<OCP_DBL>& phi)
{
    CalAjBjZj(P, T, x);

    for (USI i = 0; i < nc; i++) {
        OCP_DBL tmp = 0;
        for (USI k = 0; k < nc; k++) {
            tmp += 2 * (1 - BIC[i * nc + k]) * sqrt(Ai[i] * Ai[k]) * x[k];
        }
        phi[i] = exp(Bi[i] / Bj * (Zj - 1) - log(Zj - Bj) -
            Aj / (delta1 - delta2) / Bj * (tmp / Aj - Bi[i] / Bj) *
            log((Zj + delta1 * Bj) / (Zj + delta2 * Bj)));
        fug[i] = phi[i] * x[i] * P;
    }
}


void EoS_PR::CalFugX(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x, vector<OCP_DBL>& fugx)
{

    CalAjBjZj(P, T, x);
    CalAxBxZx(P, T, x);

    OCP_DBL C, E, G;
    OCP_DBL Cxk, Dxk, Exk, Gxk;
    OCP_DBL aik;
    G = (Zj + delta1 * Bj) / (Zj + delta2 * Bj);

    for (USI i = 0; i < nc; i++) {
        // i th fugacity
        C = x[i] * P / (Zj - Bj);
        E = -Aj / ((delta1 - delta2) * Bj) * (Ax[i] / Aj - Bx[i] / Bj);

        for (USI k = 0; k < nc; k++) {
            // k th components
            aik = (1 - BIC[i * nc + k]) * sqrt(Ai[i] * Ai[k]);

            Cxk = ((Zj - Bj) * delta(i, k) - x[i] * (Zx[k] - Bx[k])) * P / ((Zj - Bj) * (Zj - Bj));
            Dxk = Bx[i] / Bj * (Zx[k] - Bx[k] * (Zj - 1) / Bj);
            Exk = (2 * (Aj / Bj * Bx[k] * Bx[i] + Bj * aik) - Ax[i] * Bx[k] - Ax[k] * Bi[i]) / (Bj * Bj);
            Exk /= -(delta1 - delta2);
            Gxk = (delta1 - delta2) / (Zj + delta2 * Bj) / (Zj + delta2 * Bj) * (Zj * Bx[k] - Zx[k] * Bj);
            fugx[i * nc + k] = 1 / C * Cxk + Dxk + Exk * log(G) + E / G * Gxk;
        }
    }

#ifdef OCP_NANCHECK
	if (!CheckNan(fugx.size(), &fugx[0])) {
		OCP_ABORT("INF or NAN in fugX !");
	}
#endif // NANCHECK
}


void EoS_PR::CalFugN(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& n, vector<OCP_DBL>& fugn)
{
    OCP_DBL nt = 0;
    for (USI i = 0; i < nc; i++) nt += n[i];
    for (USI i = 0; i < nc; i++) x[i] = n[i] / nt;

    CalAjBjZj(P, T, x);
    CalAnBnZn(P, T, n, nt);

    OCP_DBL C, E, G;
    OCP_DBL Cnk, Dnk, Enk, Gnk;
    OCP_DBL tmp, aik;


	G = (Zj + delta1 * Bj) / (Zj + delta2 * Bj);
	for (USI i = 0; i < nc; i++) {
		// i th fugacity
		C = n[i] * P / (Zj - Bj);
		// D = Bi[i] / Bj * (Zj - 1);
		tmp = 0;
		for (USI k = 0; k < nc; k++) {
			tmp += (1 - BIC[i * nc + k]) * sqrt(Ai[i] * Ai[k]) * n[k];
		}
		E = -Aj / ((delta1 - delta2) * Bj) * (2 * tmp / Aj - Bi[i] / Bj);

		for (USI k = 0; k <= i; k++) {
			// k th components

			aik = (1 - BIC[i * nc + k]) * sqrt(Ai[i] * Ai[k]);

			Cnk = P / (Zj - Bj) / (Zj - Bj) *
				((Zj - Bj) / nt * (delta(i, k) - n[i]) -
					n[i] * (Zn[k] - Bn[k]));
			Dnk = Bi[i] / Bj * (Zn[k] - (Bi[k] - Bj) * (Zj - 1) / (nt * Bj));
			Gnk = (delta1 - delta2) / ((Zj + delta2 * Bj) * (Zj + delta2 * Bj)) *
				(Bn[k] * Zj - Zn[k] * Bj);
			Enk = -1 / (delta1 - delta2) / (Bj * Bj) *
				(2 * (Bj * aik / nt + Bn[k] * (Bi[i] * Aj / Bj - tmp)) -
					An[k] * Bi[i] - Aj * Bi[i] / nt);
			Enk -= E / nt;
			fugn[i * nc + k] = 1 / C * Cnk + Dnk + Enk * log(G) + E / G * Gnk;
			fugn[k * nc + i] = fugn[i * nc + k];
		}
	}

#ifdef OCP_NANCHECK
    if (!CheckNan(fugn.size(), &fugn[0])) {
        OCP_ABORT("INF or NAN in fugN !");
    }
#endif // NANCHECK
}


void EoS_PR::CalFugP(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x, vector<OCP_DBL>& fugP)
{

    CalAjBjZj(P, T, x);
    CalApBpZp(P, T, x);

    OCP_DBL C, E, G;
    OCP_DBL Cp, Dp, Gp;

	G  = (Zj + delta1 * Bj) / (Zj + delta2 * Bj);
	Gp = (delta1 - delta2) / ((Zj + delta2 * Bj) * (Zj + delta2 * Bj)) *
		(Bp * Zj - Zp * Bj);
	for (USI i = 0; i < nc; i++) {

		C = P / (Zj - Bj);
		// D = Bi[i] / Bj * (Zj - 1);

        OCP_DBL tmp = 0;
		for (USI m = 0; m < nc; m++) {
			tmp += (1 - BIC[i * nc + m]) * sqrt(Ai[i] * Ai[m]) * x[m];
		}

		E  = -Aj / ((delta1 - delta2) * Bj) * (2 * tmp / Aj - Bi[i] / Bj);
		Cp = ((Zj - Bj) - P * (Zp - Bp)) / ((Zj - Bj) * (Zj - Bj));
		Dp = Bi[i] / Bj * Zp;
		// Ep = 0;

		// Attention that if x[i] = 0, then fugp[i] = nan
		// but Cp also contains x[i], so eliminate it first
		fugP[i] = Cp / C + Dp + E / G * Gp;
	}

#ifdef OCP_NANCHECK
    if (!CheckNan(fugP.size(), &fugP[0])) {
    	OCP_ABORT("INF or NAN in fugP !");
    }
#endif // NANCHECK
}


OCP_DBL EoS_PR::CalVj(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x)
{

}


OCP_DBL EoS_PR::CalVjDer(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x, OCP_DBL& xiP, OCP_DBL* xix)
{
    const OCP_DBL CgTP = GAS_CONSTANT * T / P;
    OCP_DBL tmp;

    CalAjBjZj(P, T, x);
    CalAxBxZx(P, T, x);
    CalApBpZp(P, T, x);

    OCP_DBL v = Zj * GAS_CONSTANT * T / P;
    for (USI i = 0; i < nc; i++) {
        v -= x[i] * Vshift[i];
    }
    const OCP_DBL xi = 1 / v;

	tmp = -xi * xi;
	for (USI i = 0; i < nc; i++) {
		xix[i] = tmp * (CgTP * Zx[i] - Vshift[i]);
	}
	xiP = tmp * CgTP * (Zp - Zj / P);

    return xi;
}


void EoSCalculation::Setup(const ComponentParam& param, const USI& tarId)
{
    eos = new EoS_PR(param, tarId);
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/23/2023      Create file                          */
/*----------------------------------------------------------------------------*/