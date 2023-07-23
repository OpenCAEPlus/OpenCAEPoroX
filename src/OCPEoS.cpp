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
// OCPEoS_PR
/////////////////////////////////////////////////////


OCPEoS_PR::OCPEoS_PR(const ComponentParam& param, const USI& tarId)
{
    nc = param.numCom;
    if (param.Tc.activity)      Tc = param.Tc.data[tarId];
    else                        OCP_ABORT("TCRIT is Missing!");
    if (param.Pc.activity)      Pc = param.Pc.data[tarId];
    else                        OCP_ABORT("PCRIT is Missing!");
    if (param.Acf.activity)     acf = param.Acf.data[tarId];
    else                        OCP_ABORT("ACF is Missing!");
    if (param.OmegaA.activity)  OmegaA = param.OmegaA.data[tarId];
    else                        OmegaA.resize(nc, 0.457235529);
    if (param.OmegaB.activity)  OmegaB = param.OmegaB.data[tarId];
    else                        OmegaB.resize(nc, 0.077796074);
}


OCP_DBL OCPEoS_PR::SolEoS(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi)
{
    //OCP_DBL Aj = 0;
    //OCP_DBL Bj = 0;

    //for (USI i1 = 0; i1 < nc; i1++) {
    //    Bj += Bi[i1] * zi[i1];
    //    Aj += zi[i1] * zi[i1] * Ai[i1] * (1 - BIC[i1 * nc + i1]);

    //    for (USI i2 = 0; i2 < i1; i2++) {
    //        Aj += 2 * zi[i1] * zi[i2] * sqrt(Ai[i1] * Ai[i2]) * (1 - BIC[i1 * nc + i2]);
    //    }
    //}


    //const OCP_DBL a = (delta1 + delta2 - 1) * Bj - 1;
    //const OCP_DBL b =
    //    (Aj + delta1 * delta2 * Bj * Bj - (delta1 + delta2) * Bj * (Bj + 1));
    //const OCP_DBL c = -(Aj * Bj + delta1 * delta2 * Bj * Bj * (Bj + 1));

    //USI flag = CubicRoot(a, b, c, OCP_TRUE); // True with NT
    //if (flag == 1) {
    //    ZjT = Ztmp[0];
    //}
    //else {
    //    OCP_DBL zj1 = Ztmp[0];
    //    OCP_DBL zj2 = Ztmp[2];
    //    OCP_DBL dG = (zj2 - zj1) + log((zj1 - Bj) / (zj2 - Bj)) -
    //        Aj / (Bj * (delta2 - delta1)) *
    //        log((zj1 + delta1 * Bj) * (zj2 + delta2 * Bj) /
    //            ((zj1 + delta2 * Bj) * (zj2 + delta1 * Bj)));
    //    if (dG > 0)
    //        ZjT = zj1;
    //    else
    //        ZjT = zj2;
    //}
}


void OCPEoS_PR::CalAiBi(const OCP_DBL& P, const OCP_DBL& T)
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
        Pri   = P / Pc[i];
        Tri   = T / Tc[i];
        Ai[i] = OmegaA[i] * Pri / pow(Tri, 2) * pow((1 + mwi * (1 - sqrt(Tri))), 2);
        Bi[i] = OmegaB[i] * Pri / Tri;
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/23/2023      Create file                          */
/*----------------------------------------------------------------------------*/