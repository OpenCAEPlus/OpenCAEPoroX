/*! \file    OCPMiscible.cpp
 *  \brief   OCPMiscible class declaration
 *  \author  Shizhe Li
 *  \date    Dec/26/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPMiscible.hpp"


MisFacMethod01::MisFacMethod01(const Miscstr& param, MisFacVarSet& mfvs)
{
    stref = param.surTenRef;
    fkExp = param.surTenExp;
    stPcf = param.surTenPc;
    if (stPcf < 0) stPcf = stref;

    mfvs.Fk.resize(mfvs.nb, 0);
    mfvs.Fp.resize(mfvs.nb, 0);
}

void MisFacMethod01::CalculateMiscibleFactor(const OCP_USI& bId, const SurTenVarSet& stvs, MisFacVarSet& mfvs) const
{
    const OCP_DBL& st = stvs.surTen[bId];
    if (st >= stref) {  // InMiscible
        mfvs.Fk[bId] = -1.0;
        mfvs.Fp[bId] = -1.0;
    }
    else {             // Miscible
        mfvs.Fk[bId] = min(static_cast<OCP_DBL>(1.0), pow(st / stref, fkExp));
        mfvs.Fp[bId] = min(stPcf, st / stref);
    }
}


USI MiscibleFactor::Setup(const ParamReservoir& param, const USI& i, const OCP_USI& nb, const SurfaceTension* st)
{
    if (param.miscstr.ifMiscible) {
        if (!st->IfUse()) {
            OCP_ABORT("Surface Tension is not set!");
        }
        else {
            surTen = st;
        }
        ifUse = OCP_TRUE;
        vs.SetNb(nb);
        mfMethod.push_back(new MisFacMethod01(param.miscstr, vs));
    }

    return mfMethod.size() - 1;
}


void MisCurveMethod01::CurveCorrect(const OCP_USI& bId, const MisFacVarSet& mfvs)
{
    const OCP_DBL& Fk = mfvs.Fk[bId];
    const OCP_DBL& Fp = mfvs.Fp[bId];

    if (Fk > -TINY) {
        // miscible

        OCPFlowVarSet& vs = flow->GetVarSet();
        const INT&     o  = vs.o;
        const INT&     g  = vs.g;
        const INT&     w  = vs.w;

        // for permeability
        OCP_DBL tmp;
        const OCP_DBL krgt = flow->CalKrg((1 - vs.S[w]), tmp);
        const OCP_DBL krh = 0.5 * (vs.krow + krgt);

        vs.kr[o] = Fk * vs.kr[o] + (1 - Fk) * krh * vs.S[o] / (1 - vs.S[w]);
        vs.kr[g] = Fk * vs.kr[g] + (1 - Fk) * krh * vs.S[g] / (1 - vs.S[w]);

        // for capillary pressure
        vs.Pc[g] *= Fp;
    }
}


void MisCurveMethod01::CurveCorrectDer(const OCP_USI& bId, const MisFacVarSet& mfvs)
{
    const OCP_DBL& Fk = mfvs.Fk[bId];
    const OCP_DBL& Fp = mfvs.Fp[bId];

    if (Fk > -TINY) {
        // miscible

        OCPFlowVarSet& vs = flow->GetVarSet();
        const INT&     o  = vs.o;
        const INT&     g  = vs.g;
        const INT&     w  = vs.w;

        // for permeability
        OCP_DBL       dkrgd1_Sw = 0;
        const OCP_DBL krgt = flow->CalKrg((1 - vs.S[w]), dkrgd1_Sw);
        const OCP_DBL krh = 0.5 * (vs.krow + krgt);

        vs.kr[o] = Fk * vs.kr[o] + (1 - Fk) * krh * vs.S[o] / (1 - vs.S[w]);
        vs.kr[g] = Fk * vs.kr[g] + (1 - Fk) * krh * vs.S[g] / (1 - vs.S[w]);

        const OCP_DBL dKrhdSo = 0.5 * vs.dKrowdSo;
        const OCP_DBL dKrhdSw = 0.5 * (vs.dKrowdSw - dkrgd1_Sw);      

        vs.dKrdS[vs.oo] = Fk * vs.dKrdS[vs.oo] + (1 - Fk) * (dKrhdSo * vs.S[o] + krh) / (1 - vs.S[w]);
        vs.dKrdS[vs.og] = Fk * vs.dKrdS[vs.og];
        vs.dKrdS[vs.ow] = Fk * vs.dKrdS[vs.ow] + (1 - Fk) * vs.S[o] * (dKrhdSw * (1 - vs.S[w]) + krh) / ((1 - vs.S[w]) * (1 - vs.S[w]));
        vs.dKrdS[vs.go] = (1 - Fk) * dKrhdSo * vs.S[g] / (1 - vs.S[w]);
        vs.dKrdS[vs.gg] = Fk * vs.dKrdS[vs.gg] + (1 - Fk) * krh / (1 - vs.S[w]);
        vs.dKrdS[vs.gw] = (1 - Fk) * vs.S[g] * (dKrhdSw * (1 - vs.S[w]) + krh) / ((1 - vs.S[w]) * (1 - vs.S[w]));

        // for capillary pressure
        vs.Pc[g]        *= Fp;
        vs.dPcdS[vs.go] *= Fp;
        vs.dPcdS[vs.gg] *= Fp;
        vs.dPcdS[vs.gw] *= Fp;
    }
}


USI MiscibleCurve::Setup(OCPFlow* flowin, const MiscibleFactor* misfactor)
{
    if (!misfactor->IfUse()) {
        ifUse = OCP_FALSE;
        return 0;
    }
    else {
        misFac = misfactor;
        ifUse  = OCP_TRUE;
        mcMethod.push_back(new MisCurveMethod01(flowin));
        return mcMethod.size() - 1;
    }  
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/26/2022      Create file                          */
/*----------------------------------------------------------------------------*/