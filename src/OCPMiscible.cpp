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


void MisFacMethod01::CalculateMiscibleFactor(const OCP_DBL& st, OCP_DBL& fk, OCP_DBL& fp) const
{
    if (st >= stref) {  // InMiscible
        fk = -1.0;
        fp = -1.0;
    }
    else {             // Miscible
        fk = min(1.0, pow(st / stref, fkExp));
        fp = min(stPcf, st / stref);
    }
}


USI MiscibleFcator::Setup(const MisFacMethodParams& param)
{
    if (OCP_FALSE) {

    }
    else {
        // Default option
        if (param.param01.IfUse()) {
            mfMethod.push_back(new MisFacMethod01(param.param01));
        }
        else {
            OCP_ABORT("NO MATHCHED METHOD IN MISCIBLE FACTOR METHOD!");
        }      
    }
    
    return mfMethod.size() - 1;
}


void MisCurveMethod01::CurveCorrect(const OCP_DBL& Fk, const OCP_DBL& Fp)
{
    if (Fk > -TINY) {
        // miscible

        OCP3PFVarSet& vs = PF3->GetVarSet();

        // for permeability
        OCP_DBL tmp;
        const OCP_DBL krgt = PF3->CalKrg((1 - vs.Sw), tmp);
        const OCP_DBL krh = 0.5 * (vs.krow + krgt);

        vs.kro = Fk * vs.kro + (1 - Fk) * krh * vs.So / (1 - vs.Sw);
        vs.krg = Fk * vs.krg + (1 - Fk) * krh * vs.Sg / (1 - vs.Sw);

        // for capillary pressure
        vs.Pcgo *= Fp;
    }
}


void MisCurveMethod01::CurveCorrectDer(const OCP_DBL& Fk, const OCP_DBL& Fp)
{
    if (Fk > -TINY) {
        // miscible

        OCP3PFVarSet& vs = PF3->GetVarSet();

        // for permeability
        OCP_DBL       dkrgd1_Sw = 0;
        const OCP_DBL krgt = PF3->CalKrg((1 - vs.Sw), dkrgd1_Sw);
        const OCP_DBL krh = 0.5 * (vs.krow + krgt);

        vs.kro = Fk * vs.kro + (1 - Fk) * krh * vs.So / (1 - vs.Sw);
        vs.krg = Fk * vs.krg + (1 - Fk) * krh * vs.Sg / (1 - vs.Sw);

        const OCP_DBL dKrhdSo = 0.5 * vs.dKrowdSo;
        const OCP_DBL dKrhdSw = 0.5 * (vs.dKrowdSw - dkrgd1_Sw);      

        vs.dKrodSo = Fk * vs.dKrodSo + (1 - Fk) * (dKrhdSo * vs.So + krh) / (1 - vs.Sw);
        vs.dKrodSg = Fk * vs.dKrodSg;
        vs.dKrodSw = Fk * vs.dKrodSw + (1 - Fk) * vs.So * (dKrhdSw * (1 - vs.Sw) + krh) / ((1 - vs.Sw) * (1 - vs.Sw));
        vs.dKrgdSo = (1 - Fk) * dKrhdSo * vs.Sg / (1 - vs.Sw);
        vs.dKrgdSg = Fk * vs.dKrgdSg + (1 - Fk) * krh / (1 - vs.Sw);
        vs.dKrgdSw = (1 - Fk) * vs.Sg * (dKrhdSw * (1 - vs.Sw) + krh) / ((1 - vs.Sw) * (1 - vs.Sw));

        // for capillary pressure
        vs.Pcgo     *= Fp;
        vs.dPcgodSo *= Fp;
        vs.dPcgodSg *= Fp;
        vs.dPcgodSw *= Fp;
    }
}


USI MiscibleCurve::Setup(OCP3PhaseFlow* pf3)
{
    mcMethod.push_back(new MisCurveMethod01(pf3));
    
    return mcMethod.size() - 1;
}


void Miscible::InputParam(const OCP_BOOL& ifmiscible)
{
    ifMiscible  = ifmiscible;
}

USI Miscible::Setup(const OCP_USI& numBulk, const SurTenMethodParams& stparams, const MisFacMethodParams& mfparams)
{
    if (ifMiscible) {
        surTen.resize(numBulk);
        Fk.resize(numBulk);
        Fp.resize(numBulk);

        // Last time step
        lsurTen.resize(numBulk);


        // Setup Method
        // Setup surface tension method
        USI stmIndex = sT.Setup(stparams);
        // Setup miscible factor method
        USI mfmIndex = mF.Setup(mfparams);

        OCP_ASSERT(stmIndex == mfmIndex, "Mismatched Methods!");

        return stmIndex;
    }
    return 0;
}


USI Miscible::Setup(OCP3PhaseFlow* pf3)
{
    if (ifMiscible) {
        return mC.Setup(pf3);
    }
    return 0;
}


void Miscible::CalMiscibleFactor(const OCP_USI& bId, const USI& mIndex)
{
    if (ifMiscible) {
        surTen[bId] = sT.CalSurfaceTension(mIndex);
        mF.CalculateMiscibleFactor(mIndex, surTen[bId], Fk[bId], Fp[bId]);
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/26/2022      Create file                          */
/*----------------------------------------------------------------------------*/