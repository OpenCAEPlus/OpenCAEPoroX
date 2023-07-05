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


void Miscible::InputParam(const Miscstr& misterm)
{
    ifMiscible  = misterm.ifMiscible;
    surTenRef   = misterm.surTenRef;
    surTenMaxPc = misterm.surTenPc;
    if (surTenMaxPc < 0) surTenMaxPc = surTenRef;
    surTenPc    = 1;
    Fkexp       = misterm.surTenExp;
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


void Miscible::CalMiscibleFactor(const OCP_USI& bId, const USI& mIndex)
{
    if (ifMiscible) {
        surTen[bId] = sT.CalSurfaceTension(mIndex);
        mF.CalculateMiscibleFactor(mIndex, surTen[bId], Fk[bId], Fp[bId]);
    }

}

/// Coats expression
OCP_BOOL Miscible::CalFkFp(const OCP_USI& n, OCP_DBL& fk, OCP_DBL& fp)
{
    if (Fk[n] < 0) {
        return OCP_FALSE;
    }
    else {
        fk = Fk[n];
        fp = Fp[n];
        return OCP_TRUE;
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/26/2022      Create file                          */
/*----------------------------------------------------------------------------*/