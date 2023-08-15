/*! \file    MixtureUnitComp.cpp
 *  \brief   MixtureUnitComp class definition for compositional models
 *  \author  Shizhe Li
 *  \date    Jan/05/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "MixtureUnitComp.hpp"


MixtureUnitComp::MixtureUnitComp(const ParamReservoir& rs_param, const USI& i, OptionalFeatures& opts)
{
    compM.Setup(rs_param, i);
    mixtureType = compM.MixtureType();
    vs          = &compM.GetVarSet();

    /// Optional Features
    // Skip stability analysis
    skipPSA         = &opts.skipPSA;
    skipMethodIndex = skipPSA->Setup(opts.numBulk, &compM);
    // Calculate surface tension
    surTen          = &opts.surTen;
    stMethodIndex   = surTen->Setup(rs_param, i, opts.numBulk, &compM);
    // Miscible Factor
    misFac          = &opts.misFac;
    mfMethodIndex   = misFac->Setup(rs_param, i, opts.numBulk, surTen);
}

void MixtureUnitComp::Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{
    compM.Flash(Pin, Tin, Niin);
}


void MixtureUnitComp::InitFlashIMPEC(const OCP_DBL& Pin,
                                 const OCP_DBL& Pbbin,
                                 const OCP_DBL& Tin,
                                 const OCP_DBL* Sjin,
                                 const OCP_DBL& Vpore,
                                 const OCP_DBL* Ziin,
                                 const OCP_USI& bId)
{
    compM.InitFlash(Pin, Tin, Sjin, Ziin, Vpore);
    skipPSA->CalSkipForNextStep(bId, skipMethodIndex);
    surTen->CalSurfaceTension(bId, stMethodIndex);
    misFac->CalMiscibleFactor(bId, mfMethodIndex);
}

void MixtureUnitComp::InitFlashFIM(const OCP_DBL& Pin,
                               const OCP_DBL& Pbbin,
                               const OCP_DBL& Tin,
                               const OCP_DBL* Sjin,
                               const OCP_DBL& Vpore,
                               const OCP_DBL* Ziin,
                               const OCP_USI& bId)
{
    compM.InitFlashDer(Pin, Tin, Sjin, Ziin, Vpore);
    skipPSA->CalSkipForNextStep(bId, skipMethodIndex);
    surTen->CalSurfaceTension(bId, stMethodIndex);
    misFac->CalMiscibleFactor(bId, mfMethodIndex);
}

void MixtureUnitComp::FlashIMPEC(const OCP_DBL& Pin,
                             const OCP_DBL& Tin,
                             const OCP_DBL* Niin,
                             const USI&     lastNP,
                             const OCP_DBL* xijin,
                             const OCP_USI& bId)
{
    const USI ftype = skipPSA->CalFtype(Pin, Tin, Niin, bId, skipMethodIndex);
    compM.Flash(Pin, Tin, Niin, ftype, lastNP, xijin);
    skipPSA->CalSkipForNextStep(bId, skipMethodIndex);
    surTen->CalSurfaceTension(bId, stMethodIndex);
    misFac->CalMiscibleFactor(bId, mfMethodIndex);
}

void MixtureUnitComp::FlashFIM(const OCP_DBL& Pin,
                           const OCP_DBL& Tin,
                           const OCP_DBL* Niin,
                           const OCP_DBL* Sjin,
                           const USI&     lastNP,
                           const OCP_DBL* xijin,
                           const OCP_USI& bId)
{
    const USI ftype = skipPSA->CalFtype(Pin, Tin, Niin, Sjin, compM.GetNumPhasePE(lastNP), bId, skipMethodIndex);
    compM.FlashDer(Pin, Tin, Niin, ftype, lastNP, xijin);
    skipPSA->CalSkipForNextStep(bId, skipMethodIndex);
    surTen->CalSurfaceTension(bId, stMethodIndex);
    misFac->CalMiscibleFactor(bId, mfMethodIndex);
}


OCP_DBL MixtureUnitComp::XiPhase(const OCP_DBL& Pin,
                         const OCP_DBL& Tin,
                         const vector<OCP_DBL>& Ziin,
                         const USI&     tarPhase)
{
    return compM.CalXi(Pin, Tin, &Ziin[0], tarPhase);
}

OCP_DBL MixtureUnitComp::RhoPhase(const OCP_DBL& Pin,
                      const OCP_DBL& Pbb,
                      const OCP_DBL& Tin,
                      const vector<OCP_DBL>& Ziin,
                      const USI&     tarPhase)
{
    return compM.CalRho(Pin, Tin, &Ziin[0], tarPhase);
}


void MixtureUnitComp::CalProdWeight(const OCP_DBL&     Pin,
                                const OCP_DBL&         Tin,
                                const OCP_DBL*         Niin,
                                const vector<OCP_DBL>& prodPhase,
                                vector<OCP_DBL>&       prodWeight)
{
    compM.Flash(Pin, Tin, Niin);

    const OCP_DBL   qt = vs->Nt;
    vector<OCP_DBL> factor(vs->np, 0);

    factor[0] = vs->vj[0] / qt / CONV1; // stb / lbmol
    factor[1] = vs->vj[1] / qt / 1000;  // Mscf / lbmol
    factor[2] = vs->xi[2] * vs->vj[2] / qt; // stb  / stb

    OCP_DBL tmp = 0;
    for (USI i = 0; i < 3; i++) {
        tmp += factor[i] * prodPhase[i];
    }
    if (tmp < 1E-12 || !isfinite(tmp)) {
        OCP_ABORT("Wrong Condition!");
    }
    fill(prodWeight.begin(), prodWeight.end(), tmp);
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jan/05/2022      Create file                          */
/*----------------------------------------------------------------------------*/