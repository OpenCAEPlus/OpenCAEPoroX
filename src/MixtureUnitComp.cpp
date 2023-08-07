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


MixtureUnitComp::MixtureUnitComp(const ParamReservoir& rs_param, const USI& i)
{
    compM.Setup(rs_param, i);
    mixtureType = compM.MixtureType();
    vs          = &compM.GetVarSet();
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
    skipPSA.CalSkipForNextStep(bId, skipMethodIndex, compM);
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
    skipPSA.CalSkipForNextStep(bId, skipMethodIndex, compM);
}

void MixtureUnitComp::FlashIMPEC(const OCP_DBL& Pin,
                             const OCP_DBL& Tin,
                             const OCP_DBL* Niin,
                             const USI&     lastNP,
                             const OCP_DBL* xijin,
                             const OCP_USI& bId)
{
    const USI ftype = skipPSA.CalFtype(Pin, Tin, Niin, bId, skipMethodIndex);
    compM.Flash(Pin, Tin, Niin, ftype, lastNP, xijin);
    skipPSA.CalSkipForNextStep(bId, skipMethodIndex, compM);
}

void MixtureUnitComp::FlashFIM(const OCP_DBL& Pin,
                           const OCP_DBL& Tin,
                           const OCP_DBL* Niin,
                           const OCP_DBL* Sjin,
                           const USI&     lastNP,
                           const OCP_DBL* xijin,
                           const OCP_USI& bId)
{
    const USI ftype = skipPSA.CalFtype(Pin, Tin, Niin, Sjin, compM.GetNumPhasePE(lastNP), bId, skipMethodIndex);
    compM.FlashDer(Pin, Tin, Niin, ftype, lastNP, xijin);
    skipPSA.CalSkipForNextStep(bId, skipMethodIndex, compM);
}


OCP_DBL
MixtureUnitComp::XiPhase(const OCP_DBL& Pin,
                         const OCP_DBL& Tin,
                         const vector<OCP_DBL>& Ziin,
                         const USI&     tarPhase)
{
    return compM.CalXi(Pin, Tin, &Ziin[0], tarPhase);
}

OCP_DBL
MixtureUnitComp::RhoPhase(const OCP_DBL& Pin,
                      const OCP_DBL& Pbb,
                      const OCP_DBL& Tin,
                      const vector<OCP_DBL>& Ziin,
                      const USI&     tarPhase)
{
    return compM.CalRho(Pin, Tin, &Ziin[0], tarPhase);
}

void MixtureUnitComp::SetupWellOpt(WellOpt&                  opt,
                               const vector<SolventINJ>& sols,
                               const OCP_DBL&            Psurf,
                               const OCP_DBL&            Tsurf)
{
    const USI wellType = opt.WellType();
    if (wellType == INJ) {
        const string    fluidName = opt.InjFluidType();
        vector<OCP_DBL> tmpZi(vs->nc, 0);
        if (fluidName == "WAT") {
            tmpZi.back() = 1;
            opt.SetInjProdPhase(WATER);
            opt.SetInjFactor(1.0);
        } else {
            // inj phase is gas
            opt.SetInjProdPhase(GAS);
            const USI len = sols.size();
            for (USI i = 0; i < len; i++) {
                if (fluidName == sols[i].name) {
                    tmpZi = sols[i].data;
                    tmpZi.resize(vs->nc);
                    // Convert volume units Mscf/stb to molar units lbmoles for
                    // injfluid Use flash in Bulk in surface condition
                    OCP_DBL tmp = 1000 * XiPhase(Psurf, Tsurf, tmpZi, GAS);
                    opt.SetInjFactor(tmp);
                    break;
                }
                if (i == len - 1) {
                    OCP_ABORT("Wrong FluidType!");
                }
            }
        }
        opt.SetInjZi(tmpZi);
    } else if (wellType == PROD) {
        vector<OCP_DBL> tmpWght(vs->np, 0);
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
                tmpWght[0] = 1;
                tmpWght[2] = 1;
                break;
            default:
                OCP_ABORT("WRONG optMode");
                break;
        }
        opt.SetProdPhaseWeight(tmpWght);
    } else {
        OCP_ABORT("Wrong Well Type!");
    }
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

void MixtureUnitComp::CalProdRate(const OCP_DBL&   Pin,
                              const OCP_DBL&   Tin,
                              const OCP_DBL*   Niin,
                              vector<OCP_DBL>& prodRate)
{
    compM.Flash(Pin, Tin, Niin);

    prodRate[0] = vs->vj[0] / CONV1; // stb
    prodRate[1] = vs->vj[1] / 1000;  // Mscf
    prodRate[2] = vs->vj[2] * vs->xi[2]; // stb
}


/////////////////////////////////////////////////////////////////////
// Optional Features
/////////////////////////////////////////////////////////////////////

void MixtureUnitComp::SetupOptionalFeatures(OptionalFeatures& optFeatures)
{
    skipMethodIndex = skipPSA.Setup(optFeatures.numBulk, compM.GetNPmaxPE(), compM.GetNCPE());
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jan/05/2022      Create file                          */
/*----------------------------------------------------------------------------*/