/*! \file    MixtureThermal_K.cpp
 *  \brief   MixtureThermal_K class declaration
 *  \author  Shizhe Li
 *  \date    Nov/10/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// OpenCAEPoroX header files
#include "MixtureThermal.hpp"

MixtureThermal_K01::MixtureThermal_K01(const ParamReservoir& param, const USI& tarId)
{
    mixtureType = THERMAL;
    numPhase    = 2;
    numCom      = 2;

    eC.Setup(param.comsParam, tarId);
    OWTM.Setup(param, tarId);
}

void MixtureThermal_K01::Flash(const OCP_DBL& Pin,
                               const OCP_DBL& Tin,
                               const OCP_DBL* Niin)
{
    OWTM.Flash(Pin, Tin, Niin);
}

void MixtureThermal_K01::InitFlashIMPEC(const OCP_DBL& Pin,
                                        const OCP_DBL& Pbbin,
                                        const OCP_DBL& Tin,
                                        const OCP_DBL* Sjin,
                                        const OCP_DBL& Vpore,
                                        const OCP_DBL* Ziin,
                                        const OCP_USI& bId)
{
    OWTM.InitFlash(Pin, Tin, Sjin[1], Vpore); 
}

void MixtureThermal_K01::InitFlashFIM(const OCP_DBL& Pin,
                                      const OCP_DBL& Pbbin,
                                      const OCP_DBL& Tin,
                                      const OCP_DBL* Sjin,
                                      const OCP_DBL& Vpore,
                                      const OCP_DBL* Ziin,
                                      const OCP_USI& bId)
{
    OWTM.InitFlashDer(Pin, Tin, Sjin[1], Vpore);
}

void MixtureThermal_K01::FlashIMPEC(const OCP_DBL& Pin,
                                    const OCP_DBL& Tin,
                                    const OCP_DBL* Niin,
                                    const USI&     lastNP,
                                    const OCP_DBL* xijin,
                                    const OCP_USI& bId)
{
    OWTM.Flash(Pin, Tin, Niin);
}

void MixtureThermal_K01::FlashFIM(const OCP_DBL& Pin,
                                  const OCP_DBL& Tin,
                                  const OCP_DBL* Niin,
                                  const OCP_DBL* Sjin,
                                  const USI&     lastNP,
                                  const OCP_DBL* xijin,
                                  const OCP_USI& bId)
{
    OWTM.FlashDer(Pin, Tin, Niin);
}


OCP_DBL MixtureThermal_K01::CalInjWellEnthalpy(const OCP_DBL& Tin, const OCP_DBL* Ziin)
{
    return OWTM.CalEnthalpy(Tin, Ziin);
}

OCP_DBL MixtureThermal_K01::XiPhase(const OCP_DBL& Pin,
                                    const OCP_DBL& Tin,
                                    const OCP_DBL* Ziin,
                                    const USI&     tarPhase)
{
    return OWTM.CalXi(Pin, Tin, tarPhase);
}

OCP_DBL
MixtureThermal_K01::RhoPhase(const OCP_DBL& Pin,
                             const OCP_DBL& Pbb,
                             const OCP_DBL& Tin,
                             const OCP_DBL* Ziin,
                             const USI&     tarPhase)
{
    return OWTM.CalRho(Pin, Tin, tarPhase);
}

void MixtureThermal_K01::CalProdWeight(const OCP_DBL&         Pin,
                                       const OCP_DBL&         Tin,
                                       const OCP_DBL*         Niin,
                                       const vector<OCP_DBL>& prodPhase,
                                       vector<OCP_DBL>&       prodWeight)
{

    vector<OCP_DBL> factor(numPhase, 0);
    factor[0] = Niin[0] / (Niin[0] + Niin[1]) / OWTM.CalXi(Pin, Tin, OIL) / CONV1; // stb / lbmol
    factor[1] = Niin[1] / (Niin[0] + Niin[1]) / OWTM.CalXi(Pin, Tin, WATER) / CONV1; // stb  / lbmol

    OCP_DBL tmp = 0;
    for (USI i = 0; i < 2; i++) {
        tmp += factor[i] * prodPhase[i];
    }
    if (tmp < 1E-12 || !isfinite(tmp)) {
        OCP_ABORT("Wrong Condition!");
    }
    fill(prodWeight.begin(), prodWeight.end(), tmp);
}

void MixtureThermal_K01::CalProdRate(const OCP_DBL&   Pin,
                                     const OCP_DBL&   Tin,
                                     const OCP_DBL*   Niin,
                                     vector<OCP_DBL>& prodRate)
{
    prodRate[0] = Niin[0] / OWTM.CalXi(Pin, Tin, OIL) / CONV1; // stb
    prodRate[1] = Niin[1] / OWTM.CalXi(Pin, Tin, WATER) / CONV1; // stb
}

void MixtureThermal_K01::SetupWellOpt(WellOpt&                  opt,
                                      const vector<SolventINJ>& sols,
                                      const OCP_DBL&            Psurf,
                                      const OCP_DBL&            Tsurf)
{
    const USI wellType = opt.WellType();
    if (wellType == INJ) {
        const string    fluidName = opt.InjFluidType();
        vector<OCP_DBL> tmpZi(numCom, 0);
        if (fluidName == "WAT") {
            tmpZi.back() = 1;
            opt.SetInjProdPhase(WATER);
            // lbmol / ft3 -> lbmol  / bbl  for
            // injfluid Use flash in Bulk in surface condition
            OCP_DBL tmp = CONV1 * XiPhase(Psurf, Tsurf, &tmpZi[0], WATER);
            opt.SetInjFactor(tmp);
        } else {
            OCP_ABORT("WRONG Fluid Type!");
        }
        opt.SetInjZi(tmpZi);
    } else if (wellType == PROD) {
        vector<OCP_DBL> tmpWght(numPhase, 0);
        switch (opt.OptMode()) {
            case BHP_MODE:
                break;
            case ORATE_MODE:
                tmpWght[0] = 1;
                break;
            case WRATE_MODE:
                tmpWght[1] = 1;
                break;
            case LRATE_MODE:
                tmpWght[0] = 1;
                tmpWght[1] = 1;
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
/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           NOV/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/
