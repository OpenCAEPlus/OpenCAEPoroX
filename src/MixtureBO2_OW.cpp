/*! \file    MixtureBO2_OW.cpp
 *  \brief   Used for Black Oil model, where only water and oil exist.
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
// BOMixture_OW
///////////////////////////////////////////////

BOMixture_OW::BOMixture_OW(const ParamReservoir& rs_param, const USI& i)
{
    mixtureType = BLKOIL_OW;
    numPhase    = 2;
    numCom      = 2;

    OWM.Setup(rs_param, i);
}

void BOMixture_OW::Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{
    OWM.Flash(Pin, Niin);  
}

void BOMixture_OW::InitFlashIMPEC(const OCP_DBL& Pin,
                                  const OCP_DBL& Pbbin,
                                  const OCP_DBL& Tin,
                                  const OCP_DBL* Sjin,
                                  const OCP_DBL& Vpore,
                                  const OCP_DBL* Ziin,
                                  const OCP_USI& bId)
{
    OWM.InitFlash(Pin, Sjin[1], Vpore);   
}

void BOMixture_OW::InitFlashFIM(const OCP_DBL& Pin,
                                const OCP_DBL& Pbbin,
                                const OCP_DBL& Tin,
                                const OCP_DBL* Sjin,
                                const OCP_DBL& Vpore,
                                const OCP_DBL* Ziin,
                                const OCP_USI& bId)
{
    OWM.InitFlashDer(Pin, Sjin[1], Vpore);   
}

void BOMixture_OW::FlashIMPEC(const OCP_DBL& Pin,
                              const OCP_DBL& Tin,
                              const OCP_DBL* Niin,
                              const USI&     lastNP,
                              const OCP_DBL* xijin,
                              const OCP_USI& bId)
{
    OWM.Flash(Pin, Niin); 
}

void BOMixture_OW::FlashFIM(const OCP_DBL& Pin,
                            const OCP_DBL& Tin,
                            const OCP_DBL* Niin,
                            const OCP_DBL* Sjin,
                            const USI&     lastNP,
                            const OCP_DBL* xijin,
                            const OCP_USI& bId)
{
    OWM.FlashDer(Pin, Niin);
}

OCP_DBL BOMixture_OW::XiPhase(const OCP_DBL& Pin,
                              const OCP_DBL& Tin,
                              const OCP_DBL* Ziin,
                              const USI&     tarPhase)
{
    return OWM.CalXi(Pin, tarPhase);
}

OCP_DBL BOMixture_OW::RhoPhase(const OCP_DBL& Pin,
                               const OCP_DBL& Pbb,
                               const OCP_DBL& Tin,
                               const OCP_DBL* Ziin,
                               const USI&     tarPhase)
{
    return OWM.CalRho(Pin, tarPhase);
}

void BOMixture_OW::SetupWellOpt(WellOpt&                  opt,
                                const vector<SolventINJ>& sols,
                                const OCP_DBL&            Psurf,
                                const OCP_DBL&            Tsurf)
{
    const USI wellType = opt.WellType();
    if (wellType == INJ) {
        const string fluidName = opt.InjFluidType();
        opt.SetInjFactor(1.0);

        if (fluidName == "WAT") {
            vector<OCP_DBL> tmpZi({0, 1});
            opt.SetInjZi(tmpZi);
            opt.SetInjProdPhase(WATER);
        } else {
            OCP_ABORT("WRONG Injecting Fluid!");
        }
    } else if (wellType == PROD) {
        vector<OCP_DBL> tmpWght(2, 0);
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
                tmpWght[0] = tmpWght[1] = 1;
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
/*----------------------------------------------------------------------------*/