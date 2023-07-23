/*! \file    MixtureBlkOilUnit.cpp
 *  \brief   MixtureBlkOilUnit class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "MixtureUnitBlkOil.hpp"

 ///////////////////////////////////////////////
 // MixtureUnitBlkOil_OW
 ///////////////////////////////////////////////

MixtureUnitBlkOil_OW::MixtureUnitBlkOil_OW(const ParamReservoir& rs_param, const USI& i)
{
    mixtureType = BLKOIL_OW;
    numPhase = 2;
    numCom = 2;

    OWM.Setup(rs_param, i);
}

void MixtureUnitBlkOil_OW::Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{
    OWM.Flash(Pin, Niin);
}

void MixtureUnitBlkOil_OW::InitFlashIMPEC(const OCP_DBL& Pin,
    const OCP_DBL& Pbbin,
    const OCP_DBL& Tin,
    const OCP_DBL* Sjin,
    const OCP_DBL& Vpore,
    const OCP_DBL* Ziin,
    const OCP_USI& bId)
{
    OWM.InitFlash(Pin, Sjin[1], Vpore);
}

void MixtureUnitBlkOil_OW::InitFlashFIM(const OCP_DBL& Pin,
    const OCP_DBL& Pbbin,
    const OCP_DBL& Tin,
    const OCP_DBL* Sjin,
    const OCP_DBL& Vpore,
    const OCP_DBL* Ziin,
    const OCP_USI& bId)
{
    OWM.InitFlashDer(Pin, Sjin[1], Vpore);
}

void MixtureUnitBlkOil_OW::FlashIMPEC(const OCP_DBL& Pin,
    const OCP_DBL& Tin,
    const OCP_DBL* Niin,
    const USI& lastNP,
    const OCP_DBL* xijin,
    const OCP_USI& bId)
{
    OWM.Flash(Pin, Niin);
}

void MixtureUnitBlkOil_OW::FlashFIM(const OCP_DBL& Pin,
    const OCP_DBL& Tin,
    const OCP_DBL* Niin,
    const OCP_DBL* Sjin,
    const USI& lastNP,
    const OCP_DBL* xijin,
    const OCP_USI& bId)
{
    OWM.FlashDer(Pin, Niin);
}

OCP_DBL MixtureUnitBlkOil_OW::XiPhase(const OCP_DBL& Pin,
    const OCP_DBL& Tin,
    const OCP_DBL* Ziin,
    const USI& tarPhase)
{
    return OWM.CalXi(Pin, tarPhase);
}

OCP_DBL MixtureUnitBlkOil_OW::RhoPhase(const OCP_DBL& Pin,
    const OCP_DBL& Pbb,
    const OCP_DBL& Tin,
    const OCP_DBL* Ziin,
    const USI& tarPhase)
{
    return OWM.CalRho(Pin, tarPhase);
}

void MixtureUnitBlkOil_OW::SetupWellOpt(WellOpt& opt,
    const vector<SolventINJ>& sols,
    const OCP_DBL& Psurf,
    const OCP_DBL& Tsurf)
{
    const USI wellType = opt.WellType();
    if (wellType == INJ) {
        const string fluidName = opt.InjFluidType();
        opt.SetInjFactor(1.0);

        if (fluidName == "WAT") {
            vector<OCP_DBL> tmpZi({ 0, 1 });
            opt.SetInjZi(tmpZi);
            opt.SetInjProdPhase(WATER);
        }
        else {
            OCP_ABORT("WRONG Injecting Fluid!");
        }
    }
    else if (wellType == PROD) {
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
    }
    else {
        OCP_ABORT("Wrong Well Type!");
    }
}


///////////////////////////////////////////////
// MixtureUnitBlkOil_OGW
///////////////////////////////////////////////

MixtureUnitBlkOil_OGW::MixtureUnitBlkOil_OGW(const ParamReservoir& rs_param, const USI& i)
{
    mixtureType = BLKOIL_ODGW;
    numPhase = 3;
    numCom = 3;

    OGWM.Setup(rs_param, i);
}

void MixtureUnitBlkOil_OGW::Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{

    OGWM.Flash(Pin, Niin);
}


void MixtureUnitBlkOil_OGW::InitFlashIMPEC(const OCP_DBL& Pin,
    const OCP_DBL& Pbbin,
    const OCP_DBL& Tin,
    const OCP_DBL* Sjin,
    const OCP_DBL& Vpore,
    const OCP_DBL* Ziin,
    const OCP_USI& bId)
{
    OGWM.InitFlash(Pin, Pbbin, Sjin[1], Sjin[2], Vpore);
}

void MixtureUnitBlkOil_OGW::InitFlashFIM(const OCP_DBL& Pin,
    const OCP_DBL& Pbbin,
    const OCP_DBL& Tin,
    const OCP_DBL* Sjin,
    const OCP_DBL& Vpore,
    const OCP_DBL* Ziin,
    const OCP_USI& bId)
{
    OGWM.InitFlashDer(Pin, Pbbin, Sjin[1], Sjin[2], Vpore);
}

void MixtureUnitBlkOil_OGW::FlashIMPEC(const OCP_DBL& Pin,
    const OCP_DBL& Tin,
    const OCP_DBL* Niin,
    const USI& lastNP,
    const OCP_DBL* xijin,
    const OCP_USI& bId)
{
    OGWM.Flash(Pin, Niin);
}

void MixtureUnitBlkOil_OGW::FlashFIM(const OCP_DBL& Pin,
    const OCP_DBL& Tin,
    const OCP_DBL* Niin,
    const OCP_DBL* Sjin,
    const USI& lastNP,
    const OCP_DBL* xijin,
    const OCP_USI& bId)
{
    OGWM.FlashDer(Pin, Niin);
}

OCP_DBL
MixtureUnitBlkOil_OGW::XiPhase(const OCP_DBL& Pin,
    const OCP_DBL& Tin,
    const OCP_DBL* Ziin,
    const USI& tarPhase)
{
    return OGWM.CalXi(Pin, Pin, tarPhase);
}

OCP_DBL
MixtureUnitBlkOil_OGW::RhoPhase(const OCP_DBL& Pin,
    const OCP_DBL& Pbbin,
    const OCP_DBL& Tin,
    const OCP_DBL* Ziin,
    const USI& tarPhase)
{
    return OGWM.CalRho(Pin, Pbbin, tarPhase);
}

void MixtureUnitBlkOil_OGW::SetupWellOpt(WellOpt& opt,
    const vector<SolventINJ>& sols,
    const OCP_DBL& Psurf,
    const OCP_DBL& Tsurf)
{
    const USI wellType = opt.WellType();
    if (wellType == INJ) {
        const string fluidName = opt.InjFluidType();
        opt.SetInjFactor(1.0);

        if (fluidName == "WAT") {
            vector<OCP_DBL> tmpZi({ 0, 0, 1 });
            opt.SetInjZi(tmpZi);
            opt.SetInjProdPhase(WATER);
        }
        else if (fluidName == "GAS") {
            vector<OCP_DBL> tmpZi({ 0, 1, 0 });
            opt.SetInjZi(tmpZi);
            opt.SetInjProdPhase(GAS);
        }
        else {
            OCP_ABORT("WRONG Injecting Fluid!");
        }

    }
    else if (wellType == PROD) {
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
    }
    else {
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