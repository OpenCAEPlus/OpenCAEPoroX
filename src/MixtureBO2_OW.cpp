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
    BOMixtureInit(rs_param);

    PVTW.Setup(rs_param.PVTW_T.data[i], std_RhoW);
    PVDO.Setup(rs_param.PVDO_T.data[i], std_RhoO); 

    phaseExist[0] = OCP_TRUE;
    phaseExist[1] = OCP_TRUE;

    xij[0 * 2 + 0] = 1;
    xij[0 * 2 + 1] = 0;
    xij[1 * 2 + 0] = 0;
    xij[1 * 2 + 1] = 1;
}

void BOMixture_OW::Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{
    FlashIMPEC(Pin, Tin, Niin, 0, 0, 0);
}

void BOMixture_OW::InitFlashIMPEC(const OCP_DBL& Pin,
                                  const OCP_DBL& Pbbin,
                                  const OCP_DBL& Tin,
                                  const OCP_DBL* Sjin,
                                  const OCP_DBL& Vpore,
                                  const OCP_DBL* Ziin,
                                  const OCP_USI& bId)
{
    P    = Pin;
    S[1] = Sjin[1];
    // Water Properties
    PVTW.CalRhoXiMuDer(P, rho[1], xi[1], mu[1], rhoP[1], xiP[1], muP[1]);

    Ni[1] = Vpore * S[1] * xi[1];

    // Oil Properties
    PVDO.CalRhoXiMuDer(P, rho[0], xi[0], mu[0], rhoP[0], xiP[0], muP[0]);

    Ni[0]  = Vpore * (1 - S[1]) * xi[0];

    vj[0]  = Ni[0] / xi[0];
    vj[1]  = Ni[1] / xi[1];
    vf     = vj[0] + vj[1];
    S[0]   = vj[0] / vf;
    S[1]   = vj[1] / vf;
    vfi[0] = 1 / xi[0];
    vfi[1] = 1 / xi[1];
    vfP    = -Ni[0] * xiP[0] / (xi[0] * xi[0]) - Ni[1] * xiP[1] / (xi[1] * xi[1]);
}

void BOMixture_OW::InitFlashFIM(const OCP_DBL& Pin,
                                const OCP_DBL& Pbbin,
                                const OCP_DBL& Tin,
                                const OCP_DBL* Sjin,
                                const OCP_DBL& Vpore,
                                const OCP_DBL* Ziin,
                                const OCP_USI& bId)
{
    P       = Pin;
    S[oil]  = Sjin[oil];
    // Water Properties
    xi[wat] = PVTW.CalXiW(P);
    Ni[WAT] = Vpore * S[wat] * xi[wat];
    // Oil Properties
    xi[oil] = PVDO.CalXiO(P);
    Ni[OIL] = Vpore * (1 - S[wat]) * xi[oil];

    FlashFIM(Pin, Tin, &Ni[0], 0, 0, 0, 0);
}

void BOMixture_OW::FlashIMPEC(const OCP_DBL& Pin,
                              const OCP_DBL& Tin,
                              const OCP_DBL* Niin,
                              const USI&     lastNP,
                              const OCP_DBL* xijin,
                              const OCP_USI& bId)
{

    P     = Pin;
    Ni[0] = Niin[0];
    Ni[1] = Niin[1];

    // Water Properties
    PVTW.CalRhoXiMuDer(P, rho[1], xi[1], mu[1], rhoP[1], xiP[1], muP[1]);

    // Oil Properties
    PVDO.CalRhoXiMuDer(P, rho[0], xi[0], mu[0], rhoP[0], xiP[0], muP[0]);

    vj[0]  = Ni[0] / xi[0];
    vj[1]  = Ni[1] / xi[1];
    vf     = vj[0] + vj[1];
    S[0]   = vj[0] / vf;
    S[1]   = vj[1] / vf;
    vfi[0] = 1 / xi[0];
    vfi[1] = 1 / xi[1];
    vfP    = -Ni[0] * xiP[0] / (xi[0] * xi[0]) - Ni[1] * xiP[1] / (xi[1] * xi[1]);
}

void BOMixture_OW::FlashFIM(const OCP_DBL& Pin,
                            const OCP_DBL& Tin,
                            const OCP_DBL* Niin,
                            const OCP_DBL* Sjin,
                            const USI&     lastNP,
                            const OCP_DBL* xijin,
                            const OCP_USI& bId)
{
    fill(dXsdXp.begin(), dXsdXp.end(), 0.0);

    P     = Pin;
    Ni[OIL] = Niin[OIL];
    Ni[WAT] = Niin[WAT];
    Nt      = Ni[OIL] + Ni[WAT];

    // Water Properties
    PVTW.CalRhoXiMuDer(P, rho[wat], xi[wat], mu[wat], rhoP[wat], xiP[wat], muP[wat]);

    // Oil Properties
    PVDO.CalRhoXiMuDer(P, rho[oil], xi[oil], mu[oil], rhoP[oil], xiP[oil], muP[oil]);

    vj[oil]       = Ni[OIL] / xi[oil];
    vj[wat]       = Ni[WAT] / xi[wat];
    vf            = vj[oil] + vj[wat];
    S[oil]        = vj[oil] / vf;
    S[wat]        = vj[wat] / vf;
    vji[oil][OIL] = 1 / xi[oil];
    vji[wat][WAT] = 1 / xi[wat];
    vjP[oil]      = -Ni[OIL] * xiP[oil] / (xi[oil] * xi[oil]);
    vjP[wat]      = -Ni[WAT] * xiP[wat] / (xi[wat] * xi[wat]);




    vfi[OIL]      = vji[oil][OIL];
    vfi[WAT]      = vji[wat][WAT];
    vfP           = vjP[oil] + vjP[wat];

    dXsdXp[0] = (vjP[oil] - S[oil] * vfP) / vf;           // dSo / dP
    dXsdXp[1] = (vji[oil][OIL] - S[oil] * vfi[OIL]) / vf; // dSo / dNo
    dXsdXp[2] = -S[oil] * vfi[WAT] / vf;                  // dSo / dNw

    dXsdXp[3] = (vjP[wat] - S[wat] * vfP) / vf;           // dSw / dP
    dXsdXp[4] = -S[wat] * vfi[OIL] / vf;                  // dSw / dNo
    dXsdXp[5] = (vji[wat][WAT] - S[wat] * vfi[WAT]) / vf; // dSw / dNw
}

OCP_DBL BOMixture_OW::XiPhase(const OCP_DBL& Pin,
                              const OCP_DBL& Tin,
                              const OCP_DBL* Ziin,
                              const USI&     tarPhase)
{
    if (tarPhase == WATER) {
        // inj fluid is water
        return PVTW.CalXiW(Pin);
    } else {
        OCP_ABORT("Wrong tarPhase!");
    }
}

OCP_DBL BOMixture_OW::RhoPhase(const OCP_DBL& Pin,
                               const OCP_DBL& Pbb,
                               const OCP_DBL& Tin,
                               const OCP_DBL* Ziin,
                               const USI&     tarPhase)
{
    if (tarPhase == OIL) {
        return PVDO.CalRhoO(Pin);
    } else if (tarPhase == WATER) {
        return PVTW.CalRhoW(Pin);
    } else {
        OCP_ABORT("Wrong tarPhase!");
    }
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