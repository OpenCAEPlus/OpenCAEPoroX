/*! \file    MixtureComp.cpp
 *  \brief   MixtureComp class definition for compositional models
 *  \author  Shizhe Li
 *  \date    Jan/05/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "MixtureComp.hpp"


MixtureComp::MixtureComp(const ComponentParam& param, const USI& tarId)
{
    mixtureType = EOS_PVTW;
    // for Mixture class
    // Now, only one case is considered: oil, gas, water could exist
    numPhase = param.numPhase + 1;
    numCom   = param.numCom + 1;
    NPmax    = param.numPhase;
    NC       = param.numCom;

    if (param.Tc.activity)   Tc = param.Tc.data[tarId];
    else                     OCP_ABORT("TCRIT hasn't been input!");
    if (param.Pc.activity)   Pc = param.Pc.data[tarId];
    else                     OCP_ABORT("PCRIT hasn't been input!");

    if (param.Vc.activity)   Vc = param.Vc.data[tarId];
    else if (param.Zc.activity) {
        const vector<OCP_DBL>& Zc = param.Zc.data[tarId];
        Vc.resize(NC);
        for (USI i = 0; i < NC; i++) {
            Vc[i] = GAS_CONSTANT * Zc[i] * Tc[i] / Pc[i];
        }
    }
    else                     OCP_ABORT("VCRIT or ZCRIT hasn't been input!");

    if (param.MW.activity)   MWC = param.MW.data[tarId];
    else                     OCP_ABORT("MW hasn't been input!");


    MW.resize(NPmax);
    zi.resize(NC);
    vm.resize(NPmax);
    vmP.resize(NPmax);
    vmx.resize(NPmax);
    for (auto& v : vmx) v.resize(NC);
    Allocate();
    vjp.resize(numPhase, 0);
    vji.resize(numPhase);
    for (auto& v : vji)  v.resize(numCom, 0);

    phaseLabel.resize(NPmax);
    epIndexH.resize(NPmax);
    epIndex.resize(numPhase);


    AllocateOthers();
    rowork.resize(NC);
    eos.Setup(param, tarId);
    visCal.Setup(param, tarId);
    PE.Setup(param, tarId, &eos);
}

void MixtureComp::Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{

    InitPTN(Pin, Tin + CONV5, Niin);
    PE.PhaseEquilibrium(P, T, &zi[0], 0, 0, 0);
    CalProperty();

    // Water Properties
    const USI Wpid            = numPhase - 1;
    const USI Wcid            = numCom - 1;
    phaseExist[Wpid]          = OCP_TRUE;
    x[Wpid * numCom + Wcid]   = 1.0;
    Nt                        = Nh + Ni[Wcid];

    PVTW.CalRhoXiMuDer(P, rho[Wpid], xi[Wpid], mu[Wpid], rhoP[Wpid], xiP[Wpid], muP[Wpid]);
    vj[Wpid]    = Ni[Wcid] / xi[Wpid];
    vf = 0;
    for (USI j = 0; j < NP; j++) {
        if (phaseExist[j]) {
            vf += vj[j];
        }
    }

    // Calculate Sj
    CalSaturation();
}


void MixtureComp::InitFlashIMPEC(const OCP_DBL& Pin,
                                 const OCP_DBL& Pbbin,
                                 const OCP_DBL& Tin,
                                 const OCP_DBL* Sjin,
                                 const OCP_DBL& Vpore,
                                 const OCP_DBL* Ziin,
                                 const OCP_USI& bId)
{
    // Attention: zi[numCom - 1] = 0 here, that's Zw = 0;
    SetBulkId(bId);
    ftype = 0;
    lNP   = 0;

    InitPTZ(Pin, Tin + CONV5, Ziin);
    PE.PhaseEquilibrium(P, T, &zi[0], 0, 0, 0);
    Nh = 1;
    CalProperty();

    CalSurfaceTension();

    // Calulate Nt, water is exclued
    const OCP_DBL Sw = Sjin[numPhase - 1];
    Nh               = Vpore * (1 - Sw) / (vj[0] + vj[1]);
    vj[0]            *= Nh;
    vj[1]            *= Nh;
    // CalVfiVfp_full01();
    CalVfiVfp_full02();
    // Calculate Ni
    for (USI i = 0; i < NC; i++) {
        Ni[i] = zi[i] * Nh;
    }

    // Water Properties
    USI Wpid                  = numPhase - 1;
    USI Wcid                  = numCom - 1;
    phaseExist[Wpid]          = OCP_TRUE;
    x[Wpid * numCom + Wcid] = 1.0;
    nj[Wpid] = Ni[Wcid];

    PVTW.CalRhoXiMuDer(P, rho[Wpid], xi[Wpid], mu[Wpid], rhoP[Wpid], xiP[Wpid], muP[Wpid]);

    Ni[Wcid]    = Vpore * Sw * xi[Wpid];
    Nt          = Nh + Ni[Wcid];
    vj[Wpid]    = Ni[Wcid] / xi[Wpid];
    vf = 0;
    for (USI j = 0; j < numPhase; j++) {
        if (phaseExist[j]) {
            vf += vj[j];
        }
    }
    vfi[Wcid]   = 1 / xi[Wpid];
    vfP         += -Ni[Wcid] * xiP[Wpid] / (xi[Wpid] * xi[Wpid]);

    // Calculate Sj
    CalSaturation();

    CalSkipForNextStep();
}

void MixtureComp::InitFlashFIM(const OCP_DBL& Pin,
                               const OCP_DBL& Pbbin,
                               const OCP_DBL& Tin,
                               const OCP_DBL* Sjin,
                               const OCP_DBL& Vpore,
                               const OCP_DBL* Ziin,
                               const OCP_USI& bId)
{
    // Attention: zi[numCom - 1] = 0 here, that's Zw = 0;
    SetBulkId(bId);
    ftype = 0;
    lNP   = 0;

    InitPTZ(Pin, Tin + CONV5, Ziin);
    PE.PhaseEquilibrium(P, T, &zi[0], 0, 0, 0);
    Nh = 1;
    CalPropertyDer();
    CalSurfaceTension();


    // Calulate Nt, water is exclued
    const OCP_DBL Sw = Sjin[numPhase - 1];
    Nh    = Vpore * (1 - Sw) / (vj[0] + vj[1]);
    nj[0] *= Nh;
    nj[1] *= Nh;
    vj[0] *= Nh;
    vj[1] *= Nh;
    // CalVfiVfp_full01();
    CalVfiVfp_full02();
    // Calculate Ni
    for (USI i = 0; i < NC; i++) {
        Ni[i] = zi[i] * Nh;
    }

    // Water Properties
    USI Wpid                  = numPhase - 1;
    USI Wcid                  = numCom - 1;
    phaseExist[Wpid]          = OCP_TRUE;
    x[Wpid * numCom + Wcid] = 1.0;

    PVTW.CalRhoXiMuDer(P, rho[Wpid], xi[Wpid], mu[Wpid], rhoP[Wpid], xiP[Wpid], muP[Wpid]);

    Ni[Wcid]          = Vpore * Sw * xi[Wpid];
    nj[Wpid]          = Ni[Wcid];
    Nt                = Nh + Ni[Wcid];
    vj[Wpid]          = Ni[Wcid] / xi[Wpid];
    vfi[Wcid]         = 1 / xi[Wpid];
    const OCP_DBL vwp = -Ni[Wcid] * xiP[Wpid] / (xi[Wpid] * xi[Wpid]);
    vf = 0;
    for (USI j = 0; j < numPhase; j++) {
        if (phaseExist[j]) {
            vf += vj[j];
        }
    }
    vfP += vwp;
    vji[numPhase - 1][numCom - 1] = vfi[Wcid];
    vjp[numPhase - 1]             = vwp;

    CalSaturation();

    CaldXsdXp02();
    CalSkipForNextStep();
}

void MixtureComp::FlashIMPEC(const OCP_DBL& Pin,
                             const OCP_DBL& Tin,
                             const OCP_DBL* Niin,
                             const USI&     lastNP,
                             const OCP_DBL* xijin,
                             const OCP_USI& bId)
{
    SetBulkId(bId);
    // Hydroncarbon phase, if lNp = 0, then strict stability analysis will be used
    lNP = lastNP > 0 ? lastNP - 1 : 0;
    vector<OCP_DBL> lKs(NC, 0);
    if (lNP == 2) {
        for (USI i = 0; i < NC; i++) {
            lKs[i] = xijin[i] / xijin[i + numCom];
        }
    }

    InitPTN(Pin, Tin + CONV5, Niin);
    CalFtypeIMPEC();
    PE.PhaseEquilibrium(P, T, &zi[0], ftype, lNP, &lKs[0]);
    CalProperty();
    
    // Calculate derivates for hydrocarbon phase and components
    // d vf / d Ni, d vf / d P
    // CalVfiVfp_full01();
    CalVfiVfp_full02();

    // Water Properties
    const USI Wpid            = numPhase - 1;
    const USI Wcid            = numCom - 1;
    phaseExist[Wpid]          = OCP_TRUE;
    x[Wpid * numCom + Wcid] = 1.0;
    nj[Wpid] = Ni[Wcid];
    Nt                        = Nh + Ni[Wcid];

    PVTW.CalRhoXiMuDer(P, rho[Wpid], xi[Wpid], mu[Wpid], rhoP[Wpid], xiP[Wpid], muP[Wpid]);

    vj[Wpid]    = Ni[Wcid] / xi[Wpid];
    vf = 0;
    for (USI j = 0; j < numPhase; j++) {
        if (phaseExist[j]) {
            vf += vj[j];
        }
    }
    vfi[Wcid]   = 1 / xi[Wpid];
    vfP        += -Ni[Wcid] * xiP[Wpid] / (xi[Wpid] * xi[Wpid]);

    // Calculate Sj
    CalSaturation();

    CalSkipForNextStep();
}

void MixtureComp::FlashFIM(const OCP_DBL& Pin,
                           const OCP_DBL& Tin,
                           const OCP_DBL* Niin,
                           const OCP_DBL* Sjin,
                           const USI&     lastNP,
                           const OCP_DBL* xijin,
                           const OCP_USI& bId)
{

    SetBulkId(bId);

    // Hydroncarbon phase, if lNp = 0, then strict stability analysis will be used
    vector<OCP_DBL> lKs(NC, 0);
    lNP = lastNP > 0 ? lastNP - 1 : 0;
    if (lNP == 2) {
        for (USI i = 0; i < NC; i++) {
            lKs[i] = xijin[i] / xijin[i + numCom];
        }
    }

    InitPTN(Pin, Tin + CONV5, Niin);
    CalFtypeFIM(Sjin);
    PE.PhaseEquilibrium(P, T, &zi[0], ftype, lNP, &lKs[0]);
    CalPropertyDer();

    // Calculate derivates for hydrocarbon phase and components
    // d vf / d Ni, d vf / d P 
    CalVfiVfp_full02();

    // Water Properties
    USI Wpid                  = numPhase - 1;
    USI Wcid                  = numCom - 1;
    phaseExist[Wpid]          = OCP_TRUE;
    x[Wpid * numCom + Wcid] = 1.0;
    nj[Wpid]                  = Ni[Wcid];
    Nt                        = Nh + Ni[Wcid];

    PVTW.CalRhoXiMuDer(P, rho[Wpid], xi[Wpid], mu[Wpid], rhoP[Wpid], xiP[Wpid], muP[Wpid]);
    vj[Wpid]          = Ni[Wcid] / xi[Wpid];
    vfi[Wcid]         = 1 / xi[Wpid];
    const OCP_DBL vwp = -Ni[Wcid] * xiP[Wpid] / (xi[Wpid] * xi[Wpid]);

    vf = 0;
    for (USI j = 0; j < numPhase; j++) {
        if (phaseExist[j]) {
            vf += vj[j];
        }
    }
    vfP              += vwp;
    vji[numPhase - 1][numCom - 1] = vfi[Wcid];
    vjp[numPhase - 1]             = vwp;

    CalSaturation();

    CaldXsdXp02();

    CalSkipForNextStep();
}


OCP_DBL
MixtureComp::XiPhase(const OCP_DBL& Pin,
                     const OCP_DBL& Tin,
                     const vector<OCP_DBL>& Ziin,
                     const USI&     tarPhase)
{
    // assume that only single phase exists here
    if (tarPhase == WATER) {
        // water phase
        return PVTW.CalXiW(Pin);
    } else {
        // oil phase
        return 1 / eos.CalVm(Pin, Tin + CONV5, &Ziin[0]);
    }
}

OCP_DBL
MixtureComp::RhoPhase(const OCP_DBL& Pin,
                      const OCP_DBL& Pbb,
                      const OCP_DBL& Tin,
                      const vector<OCP_DBL>& Ziin,
                      const USI&     tarPhase)
{

    // assume that only single phase exists here
    if (tarPhase == WATER) {
        // water phase
        return PVTW.CalRhoW(Pin);
    } else {
        // hydrocarbon phase
        OCP_DBL xitmp = XiPhase(Pin, Tin, Ziin, tarPhase);
        OCP_DBL MWtmp = 0;
        for (USI i = 0; i < NC; i++) MWtmp += zi[i] * MWC[i];
        return MWtmp * xitmp;
    }
}

void MixtureComp::SetupWellOpt(WellOpt&                  opt,
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
            opt.SetInjFactor(1.0);
        } else {
            // inj phase is gas
            opt.SetInjProdPhase(GAS);
            const USI len = sols.size();
            for (USI i = 0; i < len; i++) {
                if (fluidName == sols[i].name) {
                    tmpZi = sols[i].data;
                    tmpZi.resize(numCom);
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
        vector<OCP_DBL> tmpWght(numPhase, 0);
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

void MixtureComp::CalProdWeight(const OCP_DBL&         Pin,
                                const OCP_DBL&         Tin,
                                const OCP_DBL*         Niin,
                                const vector<OCP_DBL>& prodPhase,
                                vector<OCP_DBL>&       prodWeight)
{
    Flash(Pin, Tin, Niin);

    OCP_DBL         qt = Nt;
    vector<OCP_DBL> factor(numPhase, 0);

    factor[0] = vj[0] / qt / CONV1; // stb / lbmol
    factor[1] = vj[1] / qt / 1000;  // Mscf / lbmol
    factor[2] = xi[2] * vj[2] / qt; // stb  / stb

    OCP_DBL tmp = 0;
    for (USI i = 0; i < 3; i++) {
        tmp += factor[i] * prodPhase[i];
    }
    if (tmp < 1E-12 || !isfinite(tmp)) {
        OCP_ABORT("Wrong Condition!");
    }
    fill(prodWeight.begin(), prodWeight.end(), tmp);
}

void MixtureComp::CalProdRate(const OCP_DBL&   Pin,
                              const OCP_DBL&   Tin,
                              const OCP_DBL*   Niin,
                              vector<OCP_DBL>& prodRate)
{
    Flash(Pin, Tin, Niin);

    prodRate[0] = vj[0] / CONV1; // stb
    prodRate[1] = vj[1] / 1000;  // Mscf
    prodRate[2] = vj[2] * xi[2]; // stb
}


void MixtureComp::CalSaturation()
{
    for (USI j = 0; j < numPhase; j++) {
        S[j] = 0;
        if (phaseExist[j]) {
            S[j] = vj[j] / vf;
        }
    }
}


void MixtureComp::AllocateOthers()
{
    pivot.resize(NPmax * NC);
    lnfugP.resize(NPmax);
    lnfugN.resize(NPmax);
    for (USI j = 0; j < NPmax; j++) {
        lnfugP[j].resize(NC);
        lnfugN[j].resize(NC * NC);
    }
    JmatDer.resize(NPmax * NPmax * (NC + 1) * (NC + 1));
    JmatTmp = JmatDer;
    rhsDer.resize(NPmax * (NC + 1) * (NC + 1));

    // new
    JmatDer.resize((numPhase + NPmax * NC) * (numPhase + NPmax * NC));
    JmatTmp = JmatDer;
    rhsDer.resize((numPhase + NPmax * NC) * (numCom + 1 + 1));
}


void MixtureComp::CopyPhaseFromPE()
{
    NP       = PE.GetNP();
    ftype    = PE.GetFtype();
    flagSkip = PE.GetFlagSkip();
    for (USI j = 0; j < NP; j++) {
        nj[j] = Nh * PE.GetNu(j);
        copy(PE.GetX(j).begin(), PE.GetX(j).end(), &x[j * numCom]);
    }
}


void MixtureComp::CalMW()
{
    // Calculate Molecular Weight of phase
    for (USI j = 0; j < NP; j++) {
        MW[j] = 0;
        for (USI i = 0; i < NC; i++) {
            MW[j] += x[j * numCom + i] * MWC[i];
        }
    }
}


void MixtureComp::CalVmVj()
{
    for (USI j = 0; j < NP; j++) {
        vm[j] = eos.CalVmDer(P, T, &x[j * numCom], vmP[j], &vmx[j][0]);
        vj[j] = nj[j] * vm[j];
    }
}


void MixtureComp::CalProperty()
{
    CopyPhaseFromPE();
    CalMW();
    CalVmVj();
    CalXiRhoMu();
    IdentifyPhase();
    ReOrderPhase();
}


void MixtureComp::CalXiRhoMu()
{
    for (USI j = 0; j < NP; j++) {
        xi[j]  = 1 / vm[j];
        rho[j] = MW[j] * xi[j];
        mu[j]  = visCal.CalViscosity(ViscosityParams(&P, &T, &x[j + numCom], &xi[j]));
    }
}


void MixtureComp::CalPropertyDer()
{
    CopyPhaseFromPE();
    CalMW();
    CalVmVj();
    CalXiRhoMuDer();
    IdentifyPhase();
    ReOrderPhaseDer();
}


void MixtureComp::CalXiRhoMuDer()
{
    OCP_DBL   dummy;
    for (USI j = 0; j < NP; j++) {
        // molar density
        xi[j]  = 1 / vm[j];
        xiP[j] = -xi[j] * xi[j] * vmP[j];
        for (USI i = 0; i < NC; i++) {
            xix[j * numCom + i] = -xi[j] * xi[j] * vmx[j][i];
        }
        
        // mass density
        rho[j]  = MW[j] * xi[j];
        rhoP[j] = MW[j] * xiP[j];
        for (USI i = 0; i < NC; i++) {
            rhox[j * numCom + i] = xix[j * numCom + i] * MW[j] + xi[j] * MWC[i];
        }

        // viscosity
        mu[j] = visCal.CalViscosity(ViscosityParams(&P, &T, &x[j * numCom], &xi[j], &xiP[j], nullptr, &xix[j * numCom]),
            muP[j], dummy, &mux[j * numCom]);
    }
}


void MixtureComp::IdentifyPhase()
{
    phaseExist[0] = OCP_FALSE;
    phaseExist[1] = OCP_FALSE;
    if (NP == 1) {
        // Critical Temperature Method
        OCP_DBL A = 0;
        OCP_DBL B = 0;
        for (USI i = 0; i < NC; i++) {
            A += zi[i] * Vc[i] * Tc[i];
            B += zi[i] * Vc[i];
        }
        OCP_DBL Tc = A / B;
        if (T > Tc) {
            phaseLabel[0]   = GAS;
            phaseExist[GAS] = OCP_TRUE;
        }
        else {
            phaseLabel[0]   = OIL;
            phaseExist[OIL] = OCP_TRUE;
        }
    }
    else {
        // Compare MW
        if (MW[0] > MW[1]) {
            phaseLabel[0] = OIL;
            phaseLabel[1] = GAS;
        }
        else {
            phaseLabel[0] = GAS;
            phaseLabel[1] = OIL;
        }
        phaseExist[OIL] = OCP_TRUE;
        phaseExist[GAS] = OCP_TRUE;
    }

    epIndexH.clear();
    for (USI j = 0; j < numPhase - 1; j++) {
        if (phaseExist[j]) epIndexH.push_back(j);
    }
    epIndex = epIndexH;
    epIndex.push_back(WATER);
}


void MixtureComp::ReOrderPhase()
{
    // for NP <= 2 Now
    if (phaseLabel[0] != OIL) {
        OCPSwap(&nj[0], &nj[1], 1, &rowork[0]);
        OCPSwap(&x[0 * numCom], &x[1 * numCom], NC, &rowork[0]);
        OCPSwap(&MW[0], &MW[1], 1, &rowork[0]);
        OCPSwap(&vj[0], &vj[1], 1, &rowork[0]);
        OCPSwap(&vm[0], &vm[1], 1, &rowork[0]);
        OCPSwap(&vmP[0], &vmP[1], 1, &rowork[0]);
        OCPSwap(&vmx[0][0], &vmx[1][0], NC, &rowork[0]);
        OCPSwap(&xi[0], &xi[1], 1, &rowork[0]);
        OCPSwap(&rho[0], &rho[1], 1, &rowork[0]);
        OCPSwap(&mu[0], &mu[1], 1, &rowork[0]);
    }
}


void MixtureComp::ReOrderPhaseDer()
{
    // for NP <= 2 Now
    if (phaseLabel[0] != OIL) {
        OCPSwap(&nj[0], &nj[1], 1, &rowork[0]);
        OCPSwap(&x[0 * numCom], &x[1 * numCom], NC, &rowork[0]);
        OCPSwap(&MW[0], &MW[1], 1, &rowork[0]);
        OCPSwap(&vj[0], &vj[1], 1, &rowork[0]);
        OCPSwap(&vm[0], &vm[1], 1, &rowork[0]);
        OCPSwap(&vmP[0], &vmP[1], 1, &rowork[0]);
        OCPSwap(&vmx[0][0], &vmx[1][0], NC, &rowork[0]);
        OCPSwap(&xi[0], &xi[1], 1, &rowork[0]);
        OCPSwap(&xiP[0], &xiP[1], 1, &rowork[0]);
        OCPSwap(&xix[0 * numCom], &xix[1 * numCom], NC, &rowork[0]);
        OCPSwap(&rho[0], &rho[1], 1, &rowork[0]);
        OCPSwap(&rhoP[0], &rhoP[1], 1, &rowork[0]);
        OCPSwap(&rhox[0 * numCom], &rhox[1 * numCom], NC, &rowork[0]);
        OCPSwap(&mu[0], &mu[1], 1, &rowork[0]);
        OCPSwap(&muP[0], &muP[1], 1, &rowork[0]);
        OCPSwap(&mux[0 * numCom], &mux[1 * numCom], NC, &rowork[0]);
    }
}


void MixtureComp::CalVfiVfp_full01()
{

    if (NP == 1) {
        const USI j = epIndexH[0];
        vjp[j]      = nj[j] * vmP[j];
        vfP         = vjp[j];
        for (USI i = 0; i < NC; i++) {
            vji[j][i] = vm[j] + vmx[j][i];
            for (USI k = 0; k < NC; k++) {
                vji[j][i] -= vmx[j][k] * x[j * numCom + k];
            }
            vfi[i] = vji[j][i];
        }
    } else {
        // NP > 1
        for (const auto& j : epIndexH) {
            eos.CalLnFugN(P, T, &x[j * numCom], nj[j], &lnfugN[j][0]);
            eos.CalLnFugP(P, T, &x[j * numCom], &lnfugP[j][0]);
        }
        AssembleMatVfiVfp_full01();
        AssembleRhsVfiVfp_full01();
        LUSolve(NC + 1, NC * NP, &JmatDer[0], &rhsDer[0], &pivot[0]);

        // now d nm0 / dP(dNk) has been available
        const OCP_DBL* dnkjdP = &rhsDer[0];
        // Calculate Vfp
        vfP = 0;
        vector<OCP_DBL> tmp(NP, 0);
        for (const auto& j : epIndexH) {
            for (USI i = 0; i < NC; i++) {
                tmp[j] -= vmx[j][i] * x[j * numCom + i];
            }          
            vjp[j] = nj[j] * vmP[j];
            for (USI k = 0; k < NC; k++) {
                vjp[j] += (vm[j] + tmp[j] + vmx[j][k]) * dnkjdP[k];
            }
            dnkjdP += NC;
            vfP += vjp[j];
        }

        // Calculate Vfi
        const OCP_DBL* dnkjdN = dnkjdP;
        for (USI i = 0; i < NC; i++) {
            vfi[i] = 0;
            for (const auto& j : epIndexH) {
                vji[j][i] = 0;
                for (USI k = 0; k < NC; k++) {
                    vji[j][i] += (vm[j] + tmp[j] + vmx[j][k]) * dnkjdN[k];
                }
                vfi[i] += vji[j][i];
                dnkjdN += NC;
            }
        }
    }

#ifdef OCP_NANCHECK
    if (!CheckNan(vfi.size(), &vfi[0])) {
        OCP_ABORT("INF or NAN in vfi !");
    }
    if (!CheckNan(1, &vfP)) {
        OCP_ABORT("INF or NAN in vfP !");
    }
#endif // NANCHECK
}

void MixtureComp::AssembleMatVfiVfp_full01()
{
    fill(JmatDer.begin(), JmatDer.end(), 0.0);
    // Attention 1: JmatDer should be sorted by column
    // Attention 2: d ln fij / d nkj is symetric for each j;
    OCP_DBL* bId = &JmatDer[0];
    for (USI j = 0; j < NP - 1; j++) {
        // for jth phase
        OCP_DBL* fugNj = &lnfugN[epIndexH[j]][0];
        for (USI i = 0; i < NC; i++) {
            // for ith components
            Dcopy(NC, bId, fugNj);
            bId    += (NP - 1 - j) * NC;
            bId[i] = 1.0;
            bId    += (1 + j) * NC;
            fugNj  += NC;
        }
        bId += NC;
    }
    // NP - 1 phase
    bId            = &JmatDer[(NP - 1) * (NP * NC * NC)];
    OCP_DBL* fugNj = &lnfugN[epIndexH[NP - 1]][0];
    for (USI i = 0; i < NC; i++) {
        for (USI j = 0; j < NP - 1; j++) {
            Dcopy(NC, bId, fugNj);
            bId += NC;
        }
        Dscalar((NP - 1) * NC, -1.0, bId - (NP - 1) * NC);
        fugNj  += NC;
        bId[i] = 1.0;
        bId    += NC;
    }
}

void MixtureComp::AssembleRhsVfiVfp_full01()
{
    fill(rhsDer.begin(), rhsDer.end(), 0.0);
    OCP_DBL* rhstmp = &rhsDer[0];
    const USI refP  = epIndexH[NP - 1];
    // d P
    for (USI j = 0; j < NP - 1; j++) {
        for (USI i = 0; i < NC; i++) {
            rhstmp[i] = lnfugP[refP][i] - lnfugP[epIndexH[j]][i];
        }
        rhstmp += NC;
    }
    // d Nk
    for (USI k = 0; k < NC; k++) {       
        rhstmp[NC * (NP - 1) + k] = 1;
        rhstmp += NP * NC;
    }

#ifdef OCP_NANCHECK
    if (!CheckNan(rhsDer.size(), &rhsDer[0])) {
        OCP_ABORT("INF or NAN in rhsDer !");
    }
#endif // NANCHECK
}

void MixtureComp::CaldXsdXp01()
{
    // dS / dP
    // S = Sj, xij
    // P = P, Ni
    // water is included
    fill(dXsdXp.begin(), dXsdXp.end(), 0);
    const USI     ncol = numCom + 1;
    const OCP_DBL vf2  = vf * vf;

    // dS / dP, dS / dN
    for (const auto& j : epIndex) {
		OCP_DBL* bId = &dXsdXp[j * ncol];
		// dS / dP
		bId[0] = (vjp[j] * vf - vfP * vj[j]) / vf2;
		bId++;
		// dS / dN
		for (USI m = 0; m < numCom; m++) {
			bId[m] = (vji[j][m] * vf - vfi[m] * vj[j]) / vf2;
		}
    }

    if (NP == 1) {
        // dxij / dNm
        OCP_DBL* bId = &dXsdXp[(numPhase + epIndexH[0] * numCom) * ncol + 1];
        for (USI i = 0; i < NC; i++) {
            for (USI m = 0; m < NC; m++) {
                bId[m] = (delta(i, m) * Nh - Ni[i]) / (Nh * Nh);
            }
            bId += ncol;
        }
    } 
    else {
        // NP >= 2
        OCP_DBL* const bId = &dXsdXp[numPhase * ncol];
        // water component is only in water phase in current case
        // dxij / dP
        const OCP_DBL* dnkjdP = &rhsDer[0];
        for (const auto& j : epIndexH) {
            OCP_DBL*  bId0     = bId + j * numCom * ncol;                      
            OCP_DBL   njDerSum = 0;
            for (USI i = 0; i < NC; i++) njDerSum += dnkjdP[i];
            for (USI i = 0; i < NC; i++) {
                bId0[0] = (dnkjdP[i]  - njDerSum * x[j * numCom + i]) / nj[j];
                bId0    += ncol;
            }
            dnkjdP += NC;
        }
        // dxij / dNm
        const OCP_DBL* dnkjdN = dnkjdP;
        for (USI m = 0; m < NC; m++) {
            for (const auto& j : epIndexH) {
                OCP_DBL   njDerSum = 0;
                for (USI i = 0; i < NC; i++) njDerSum += dnkjdN[i];
                OCP_DBL* bId0 = bId + j * numCom * ncol + m + 1;
                for (USI i = 0; i < NC; i++) {
                    bId0[0] = (dnkjdN[i] - njDerSum * x[j * numCom + i]) / nj[j];
                    bId0    += ncol;
                }
                dnkjdN += NC;
            }
        }
    }
}


void MixtureComp::CalVfiVfp_full02()
{
    // Attention!
    // NP = 1 or NP = 2

    if (NP == 1) {
        const USI j = epIndexH[0];
		vjp[j]      = nj[j] * vmP[j];
		vfP         = vjp[j];
		for (USI i = 0; i < NC; i++) {
			vji[j][i] = vm[j] + vmx[j][i];
			for (USI k = 0; k < NC; k++) {
				vji[j][i] -= vmx[j][k] * x[j * numCom + k];
			}
			vfi[i] = vji[j][i];
		}
    } 
    else if (NP == 2) {
        // JUST FOR NP = 2
        for (const auto& j : epIndexH) {
            eos.CalLnFugN(P, T, &x[j * numCom], nj[j], &lnfugN[j][0]);
            eos.CalLnFugP(P, T, &x[j * numCom], &lnfugP[j][0]);
        }

        AssembleMatVfiVfp_full02();
        AssembleRhsVfiVfp_full02();
        LUSolve(NC + 1, NC, &JmatDer[0], &rhsDer[0], &pivot[0]);

        // now d nm0 / dP(dNk) has been available
        OCP_DBL tmp0 = 0;
        OCP_DBL tmp1 = 0;
        for (USI i = 0; i < NC; i++) {
            tmp0 -= vmx[0][i] * x[0 * numCom + i];
            tmp1 -= vmx[1][i] * x[1 * numCom + i];
        }

        // vfP
        const OCP_DBL* dnkjdP = &rhsDer[0];
        vjp[0] = nj[0] * vmP[0];
        vjp[1] = nj[1] * vmP[1];
        for (USI k = 0; k < NC; k++) {
            vjp[0] += (vm[0] + tmp0 + vmx[0][k]) * dnkjdP[k];
            vjp[1] -= (vm[1] + tmp1 + vmx[1][k]) * dnkjdP[k];
        }
        vfP = vjp[0] + vjp[1];

        // Calculate Vfi
        const OCP_DBL* dnkjdN = dnkjdP + NC;
        for (USI i = 0; i < NC; i++) {
            vji[0][i] = 0;
            vji[1][i] = 0;
            for (USI k = 0; k < NC; k++) {
                vji[0][i] += (vm[0] + tmp0 + vmx[0][k]) * dnkjdN[k];
                vji[1][i] += (vm[1] + tmp1 + vmx[1][k]) * (delta(i, k) - dnkjdN[k]);
            }
            vfi[i] = vji[0][i] + vji[1][i];
            dnkjdN += NC;
        }
    }
    else {
        OCP_ABORT("USE CalVfiVfp_full01 !");
    }
}

void MixtureComp::AssembleMatVfiVfp_full02()
{
    // NP = 2
    fill(JmatDer.begin(), JmatDer.end(), 0.0);
    // Attention 1: JmatDer should be sorted by column
    // Attention 2: d ln fij / d nkj is symetric for each j;
    OCP_DBL* tmpMat = &JmatDer[0];
    for (USI i = 0; i < NC; i++) {
        for (USI m = 0; m < NC; m++) {
            tmpMat[m] = lnfugN[0][i * NC + m] + lnfugN[1][i * NC + m];
        }
        tmpMat += NC;
    }
}

void MixtureComp::AssembleRhsVfiVfp_full02()
{
    // NP = 2
    fill(rhsDer.begin(), rhsDer.end(), 0.0);
    OCP_DBL* rhstmp = &rhsDer[0];

    // dP
    for (USI i = 0; i < NC; i++) {
        rhstmp[i] = lnfugP[1][i] - lnfugP[0][i];
    }
    rhstmp += NC;

    // dNk
    for (USI k = 0; k < NC; k++) {
        for (USI i = 0; i < NC; i++) {
            rhstmp[i] = lnfugN[1][k * NC + i]; // d lnfij / d nkj = d lnfkj / d nij
        }
        rhstmp += NC;
    }
}


void MixtureComp::CaldXsdXp02()
{
    // Attention!
    // NP = 1 or NP = 2

    // dS / dP
    // S = Sj, xij
    // P = P, Ni
    // water is included
    fill(dXsdXp.begin(), dXsdXp.end(), 0);
    const USI     ncol = numCom + 1;
    const OCP_DBL vf2  = vf * vf;

    for (const auto& j : epIndex) {
        OCP_DBL* bId = &dXsdXp[j * ncol];
        // dS / dP
        bId[0] = (vjp[j] * vf - vfP * vj[j]) / vf2;
        bId++;
        // dS / dN
        for (USI m = 0; m < numCom; m++) {
            bId[m] = (vji[j][m] * vf - vfi[m] * vj[j]) / vf2;
        }
    }

    if (NP == 1) {
        // dxij / dNm
		OCP_DBL* bId = &dXsdXp[(numPhase + epIndexH[0] * numCom) * ncol + 1];
		for (USI i = 0; i < NC; i++) {
			for (USI m = 0; m < NC; m++) {
				bId[m] = (delta(i, m) * Nh - Ni[i]) / (Nh * Nh);
			}
			bId += ncol;
		}
    } 
    else if (NP == 2)  {
        // NP = 2
        OCP_DBL* bId            = &dXsdXp[numPhase * ncol];
        // dxij / dP, dxij / dNm
        const OCP_DBL* dnkjdNP  = &rhsDer[0];
        const OCP_DBL* x0       = &x[0 * numCom];
        const OCP_DBL* x1       = &x[1 * numCom];      
        OCP_DBL*       bId0     = bId + 0 * numCom * ncol;
        OCP_DBL*       bId1     = bId + 1 * numCom * ncol;
        OCP_DBL        njDerSum = 0;
        // dxij / dP
        for (USI i = 0; i < NC; i++) njDerSum += dnkjdNP[i];
        for (USI i = 0; i < NC; i++) {
            bId0[0] = (dnkjdNP[i] - njDerSum * x0[i]) / nj[0];
            bId1[0] = (-dnkjdNP[i] + njDerSum * x1[i]) / nj[1];
            bId0    += ncol;
            bId1    += ncol;
        }
        dnkjdNP += NC;
       
        // dxij / dNm
        for (USI m = 0; m < NC; m++) {
            njDerSum = 0;
            for (USI i = 0; i < NC; i++) njDerSum += dnkjdNP[i];
            bId0 = bId + 0 * numCom * ncol + m + 1;
            bId1 = bId + 1 * numCom * ncol + m + 1;
            for (USI i = 0; i < NC; i++) {
                bId0[0] = (dnkjdNP[i] - njDerSum * x0[i]) / nj[0];
                bId1[0] = ((delta(i, m) - dnkjdNP[i]) - (1 - njDerSum) * x1[i]) / nj[1];
                bId0    += ncol;
                bId1    += ncol;
            }
            dnkjdNP += NC;
        }
    }
    else {
        OCP_ABORT("USE CaldXsdXp01 !");
    }
}


/////////////////////////////////////////////////////////////////////
// Optional Features
/////////////////////////////////////////////////////////////////////

void MixtureComp::SetupOptionalFeatures(OptionalFeatures& optFeatures)
{
    skipSta = &optFeatures.skipStaAnaly;
    if (skipSta->IfUseSkip()) {
        AllocateSkip();
        skipSta->Setup(optFeatures.numBulk, numPhase - 1, numCom - 1);       
    }
    miscible = &optFeatures.miscible;
    mIndex   =  miscible->Setup(optFeatures.numBulk, 
                     SurTenMethodParams(stm01Params), 
                     MisFacMethodParams(mfm01Params));
}


void MixtureComp::AllocateSkip()
{
    lnphiN.resize(NC * NC);
    skipMatSTA.resize(NC * NC);
    eigenSkip.resize(NC);
    eigenWork.resize(2 * NC + 1);
}


void MixtureComp::CalPhiNSkip()
{
    eos.CalLnPhiN(P, T, &zi[0], Nh, &lnphiN[0]);
}


void MixtureComp::AssembleSkipMatSTA()
{
    // Sysmetric Matrix
    // stored by colum
    const vector<OCP_DBL>& xj = zi;

    for (USI i = 0; i < NC; i++) {
        for (USI j = 0; j <= i; j++) {
            skipMatSTA[i * NC + j] = delta(i, j) + Nh * sqrt(xj[i] * xj[j]) * lnphiN[i * NC + j];
            skipMatSTA[j * NC + i] = skipMatSTA[i * NC + j];
        }
    }
}

void MixtureComp::CalSkipForNextStep()
{
    if (flagSkip && ftype == 0) {
        // 1. Np == 1 is base for Skipping
        // 2. If flagSkip == true, then next stablity analysis is possible to be
        // skipped, it depends on if conditions are met
        // 3. If ftype == 0, then the range should be calculated, which also means last
        // skip is unsatisfied
        CalPhiNSkip();
        AssembleSkipMatSTA();
#ifdef DEBUG
        if (!CheckNan(skipMatSTA.size(), &skipMatSTA[0])) {
            OCP_WARNING("Nan in skipMatSTA!");
        }
#endif // DEBUG

        CalEigenSY(NC, &skipMatSTA[0], &eigenSkip[0], &eigenWork[0], 2 * NC + 1);
        skipSta->AssignValue(bulkId, eigenSkip[0], P, T, zi);
    }
    skipSta->SetFlagSkip(bulkId, flagSkip);
}

/////////////////////////////////////////////////////////////////////
// Miscible
/////////////////////////////////////////////////////////////////////

void MixtureComp::InputMiscibleParam(const ParamReservoir& rs_param, const USI& tarId)
{
    // Input parachor
    const Miscstr&        miscstr   = rs_param.miscstr;
    const ComponentParam& compParam = rs_param.comsParam;
    if (miscstr.ifMiscible) {
        if (compParam.parachor.activity) {
            if (compParam.parachor.data[tarId].size() == NC) {
                stm01Params = SurTenMethod01Params(compParam.parachor.data[tarId], &x[0], &x[numCom],
                                               &xi[0], &xi[1], &NP);
            }
            else {
                OCP_ABORT("PARACHOR has not been Input Correctly!");
            }
        }
        // Input surface tension term
        mfm01Params = MisFacMethod01Params(miscstr.surTenRef, miscstr.surTenExp, miscstr.surTenPc);
    }
}

void MixtureComp::CalSurfaceTension()
{
    miscible->CalMiscibleFactor(bulkId, mIndex);
}


/////////////////////////////////////////////////////////////////////
// For Output
/////////////////////////////////////////////////////////////////////

void MixtureComp::OutMixtureIters() const
{
    PE.OutMixtureIters();
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jan/05/2022      Create file                          */
/*----------------------------------------------------------------------------*/