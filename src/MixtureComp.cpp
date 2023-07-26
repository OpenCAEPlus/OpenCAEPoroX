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
    // if Water don't exist?
    // for Mixture class
    // Now, only one case is considered: oil, gas, water could exist
    numPhase = param.numPhase + 1;
    numCom   = param.numCom + 1;
    Allocate();
    vjp.resize(numPhase, 0);
    vji.resize(numPhase);
    for (auto& v : vji) {
        v.resize(numCom, 0);
    }

    // for MixtureComp class
    NC    = param.numCom;
    NPmax = param.numPhase;

    zi.resize(NC);

    Cname = param.Cname;
    if (param.Tc.activity)
        Tc = param.Tc.data[tarId];
    else
        OCP_ABORT("TCRIT hasn't been input!");
    if (param.Pc.activity)
        Pc = param.Pc.data[tarId];
    else
        OCP_ABORT("PCRIT hasn't been input!");

    if (param.Vc.activity)
        Vc = param.Vc.data[tarId];
    else if (param.Zc.activity) {
        const vector<OCP_DBL>& Zc = param.Zc.data[tarId];
        Vc.resize(NC);
        for (USI i = 0; i < NC; i++) {
            Vc[i] = 10.73159 * Zc[i] * Tc[i] / Pc[i];
        }
    } else
        OCP_ABORT("VCRIT or ZCRIT hasn't been input!");

    if (param.MW.activity)
        MWC = param.MW.data[tarId];
    else
        OCP_ABORT("MW hasn't been input!");
    if (param.Acf.activity)
        Acf = param.Acf.data[tarId];
    else
        OCP_ABORT("ACF hasn't been input!");

    if (param.Vcvis.activity)
        Vcvis = param.Vcvis.data[tarId];
    else if (param.Zcvis.activity) {
        Zcvis = param.Zcvis.data[tarId];
        Vcvis.resize(NC);
        for (USI i = 0; i < NC; i++) {
            Vcvis[i] = GAS_CONSTANT * Zcvis[i] * Tc[i] / Pc[i];
        }
    } else
        Vcvis = Vc;

    LBCcoef = param.LBCcoef;
    for (auto& lbc : LBCcoef) {
        lbc *= 10;
    }

    EoSctrl.SSMsta.maxIt = stoi(param.SSMparamSTA[0]);
    EoSctrl.SSMsta.tol   = stod(param.SSMparamSTA[1]);
    EoSctrl.SSMsta.eYt   = stod(param.SSMparamSTA[2]);
    EoSctrl.SSMsta.tol2  = EoSctrl.SSMsta.tol * EoSctrl.SSMsta.tol;

    EoSctrl.NRsta.maxIt = stoi(param.NRparamSTA[0]);
    EoSctrl.NRsta.tol   = stod(param.NRparamSTA[1]);
    EoSctrl.NRsta.tol2  = EoSctrl.NRsta.tol * EoSctrl.NRsta.tol;

    EoSctrl.SSMsp.maxIt = stoi(param.SSMparamSP[0]);
    EoSctrl.SSMsp.tol   = stod(param.SSMparamSP[1]);
    EoSctrl.SSMsp.tol2  = EoSctrl.SSMsp.tol * EoSctrl.SSMsp.tol;

    EoSctrl.NRsp.maxIt = stoi(param.NRparamSP[0]);
    EoSctrl.NRsp.tol   = stod(param.NRparamSP[1]);
    EoSctrl.NRsp.tol2  = EoSctrl.NRsp.tol * EoSctrl.NRsp.tol;

    EoSctrl.RR.maxIt = stoi(param.RRparam[0]);
    EoSctrl.RR.tol   = stod(param.RRparam[1]);
    EoSctrl.RR.tol2  = EoSctrl.RR.tol * EoSctrl.RR.tol;

    AllocatePhase();
    AllocateMethod();
    AllocateOthers();

    eos.Setup(param, tarId);
}

void MixtureComp::Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{
    ftype = 0;
    lNP   = 0;
    InitPTN(Pin, Tin + CONV5, Niin);
    CalFlash();

    // Water Properties
    const USI Wpid            = numPhase - 1;
    const USI Wcid            = numCom - 1;
    phaseExist[Wpid]          = OCP_TRUE;
    xij[Wpid * numCom + Wcid] = 1.0;
    Nt                        = Nh + Ni[Wcid];

    PVTW.CalRhoXiMuDer(P, rho[Wpid], xi[Wpid], mu[Wpid], rhoP[Wpid], xiP[Wpid], muP[Wpid]);
    vj[Wpid]    = Ni[Wcid] / xi[Wpid];
    vf += vj[Wpid];

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
    PhaseEquilibrium();
    // Attention Nt = 1 now
    CalMW();
    CalVfXiRho();
    CalViscosity();
    CalSurfaceTension();
    IdentifyPhase();
    CopyPhase();

    // Calulate Nt, water is exclued
    OCP_DBL Sw = Sjin[numPhase - 1];
    Nh         = Vpore * (1 - Sw) / vf;

    // Next, nu represents moles of phase instead of molar fraction of phase
    Dscalar(NP, Nh, &nu[0]);
    // correct vj, vf with new Nt
    Dscalar(NPmax, Nh, &vj[0]);
    vf *= Nh;
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
    xij[Wpid * numCom + Wcid] = 1.0;

    PVTW.CalRhoXiMuDer(P, rho[Wpid], xi[Wpid], mu[Wpid], rhoP[Wpid], xiP[Wpid], muP[Wpid]);

    Ni[Wcid]    = Vpore * Sw * xi[Wpid];
    Nt          = Nh + Ni[Wcid];
    vj[Wpid]    = Ni[Wcid] / xi[Wpid];
    vf         += vj[Wpid];
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
    PhaseEquilibrium();
    // Attention Nt = 1 now
    CalMW();
    CalVfXiRho();
    CalViscosity();
    CalSurfaceTension();
    IdentifyPhase();
    CopyPhase();

    // Calulate Nt, water is exclued
    OCP_DBL Sw = Sjin[numPhase - 1];
    Nh         = Vpore * (1 - Sw) / vf;

    // Next, nu represents moles of phase instead of molar fraction of phase
    Dscalar(NP, Nh, &nu[0]);
    // correct vj, vf with new Nt
    Dscalar(NPmax, Nh, &vj[0]);
    vf *= Nh;
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
    xij[Wpid * numCom + Wcid] = 1.0;

    PVTW.CalRhoXiMuDer(P, rho[Wpid], xi[Wpid], mu[Wpid], rhoP[Wpid], xiP[Wpid], muP[Wpid]);

    Ni[Wcid]          = Vpore * Sw * xi[Wpid];
    nj[Wpid]          = Ni[Wcid];
    Nt                = Nh + Ni[Wcid];
    vj[Wpid]          = Ni[Wcid] / xi[Wpid];
    vfi[Wcid]         = 1 / xi[Wpid];
    const OCP_DBL vwp = -Ni[Wcid] * xiP[Wpid] / (xi[Wpid] * xi[Wpid]);
    vf += vj[Wpid];
    vfP += vwp;
    vji[numPhase - 1][numCom - 1] = vfi[Wcid];
    vjp[numPhase - 1]             = vwp;

    CalSaturation();

    CaldXsdXp02();
    CalXiPNX_partial();
    CalRhoPX_partial();
    CalMuPX_partial();

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
    if (lNP == 2) {
        for (USI i = 0; i < NC; i++) {
            lKs[i] = xijin[i] / xijin[i + numCom];
        }
    }

    InitPTN(Pin, Tin + CONV5, Niin);
    CalFtypeIMPEC();
    CalFlash();
    // Calculate derivates for hydrocarbon phase and components
    // d vf / d Ni, d vf / d P
    // CalVfiVfp_full01();
    CalVfiVfp_full02();

    // Water Properties
    const USI Wpid            = numPhase - 1;
    const USI Wcid            = numCom - 1;
    phaseExist[Wpid]          = OCP_TRUE;
    xij[Wpid * numCom + Wcid] = 1.0;
    Nt                        = Nh + Ni[Wcid];

    PVTW.CalRhoXiMuDer(P, rho[Wpid], xi[Wpid], mu[Wpid], rhoP[Wpid], xiP[Wpid], muP[Wpid]);

    vj[Wpid]    = Ni[Wcid] / xi[Wpid];
    vf         += vj[Wpid];
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
    lNP = lastNP > 0 ? lastNP - 1 : 0;
    if (lNP == 2) {
        for (USI i = 0; i < NC; i++) {
            lKs[i] = xijin[i] / xijin[i + numCom];
        }
    }

    InitPTN(Pin, Tin + CONV5, Niin);
    CalFtypeFIM(Sjin);
    CalFlash();

    // Calculate derivates for hydrocarbon phase and components
    // d vf / d Ni, d vf / d P 
    CalVfiVfp_full02();

    // Water Properties
    USI Wpid                  = numPhase - 1;
    USI Wcid                  = numCom - 1;
    phaseExist[Wpid]          = OCP_TRUE;
    xij[Wpid * numCom + Wcid] = 1.0;
    nj[Wpid]                  = Ni[Wcid];
    Nt                        = Nh + Ni[Wcid];

    PVTW.CalRhoXiMuDer(P, rho[Wpid], xi[Wpid], mu[Wpid], rhoP[Wpid], xiP[Wpid], muP[Wpid]);
    vj[Wpid]          = Ni[Wcid] / xi[Wpid];
    vfi[Wcid]         = 1 / xi[Wpid];
    const OCP_DBL vwp = -Ni[Wcid] * xiP[Wpid] / (xi[Wpid] * xi[Wpid]);
    vf               += vj[Wpid];
    vfP              += vwp;
    vji[numPhase - 1][numCom - 1] = vfi[Wcid];
    vjp[numPhase - 1]             = vwp;

    CalSaturation();

    CaldXsdXp02();

    CalXiPNX_partial();
    CalRhoPX_partial();
    CalMuPX_partial();

    CalSkipForNextStep();
}


void MixtureComp::CalFlash()
{
    PhaseEquilibrium();
    // Next, nu represents moles of phase instead of molar fraction of phase
    Dscalar(NP, Nh, &nu[0]);
    CalMW();
    CalVfXiRho();
    CalViscosity();
    CalSurfaceTension();
    IdentifyPhase();
    CopyPhase();
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
        return 1 / eos.CalVm(Pin, Tin + CONV5, Ziin);
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
        NP            = 1;
        x[0]          = zi;
        CalMW();
        return MW[0] * xitmp;
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




void MixtureComp::AllocatePhase()
{
    // Allocate Memoery for Phase variables
    vC.resize(NPmax);
    nu.resize(NPmax);
    xiC.resize(NPmax);
    rhoC.resize(NPmax);
    MW.resize(NPmax);
    phaseLabel.resize(NPmax);
    x.resize(NPmax);
    phi.resize(NPmax);
    fug.resize(NPmax);
    n.resize(NPmax);
    vmj.resize(NPmax);
    vmP.resize(NPmax);
    vmx.resize(NPmax);
    for (auto& v : vmx) v.resize(numCom);
    for (USI j = 0; j < NPmax; j++) {
        x[j].resize(NC);
        phi[j].resize(NC);
        fug[j].resize(NC);
        n[j].resize(NC);
        vmx[j].resize(NC);
    }
    ln = n;
}


void MixtureComp::CalMW()
{
    // Calculate Molecular Weight of phase
    for (USI j = 0; j < NP; j++) {
        MW[j] = 0;
        for (USI i = 0; i < NC; i++) {
            MW[j] += x[j][i] * MWC[i];
        }
    }
}

void MixtureComp::CalVfXiRho()
{
    // Attention that nu is moles of phase instead of molar fraction of phase now
    vf = 0;
    for (USI j = 0; j < NP; j++) {
        vmj[j] = eos.CalVmDer(P, T, x[j], vmP[j], vmx[j]);
        vC[j] = vmj[j] * nu[j];
        vf += vC[j];
        xiC[j]  = 1 / vmj[j];
        rhoC[j] = MW[j] * xiC[j];
    }
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

USI MixtureComp::FindMWmax()
{
    // find the phase id with the highest Molecular Weight
    USI     tmpId = 0;
    OCP_DBL tmpMW = MW[0];
    for (USI j = 1; j < NP; j++) {
        if (tmpMW < MW[j]) {
            tmpMW = MW[j];
            tmpId = j;
        }
    }
    return tmpId;
}

void MixtureComp::x2n()
{
    // Total moles are supposed to be 1
    for (USI j = 0; j < NP; j++) {

        // nu[j] = fabs(nu[j]);
        for (USI i = 0; i < NC; i++) {
            n[j][i] = nu[j] * x[j][i];
        }
    }
}


void MixtureComp::AllocateMethod()
{

    Kw.resize(4);
    for (USI i = 0; i < 4; i++) {
        Kw[i].resize(NC);
    }
    Ks.resize(NPmax - 1);
    for (USI i = 0; i < NPmax - 1; i++) {
        Ks[i].resize(NC);
    }
    phiSta.resize(NC);
    fugSta.resize(NC);

    Y.resize(NC);
    di.resize(NC);
    resSTA.resize(NC);

    JmatSTA.resize(NC * NC);
    lKs.resize(NC);

    resRR.resize(NPmax - 1);
    resSP.resize(static_cast<size_t>(NC) * NPmax);
    JmatSP.resize(static_cast<size_t>(NC * NC) * NPmax * NPmax);
    fugX.resize(NPmax);
    fugN.resize(NPmax);
    for (USI j = 0; j < NPmax; j++) {
        fugX[j].resize(NC * NC);
        fugN[j].resize(NC * NC);
    }
    pivot.resize(static_cast<size_t>(NC + 1) * NPmax, 1);
    lJmatWork = NC * (NPmax - 1);
    JmatWork.resize(lJmatWork);

    pivot.resize(NPmax * static_cast<size_t>(NC) + numPhase, 1);
}

void MixtureComp::CalKwilson()
{
    for (USI i = 0; i < NC; i++) {
        Kw[0][i] = (Pc[i] / P) * exp(5.373 * (1 + Acf[i]) * (1 - Tc[i] / T));
        Kw[1][i] = 1 / Kw[0][i];
        Kw[2][i] = pow(Kw[0][i], 1.0 / 3);
        Kw[3][i] = pow(Kw[1][i], 1.0 / 3);
    }
}

void MixtureComp::PhaseEquilibrium()
{
    // Attention: sum of components' moles equals 1
    switch (ftype) {
        case 0:
            // flash from single phase
            flagSkip = OCP_TRUE;
            NP       = 1;
            nu[0]    = 1;
            x[0]     = zi;
            CalKwilson();
            while (!PhaseStable()) {
                NP++;
                PhaseSplit();
                if (NP == NPmax || NP == 1) break;
            }
            if (NP > 1) {
                flagSkip = OCP_FALSE;
            }
            // record error
            if (NP == 1) {
                if (EoSctrl.NRsta.conflag)
                    ePEC = EoSctrl.NRsta.realTol;
                else if (EoSctrl.SSMsta.conflag)
                    ePEC = EoSctrl.SSMsta.realTol;
                else
                    ePEC = 1E8;
            } else {
                if (EoSctrl.NRsp.conflag)
                    ePEC = EoSctrl.NRsp.realTol;
                else if (EoSctrl.SSMsp.conflag)
                    ePEC = EoSctrl.SSMsp.realTol;
                else
                    ePEC = 1E8;
            }

            break;
        case 1:
            // Skip Phase Stability analysis, only single phase exists
            flagSkip = OCP_TRUE;
            NP       = 1;
            nu[0]    = 1;
            x[0]     = zi;
            // record error
            ePEC = 0.0;
            break;

        case 2:
            // Skip Phase Stability analysis, two phases exist
            flagSkip = OCP_FALSE;
            NP       = 2;
            Yt       = 1.01;
            CalKwilson();
            PhaseSplit();

            if (EoSctrl.NRsp.conflag)
                ePEC = EoSctrl.NRsp.realTol;
            else if (EoSctrl.SSMsp.conflag)
                ePEC = EoSctrl.SSMsp.realTol;
            else
                ePEC = 1E8;
            break;

        default:
            OCP_ABORT("Wrong flash type!");
            break;
    }
}

OCP_BOOL MixtureComp::PhaseStable()
{
    if (NP == 1) {
        testPId = 0;
    } else {
        CalMW();
        testPId = FindMWmax();
    }

    EoSctrl.SSMsta.conflag = OCP_FALSE;
    EoSctrl.SSMsta.curIt   = 0;
    EoSctrl.NRsta.conflag  = OCP_FALSE;
    EoSctrl.NRsta.curIt    = 0;

    // Test if a phase is stable, if stable return OCP_TRUE, else return OCP_FALSE
    OCP_BOOL flag;
    USI      tmpNP = NP;

    if (lNP == 0) {
        // strict stability ananlysis
        flag = StableSSM(testPId);
    } else {
        flag = StableSSM01(testPId);
        if (!flag) tmpNP++;

        if (tmpNP != lNP) {
            flag     = StableSSM(testPId);
            flagSkip = OCP_FALSE;
        }
    }
    itersSSMSTA += EoSctrl.SSMsta.curIt;
    itersNRSTA += EoSctrl.NRsta.curIt;
    countsSSMSTA++;
    countsNRSTA++;
    return flag;
}

OCP_BOOL MixtureComp::StableSSM01(const USI& Id)
{
    // if unsatble, return OCP_FALSE
    // if stable, return OCP_TRUE

    OCP_DBL Stol  = EoSctrl.SSMsta.tol2;
    USI     maxIt = EoSctrl.SSMsta.maxIt;
    OCP_DBL eYt   = EoSctrl.SSMsta.eYt;
    OCP_DBL Ktol  = EoSctrl.SSMsta.Ktol;
    OCP_DBL dYtol = EoSctrl.SSMsta.dYtol;
    // OCP_DBL& Sk = EoSctrl.SSMsta.curSk;
    OCP_DBL Se, Sk, dY;

    OCP_BOOL flag, Tsol; // Tsol, trivial solution
    USI      iter, k;

    const vector<OCP_DBL>& xj = x[Id];

    eos.CalFug(P, T, x[Id], fug[Id]);

    const vector<OCP_DBL>& fugId = fug[Id];

    vector<OCP_DBL>& ks = Ks[0];
    for (k = 0; k < 2; k++) {

        ks   = Kw[k];
        iter = 0;
        flag = OCP_FALSE;
        Tsol = OCP_FALSE;
        while (OCP_TRUE) {
            Yt = 0;
            for (USI i = 0; i < NC; i++) {
                Y[i] = xj[i] * ks[i];
                Yt += Y[i];
            }
            Dscalar(NC, 1 / Yt, &Y[0]);

            if (iter > 0) {
                dY = 0;
                for (USI i = 0; i < NC; i++) {
                    dY += pow((Y[i] - di[i]), 2);
                }
                if (dY < dYtol) {
                    // converges
                    flag = OCP_TRUE;
                    break;
                }
            }
            eos.CalFug(P, T, Y, fugSta);
            Se = 0;
            Sk = 0;
            for (USI i = 0; i < NC; i++) {
                di[i] = fugId[i] / (fugSta[i] * Yt);
                ks[i] *= di[i];
                Se += pow(di[i] - 1, 2);
                // Sk += pow(ks[i] - 1, 2);
                Sk += pow(log(ks[i]), 2);
            }

            iter++;
            if (Se < Stol) {
                flag = OCP_TRUE;
                break;
            }
            if (Sk < Ktol) {
                // Sk < Ktol -> trivial solution
                flag = OCP_TRUE;
                Tsol = OCP_TRUE;
                break;
            }
            if (iter > maxIt) {
                break;
            }

            // Record last Y with di
            di = Y;
        }

        if (!Tsol) {
            // flag = StableNR(Id);
        }

        // cout << "Yt = " << setprecision(8) << scientific << Yt << "   " << setw(2)
        //     << "Sk = " << setprecision(3) << scientific << Sk << "   " << setw(2)
        //     << iter << "  ";
        EoSctrl.SSMsta.curIt += iter;

        if (flag && Yt > 1 - 0.1 && Sk > 1) {
            // close to phase boundary, or more than 1 phase, So don't skip at next step
            flagSkip = OCP_FALSE;
        }
        if (flag && Yt > 1 + eYt) {
            EoSctrl.SSMsta.conflag = OCP_TRUE;
            EoSctrl.SSMsta.realTol = sqrt(Se);
            return OCP_FALSE;
        }
    }
    EoSctrl.SSMsta.realTol = sqrt(Se);
    return OCP_TRUE;
}

OCP_BOOL MixtureComp::StableSSM(const USI& Id)
{
    // if unsatble, return OCP_FALSE
    // if stable, return OCP_TRUE

    const vector<OCP_DBL>& xj = x[Id];
    eos.CalFugPhi(P, T, x[Id], fug[Id], phi[Id]);
    const vector<OCP_DBL>& fugId = fug[Id];

    for (USI i = 0; i < NC; i++) {
        di[i] = phi[Id][i] * xj[i];
    }

    OCP_DBL  Stol  = EoSctrl.SSMsta.tol2;
    USI      maxIt = EoSctrl.SSMsta.maxIt;
    OCP_DBL  eYt   = EoSctrl.SSMsta.eYt;
    OCP_DBL  Se;
    OCP_BOOL flag;
    USI      iter;
    USI      k;

    // cout << "SSMBEGINS" << endl;

    for (k = 0; k < 4; k++) {

        Yt = 0;
        for (USI i = 0; i < NC; i++) {
            Y[i] = xj[i] / Kw[k][i];
            Yt += Y[i];
        }
        Dscalar(NC, 1 / Yt, &Y[0]);

        eos.CalFugPhi(P, T, Y, fugSta, phiSta);

        Se = 0;
        for (USI i = 0; i < NC; i++) {
            Se += pow(log(fugSta[i] / fugId[i] * Yt), 2);
        }

        flag = OCP_TRUE;
        iter = 0;

        // cout << "ssmbegins" << endl;

        while (Se > Stol) {

            // cout << setprecision(6) << scientific << Se;
            // cout << "         ";
            // cout << setprecision(12) << scientific << Yt;
            // cout << "         " << iter << endl;
            // PrintDX(NC, &xj[0]);
            // PrintDX(NC, &Y[0]);

            Yt = 0;
            for (USI i = 0; i < NC; i++) {
                Y[i] = di[i] / phiSta[i];
                Yt += Y[i];
            }
            Dscalar(NC, 1 / Yt, &Y[0]);

            eos.CalFugPhi(P, T, Y, fugSta, phiSta);

            Se = 0;
            for (USI i = 0; i < NC; i++) {
                Se += pow(log(fugSta[i] / fugId[i] * Yt), 2);
            }

            iter++;
            if (iter > maxIt) {
                flag = OCP_FALSE;
                break;
            }
        }
        EoSctrl.SSMsta.curIt += iter;
        flag = StableNR(Id);
        // here, a relaxation is needed, on the one hand it can prevent the influence
        // of rounding error, on the other hand, if Yt is too close to 1, phase
        // splitting calculation may get into trouble and single phase is indentified
        // finally
        if (flag && Yt > 1 + eYt) {
            // cout << "Yt = " << setprecision(12) << scientific << Yt << "    " <<
            // setw(3)
            //     << EoSctrl.SSMsta.curIt << "    " << setw(3) << EoSctrl.NRsta.curIt
            //     << "   " << flag << "   " << k << "   " << 2 << "   ";
            EoSctrl.SSMsta.conflag = OCP_TRUE;
            EoSctrl.SSMsta.realTol = sqrt(Se);
            return OCP_FALSE;
        }
    }
    // cout << "Yt = " << setprecision(12) << scientific << Yt << "    " << setw(3)
    //     << EoSctrl.SSMsta.curIt << "    " << setw(3) << EoSctrl.NRsta.curIt << "   "
    //     << flag << "   " << k << "   " << 1 << "   ";
    /*if (!flag) {
        OCP_WARNING("SSM not converged in Stability Analysis");
    }*/
    EoSctrl.SSMsta.realTol = sqrt(Se);
    return OCP_TRUE;
}

OCP_BOOL MixtureComp::StableNR(const USI& Id)
{

    for (USI i = 0; i < NC; i++) {
        resSTA[i] = log(fug[Id][i] / (fugSta[i] * Yt));
    }

    USI     maxIt = EoSctrl.NRsta.maxIt;
    OCP_DBL Stol  = EoSctrl.NRsta.tol;
    OCP_DBL Se    = Dnorm2(NC, &resSTA[0]);
    OCP_DBL alpha = 1;
    USI     iter  = 0;

    while (Se > Stol) {

        eos.CalLnFugX(P, T, Y, fugX[0]);
        AssembleJmatSTA();
        // LUSolve(1, NC, &JmatSTA[0], &resSTA[0], &pivot[0]);
        SYSSolve(1, &uplo, NC, &JmatSTA[0], &resSTA[0], &pivot[0], &JmatWork[0],
                 lJmatWork);

        const OCP_DBL lYt = Yt;
        Yt = 0;
        for (USI i = 0; i < NC; i++) {
            Y[i] = Y[i] * lYt + alpha * resSTA[i];
            Yt   += Y[i];
        }
        Dscalar(NC, 1 / Yt, &Y[0]);

        eos.CalFug(P, T, Y, fugSta);
        for (USI i = 0; i < NC; i++) {
            resSTA[i] = log(fug[Id][i] / (fugSta[i] * Yt));
        }
        Se = Dnorm2(NC, &resSTA[0]);
        iter++;
        if (iter > maxIt) {
            EoSctrl.NRsta.curIt += iter;
            EoSctrl.NRsta.conflag = OCP_FALSE;
            EoSctrl.NRsta.realTol = Se;
            return OCP_FALSE;
        }
    }
    EoSctrl.NRsta.curIt += iter;
    EoSctrl.NRsta.conflag = OCP_TRUE;
    EoSctrl.NRsta.realTol = Se;
    return OCP_TRUE;
}

void MixtureComp::AssembleJmatSTA()
{
    vector<OCP_DBL>& fugx = fugX[0];
    fill(JmatSTA.begin(), JmatSTA.end(), 0.0);
    OCP_DBL tmp;
    for (USI i = 0; i < NC; i++) {

        tmp = 0;
        for (USI k = 0; k < NC; k++) {
            tmp += Y[k] * fugx[i * NC + k];
        }

        for (USI j = 0; j < NC; j++) {
            // Symmetric
            // JmatSTA[i * NC + j] = (fugx[i * NC + j] - tmp + delta(i, j) / Y[i]) / Yt;
            JmatSTA[i * NC + j] = (fugx[i * NC + j] - tmp + 1) / Yt;
        }
    }
}

void MixtureComp::PhaseSplit()
{
    EoSctrl.SSMsp.conflag = OCP_FALSE;
    EoSctrl.SSMsp.curIt   = 0;
    EoSctrl.NRsp.conflag  = OCP_FALSE;
    EoSctrl.NRsp.curIt    = 0;
    EoSctrl.RR.curIt      = 0;

    // cout << "begin" << endl;
    SplitSSM(OCP_FALSE);
    SplitNR();
    while (!EoSctrl.NRsp.conflag) {
        SplitSSM(OCP_TRUE);
        SplitNR();
        if (!CheckSplit()) break;
        if (EoSctrl.SSMsp.conflag) break;
    }
    CheckSplit();

    itersSSMSP += EoSctrl.SSMsp.curIt;
    itersNRSP += EoSctrl.NRsp.curIt;
    itersRR += EoSctrl.RR.curIt;
    countsSSMSP++;
    countsNRSP++;
}

OCP_BOOL MixtureComp::CheckSplit()
{
    if (NP == 2) {

        OCP_DBL eX    = 0;
        OCP_DBL nuMax = max(nu[0], nu[1]);

        for (USI i = 0; i < NC; i++) {
            eX += (x[0][i] - x[1][i]) * (x[0][i] - x[1][i]);
        }

        if (OCP_TRUE) {
            // Calculate Gibbs Energy

            eos.CalFug(P, T, zi, fugSta);
            GibbsEnergyB = 0;
            GibbsEnergyE = 0;
            for (USI i = 0; i < NC; i++) {
                GibbsEnergyB += zi[i] * log(fugSta[i]);
                GibbsEnergyE += (n[0][i] * log(fug[0][i]) + n[1][i] * log(fug[1][i]));
            }

            cout << scientific << setprecision(6);
            // cout << GibbsEnergyB << "   " << GibbsEnergyE << endl;
            if (GibbsEnergyE > GibbsEnergyB) {
                cout << ftype << "   ";
                cout << GibbsEnergyB << "   ";
                cout << GibbsEnergyE << "   ";
                cout << nuMax << "   ";
                cout << eX << "   ";
                cout << EoSctrl.NRsp.conflag << "   ";
                cout << bulkId << endl;
            }
        }

        if (nuMax < 1 && EoSctrl.NRsp.conflag && isfinite(eX)) {
            // accept this result
        } else {
            if (!isfinite(eX) || 1 - nuMax < 1E-3) {
                NP    = 1;
                x[0]  = zi;
                nu[0] = 1;

                EoSctrl.SSMsta.conflag = OCP_FALSE;
                EoSctrl.NRsta.conflag  = OCP_FALSE;
                return OCP_FALSE;
            }
        }
    }
    return OCP_TRUE;
}

void MixtureComp::SplitSSM(const OCP_BOOL& flag)
{
    if (NP == 2) {
        SplitSSM2(flag);
    } else {
        SplitSSM3(flag);
    }
}

void MixtureComp::SplitSSM2(const OCP_BOOL& flag)
{
    // NP = 2 in this case
    // Ks is very IMPORTANT!
    // flag = OCP_TRUE : Restart SSM
    // flag = OCP_FALSE : New SSM
    EoSctrl.SSMsp.conflag = OCP_TRUE;
    OCP_DBL Se            = 1;
    OCP_DBL Stol          = EoSctrl.SSMsp.tol2;
    USI     maxIt         = EoSctrl.SSMsp.maxIt;

    if (!flag) {
        if (lNP == 2) {
            Ks[NP - 2] = lKs;
        } else {
            if (Yt < 1.1 || OCP_TRUE) {
                Ks[NP - 2] = Kw[0];
            } else {
                for (USI i = 0; i < NC; i++) {
                    // Ks[NP - 2][i] = phi[testPId][i] / phiSta[i];
                    Ks[NP - 2][i] = Y[i] / x[testPId][i];
                }
            }
        }
    } else {
        Stol *= 1E-1;
        maxIt *= 2;
    }

    if (Yt < 1.1) {
        Stol *= 1E-1;
        maxIt *= 2;
    }

    USI iter = 0;
    while (Se > Stol) {

        RachfordRice2();
        UpdateXRR();
        for (USI j = 0; j < NP; j++)
            eos.CalFugPhi(P, T, x[j], fug[j], phi[j]);

        Se = 0;
        for (USI i = 0; i < NC; i++) {
            Se += pow(fug[1][i] / fug[0][i] - 1, 2);
            Ks[0][i] = phi[1][i] / phi[0][i];
        }

        iter++;
        if (iter > maxIt) {
            // OCP_WARNING("SSM not converged in Phase Spliting!");
            EoSctrl.SSMsp.conflag = OCP_FALSE;
            break;
        }
    }

    EoSctrl.SSMsp.realTol = sqrt(Se);
    EoSctrl.SSMsp.curIt += iter;
}

void MixtureComp::SplitSSM3(const OCP_BOOL& flag) {}

void MixtureComp::RachfordRice2() ///< Used when NP = 2
{
    const vector<OCP_DBL>& Ktmp = Ks[0];
    OCP_DBL                Kmin = Ktmp[0];
    OCP_DBL                Kmax = Ktmp[0];

    for (USI i = 1; i < NC; i++) {
        if (Ktmp[i] < Kmin) Kmin = Ktmp[i];
        if (Ktmp[i] > Kmax) Kmax = Ktmp[i];
    }

    const OCP_DBL numin = 1 / (1 - Kmax);
    const OCP_DBL numax = 1 / (1 - Kmin);

    nu[0] = 0.5 * (numin + numax);

    // Solve RR with NR
    OCP_DBL tmp, rj, J, dnuj, tmpnu;

    USI           iter  = 0;
    const OCP_DBL tol   = EoSctrl.RR.tol;
    const OCP_DBL maxIt = EoSctrl.RR.maxIt;
    while (OCP_TRUE) {

        rj = 0;
        J  = 0;
        for (USI i = 0; i < NC; i++) {
            tmp = 1 + nu[0] * (Ktmp[i] - 1);
            rj += zi[i] * (Ktmp[i] - 1) / tmp;
            J -= zi[i] * (Ktmp[i] - 1) * (Ktmp[i] - 1) / (tmp * tmp);
        }

        if (fabs(rj) < tol || iter > maxIt) break;

        dnuj  = -rj / J;
        tmpnu = nu[0] + dnuj;
        if (tmpnu < numax && tmpnu > numin) {
            nu[0] = tmpnu;
        } else {
            if (dnuj > 0) {
                nu[0] = (nu[0] + numax) / 2;
            } else {
                nu[0] = (nu[0] + numin) / 2;
            }
        }
        iter++;
    }

    EoSctrl.RR.curIt += iter;
    nu[1] = 1 - nu[0];

    // cout << scientific << setprecision(6) << nu[0] << "   " << nu[1] << endl;
}

void MixtureComp::RachfordRice2P()
{
    // modified RachfordRice equations
    // less iterations but more divergence --- unstable!

    const vector<OCP_DBL>& Ktmp = Ks[0];
    OCP_DBL                Kmin = Ktmp[0];
    OCP_DBL                Kmax = Ktmp[0];

    for (USI i = 1; i < NC; i++) {
        if (Ktmp[i] < Kmin) Kmin = Ktmp[i];
        if (Ktmp[i] > Kmax) Kmax = Ktmp[i];
    }

    const OCP_DBL numin = 1 / (1 - Kmax);
    const OCP_DBL numax = 1 / (1 - Kmin);

    nu[0] = 0.5 * (numin + numax);

    // Solve RR with NR
    OCP_DBL tmp, rj, J, dnuj, tmpnu;
    OCP_DBL f, df;

    USI           iter  = 0;
    const OCP_DBL tol   = EoSctrl.RR.tol;
    const OCP_DBL maxIt = EoSctrl.RR.maxIt;
    while (OCP_TRUE) {

        rj = 0;
        J  = 0;
        for (USI i = 0; i < NC; i++) {
            tmp = 1 + nu[0] * (Ktmp[i] - 1);
            rj += zi[i] * (Ktmp[i] - 1) / tmp;
            J -= zi[i] * (Ktmp[i] - 1) * (Ktmp[i] - 1) / (tmp * tmp);
        }
        f  = (nu[0] - numin) * (numax - nu[0]);
        df = -2 * nu[0] + (numax + numin);
        J *= f;
        J += df * rj;
        rj *= f;

        if (fabs(rj) < tol || iter > maxIt) break;

        dnuj  = -rj / J;
        tmpnu = nu[0] + dnuj;
        if (tmpnu < numax && tmpnu > numin) {
            nu[0] = tmpnu;
        } else {
            if (dnuj > 0) {
                nu[0] = (nu[0] + numax) / 2;
            } else {
                nu[0] = (nu[0] + numin) / 2;
            }
        }
        iter++;
    }

    EoSctrl.RR.curIt += iter;

    cout << iter << "      " << scientific << setprecision(3) << fabs(rj) << "      "
         << nu[0] << "      " << numin << "      " << numax << endl;

    nu[1] = 1 - nu[0];
}

void MixtureComp::RachfordRice3() ///< Used when NP > 2
{
}

void MixtureComp::UpdateXRR()
{
    OCP_DBL tmp = 0;
    for (USI i = 0; i < NC; i++) {
        tmp = 1;
        for (USI j = 0; j < NP - 1; j++) {
            tmp += nu[j] * (Ks[j][i] - 1);
        }
        x[NP - 1][i] = zi[i] / tmp;
        for (USI j = 0; j < NP - 1; j++) {
            x[j][i] = Ks[j][i] * x[NP - 1][i];
        }
    }
}

void MixtureComp::SplitBFGS()
{
    // Use BFGS to calculate phase splitting
    // JmatSP will store the BFGS mat
    // resSP will store the resSP - lresSP if necessary
    // n will store the n - ln if necessary

    // get initial value, n and ln, resSP and lresSP, H0

    // begin BFGS
}

void MixtureComp::SplitNR()
{
    EoSctrl.NRsp.conflag = OCP_FALSE;
    // for (USI j = 0; j < NP; j++) {
    //     nu[j] = fabs(nu[j]);
    // }

    USI len = NC * (NP - 1);
    x2n();
    CalResSP();
    OCP_DBL eNR0;
    OCP_DBL eNR   = Dnorm2(len, &resSP[0]);
    OCP_DBL NRtol = EoSctrl.NRsp.tol;
    OCP_DBL alpha;

    OCP_DBL en;
    USI     iter = 0;
    eNR0         = eNR;
    while (eNR > NRtol) {
     
        // eNR0 = eNR;
        ln = n;        
        for (USI j = 0; j < NP; j++)  
            eos.CalLnFugN(P, T, x[j], nu[j], fugN[j]);

        AssembleJmatSP();

        // LUSolve(1, len, &JmatSP[0], &resSP[0], &pivot[0]);
        // PrintDX(NC, &resSP[0]);

        int info = SYSSolve(1, &uplo, len, &JmatSP[0], &resSP[0], &pivot[0],
                            &JmatWork[0], lJmatWork);
        if (info > 0 && false) {
            for (USI i = 0; i < len; i++) {
                for (USI j = 0; j < len; j++) {
                    cout << scientific << setprecision(9) << JmatSP[i * len + j]
                         << "   ";
                }
                cout << ";\n";
            }
        }
        // PrintDX(NC, &resSP[0]);

        alpha = CalStepNRsp();

        n[NP - 1] = zi;       
        nu[NP - 1] = 1;
        for (USI j = 0; j < NP - 1; j++) {
            nu[j] = 0;
            for (USI i = 0; i < NC; i++) {
                n[j][i]         += alpha * resSP[j * NC + i];                
                n[NP - 1][i]    -= n[j][i];
                nu[j]           += n[j][i];
            }
            nu[NP - 1] -= nu[j];

            for (USI i = 0; i < NC; i++) {
                // x[j][i] = n[j][i] / nu[j];
                x[j][i] = fabs(n[j][i] / nu[j]);
            }
        }
        for (USI i = 0; i < NC; i++) {
            x[NP - 1][i] = fabs(n[NP - 1][i] / nu[NP - 1]);
        }

        for (USI j = 0; j < NP; j++)
            eos.CalFug(P, T, x[j], fug[j]);

        CalResSP();
        eNR = Dnorm2(len, &resSP[0]);
        iter++;
        if (eNR > eNR0 || iter > EoSctrl.NRsp.maxIt) {
            break;
        }

        // Maybe it should be execute before "eNR > eNR0 || iter > EoSctrl.NRsp.maxIt"
        en = 0;
        for (USI j = 0; j < NP; j++) {
            for (USI i = 0; i < NC; i++) {
                en += (ln[j][i] - n[j][i]) * (ln[j][i] - n[j][i]);
            }
        }
        if (en < 1E-16 * (NP * NC) * (NP * NC)) {
            EoSctrl.NRsp.conflag = OCP_TRUE;
            break;
        }
    }
    EoSctrl.NRsp.realTol = eNR;
    if (eNR < NRtol) EoSctrl.NRsp.conflag = OCP_TRUE;
    EoSctrl.NRsp.curIt += iter;

    // cout << iter << "   " << scientific << setprecision(3) << eNR << endl;
}

void MixtureComp::CalResSP()
{
    // So it equals -res
    for (USI j = 0; j < NP - 1; j++) {
        for (USI i = 0; i < NC; i++) {
            resSP[j * NC + i] = log(fug[NP - 1][i] / fug[j][i]);
        }
    }
}


void MixtureComp::PrintFugN()
{
    for (USI j = 0; j < NP; j++) {
        for (USI i = 0; i < NC; i++) {
            for (USI k = 0; k < NC; k++) {
                cout << fugN[j][i * NC + k] << "   ";
            }
            cout << endl;
        }
        cout << "---------------" << endl;
    }
}

void MixtureComp::AssembleJmatSP()
{
    // Dim: (NP-1)*NC
    // Attention that fugNj is sysmetric
    fill(JmatSP.begin(), JmatSP.end(), 0);

    OCP_DBL*       Jtmp = &JmatSP[0];
    const OCP_DBL* fugNp;
    const OCP_DBL* fugNj;

    for (USI j = 0; j < NP - 1; j++) {
        fugNp = &fugN[NP - 1][0];
        fugNj = &fugN[j][0];

        for (USI i = 0; i < NC; i++) {
            // ith components
            for (USI k = 0; k < NC; k++) {
                // kth fugacity
                Jtmp[k] = fugNj[k] + fugNp[k];
            }
            Jtmp += NC * (NP - 1);
            fugNp += NC;
            fugNj += NC;
        }
        Jtmp += NC;
    }

    // cout << endl << "Jmat" << endl;
    // for (USI i = 0; i < NC; i++) {
    //     for (USI j = 0; j < NC; j++) {
    //         cout << scientific << setprecision(6) << JmatSP[i * NC + j] << "   ";
    //     }
    //     cout << endl;
    // }
}

OCP_DBL MixtureComp::CalStepNRsp()
{
    OCP_DBL alpha = 1;
    OCP_DBL tmp;

    for (USI j = 0; j < NP - 1; j++) {

        const OCP_DBL* nj = &n[j][0];
        const OCP_DBL* r  = &resSP[j * NC];

        for (USI i = 0; i < NC; i++) {
            tmp = nj[i] + alpha * r[i];
            if (tmp < 0) {
                alpha *= 0.9 * fabs(nj[i] / r[i]);
            }
        }
    }
    return alpha;
}

void MixtureComp::AllocateOthers()
{
    sqrtMWi.resize(NC);
    for (USI i = 0; i < NC; i++) sqrtMWi[i] = sqrt(MWC[i]);
    muC.resize(NPmax);
    muAux.resize(NPmax);
    for (USI i = 0; i < NPmax; i++) {
        muAux[i].resize(5);
    }
    muAux1I.resize(NC);
    for (USI i = 0; i < NC; i++) {
        muAux1I[i] = 5.4402 * pow(Tc[i], 1.0 / 6) / pow(Pc[i], 2.0 / 3);
    }
    fugP.resize(NPmax);
    for (USI j = 0; j < NPmax; j++) {
        fugP[j].resize(NC);
    }
    Zp.resize(NPmax);
    JmatDer.resize(NPmax * NPmax * (NC + 1) * (NC + 1));
    JmatTmp = JmatDer;
    rhsDer.resize(NPmax * (NC + 1) * (NC + 1));

    // new
    JmatDer.resize((numPhase + NPmax * NC) * (numPhase + NPmax * NC));
    JmatTmp = JmatDer;
    rhsDer.resize((numPhase + NPmax * NC) * (numCom + 1 + 1));
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
            A += x[0][i] * Vc[i] * Tc[i];
            B += x[0][i] * Vc[i];
        }
        OCP_DBL Tc = A / B;
        if (T > Tc) {
            phaseLabel[0] = GAS;
            phaseExist[GAS] = OCP_TRUE;
        } else {
            phaseLabel[0] = OIL;
            phaseExist[OIL] = OCP_TRUE;
        }
    } else {
        // Compare MW
        if (MW[0] > MW[1]) {
            phaseLabel[0] = OIL;
            phaseLabel[1] = GAS;
        } else {
            phaseLabel[0] = GAS;
            phaseLabel[1] = OIL;
        }
        phaseExist[OIL] = OCP_TRUE;
        phaseExist[GAS] = OCP_TRUE;
    }
}

void MixtureComp::CopyPhase()
{
    // copy vj, x, mu, xi, rho
    for (USI j = 0; j < NP; j++) {
        const USI j1 = phaseLabel[j];
        vj[j1]       = vC[j];
        copy(x[j].begin(), x[j].end(), &xij[j1 * numCom]);
        mu[j1]  = muC[j];
        xi[j1]  = xiC[j];
        rho[j1] = rhoC[j];
        nj[j1]  = nu[j];
    }
}

void MixtureComp::CalViscosity() { CalViscoLBC(); }

void MixtureComp::CalViscoLBC()
{
    OCP_DBL tmp;
    OCP_DBL Tri;
    OCP_DBL xijT;
    OCP_DBL xijP;
    OCP_DBL xijV;

    for (USI j = 0; j < NP; j++) {
        const vector<OCP_DBL>& xj  = x[j];
        vector<OCP_DBL>&       muA = muAux[j];
        fill(muA.begin(), muA.end(), 0.0);
        xijT = 0;
        xijP = 0;
        xijV = 0;

        for (USI i = 0; i < NC; i++) {
            tmp = 5.4402 * pow(Tc[i], 1.0 / 6) / sqrt(MW[j]) / pow(Pc[i], 2.0 / 3);
            // tmp = muAux1I[i] / sqrt(MW[j]);
            Tri = T / Tc[i];
            if (Tri <= 1.5) {
                tmp = 34 * 1E-5 * pow(Tri, 0.94) / tmp;
            } else {
                tmp = 17.78 * 1E-5 * pow((4.58 * Tri - 1.67), 0.625) / tmp;
            }
            muA[0] += xj[i] * sqrt(MWC[i]) * tmp;
            muA[1] += xj[i] * sqrt(MWC[i]);
            xijT += xj[i] * Tc[i];
            xijP += xj[i] * Pc[i];
            xijV += xj[i] * Vcvis[i];
        }
        muA[2] = 5.4402 * pow(xijT, 1.0 / 6) / sqrt(MW[j]) / pow(xijP, 2.0 / 3);
        muA[3] = xiC[j] * xijV;

        if (muA[3] <= 0.18 && OCP_FALSE) {
            muC[j] = muA[0] / muA[1] + 2.05 * 1E-4 * muA[3] / muA[2];
        } else {
            muA[4] = muA[3] * (muA[3] * (muA[3] * (LBCcoef[4] * muA[3] + LBCcoef[3]) +
                                         LBCcoef[2]) +
                               LBCcoef[1]) +
                     LBCcoef[0];
            muC[j] = muA[0] / muA[1] + 1E-4 * (pow(muA[4], 4) - 1) / muA[2];
        }
    }
}

void MixtureComp::CalViscoHZYT() {}


void MixtureComp::CalXiPNX_partial()
{
    // Calculate xix and xiP
    for (USI j = 0; j < NP; j++) {
        const USI j0 = phaseLabel[j];
        xi[j0] = 1 / vmj[j];
        for (USI i = 0; i < NC; i++) {
            xix[j0 * numCom + i] = -xi[j0] * xi[j0] * vmx[j][i];
        }
        xiP[j0] = -xi[j0] * xi[j0] * vmP[j];
    }
}

void MixtureComp::CalRhoPX_partial()
{
    for (USI j = 0; j < NP; j++) {
        const USI j1 = phaseLabel[j];
        rhoP[j1] = xiP[j1] * MW[j];
        for (USI i = 0; i < NC; i++) {
            rhox[j1 * numCom + i] = xix[j1 * numCom + i] * MW[j] + xi[j1] * MWC[i];
        }
    }
}

void MixtureComp::CalMuPX_partial()
{
    CalMuPXLBC_partial();
}

void MixtureComp::CalMuPXLBC_partial()
{
    // Here!
    // dmu/dP = dmu/dP
    // dmu/dNk = dmu/dNk = 0
    // // dxij / dP = dxij / dNk = 0
    // See MixtureComp::CalMuPXLBC_full01()

    OCP_DBL val1IJ, val2IJ;
    OCP_DBL der1IJ, der2IJ, der3J, der4J, der6J, der7J, der8J;
    OCP_DBL Tri, tmp;
    OCP_DBL xTj, xPj, xVj;
    OCP_DBL derxTj, derxPj, derMWj;

    // dxij / dP = 0, d MJ / dP = 0
    // Calculate dmuj / dP
    // der2IJ = der3J = der4J = der6J = 0;

    for (USI j = 0; j < NP; j++) {
        const USI              j1     = phaseLabel[j];
        const vector<OCP_DBL>& xj     = x[j];
        const vector<OCP_DBL>& muAuxj = muAux[j];

        xTj = xPj = xVj = 0;
        for (USI i = 0; i < NC; i++) {
            xTj += xj[i] * Tc[i];
            xPj += xj[i] * Pc[i];
            xVj += xj[i] * Vcvis[i];
        }
        der7J = xVj * xiP[j1];
        if (muAuxj[3] <= 0.18 && OCP_FALSE) {
            muP[j1] = (2.05 * 1E-4) * der7J / muAuxj[2];
        } else {
            der8J   = der7J * (LBCcoef[1] +
                             muAuxj[3] * (2 * LBCcoef[2] +
                                          muAuxj[3] * (3 * LBCcoef[3] +
                                                       muAuxj[3] * 4 * LBCcoef[4])));
            muP[j1] = (4 * 1E-4) * pow(muAuxj[4], 3) * der8J / muAuxj[2];
        }

        // Calculate dmuj / xkj
        const USI bId = numCom * j1;
        for (USI k = 0; k < NC; k++) {
            derxTj = Tc[k];
            derxPj = Pc[k];
            derMWj = MWC[k];
            der3J  = 0;
            for (USI i = 0; i < NC; i++) {
                val1IJ = muAux1I[i] / sqrt(MW[j]);
                der1IJ = -(1 / 2) * muAux1I[i] * pow(MW[j], -1.5) * derMWj;
                Tri    = T / Tc[i];
                if (Tri <= 1.5) {
                    tmp = 34 * 1E-5 * pow(Tri, 0.94);
                } else {
                    tmp = 17.78 * 1E-5 * pow(4.58 * Tri - 1.67, 0.625);
                }
                val2IJ = tmp / val1IJ;
                der2IJ = -tmp * der1IJ / (val1IJ * val1IJ);
                der3J +=
                    xj[i] * sqrtMWi[i] * der2IJ + delta(i, k) * sqrtMWi[k] * val2IJ;
            }
            der4J = sqrtMWi[k];
            der6J =
                5.4402 *
                (1.0 / 6 * pow(xTj, -5.0 / 6) * derxTj -
                 pow(xTj, 1.0 / 6) * (0.5 / MW[j] * derMWj + 2.0 / 3 / xPj * derxPj)) /
                (sqrt(MW[j]) * pow(xPj, 2.0 / 3));
            der7J = xix[j1 * numCom + k] * xVj + xi[j1] * Vcvis[k];
            if (muAuxj[3] <= 0.18 && OCP_FALSE) {
                mux[bId + k] =
                    (der3J * muAuxj[1] - muAuxj[0] * der4J) / (muAuxj[1] * muAuxj[1]) +
                    2.05 * 1E-4 * (der7J * muAuxj[2] - muAuxj[3] * der6J) /
                        (muAuxj[2] * muAuxj[2]);
            } else {
                der8J =
                    der7J * (LBCcoef[1] +
                             muAuxj[3] * (2 * LBCcoef[2] +
                                          muAuxj[3] * (3 * LBCcoef[3] +
                                                       muAuxj[3] * 4 * LBCcoef[4])));
                mux[bId + k] =
                    (der3J * muAuxj[1] - muAuxj[0] * der4J) / (muAuxj[1] * muAuxj[1]) +
                    1E-4 *
                        (4 * pow(muAuxj[4], 3) * der8J * muAuxj[2] -
                         (pow(muAuxj[4], 4) - 1) * der6J) /
                        (muAuxj[2] * muAuxj[2]);
            }
        }
    }
}


void MixtureComp::CalVfiVfp_full01()
{

    for (USI j = 0; j < NP; j++) {
        vmj[j] = eos.CalVmDer(P, T, x[j], vmP[j], vmx[j]);
    }

    if (NP == 1) {
        const USI j0 = phaseLabel[0];
        vjp[j0]      = nu[0] * vmP[0];
        vfP          = vjp[j0];
        for (USI i = 0; i < NC; i++) {
            vji[j0][i] = vmj[0] + vmx[0][i];
            for (USI k = 0; k < NC; k++) {
                vji[j0][i] -= vmx[0][k] * x[0][k];
            }
            vfi[i] = vji[j0][i];
        }
    } else {
        // NP > 1
        for (USI j = 0; j < NP; j++) {
            eos.CalLnFugN(P, T, x[j], nu[j], fugN[j]);
            eos.CalLnFugP(P, T, x[j], fugP[j]);
        }
        AssembleMatVfiVfp_full01();
        AssembleRhsVfiVfp_full01();
        LUSolve(NC + 1, NC * NP, &JmatDer[0], &rhsDer[0], &pivot[0]);

        // now d nm0 / dP(dNk) has been available
        const OCP_DBL* dnkjdP = &rhsDer[0];
        // Calculate Vfp
        vfP = 0;
        vector<OCP_DBL> tmp(NP, 0);
        for (USI j = 0; j < NP; j++) {
            const USI j0 = phaseLabel[j];
            for (USI i = 0; i < NC; i++) {
                tmp[j] -= vmx[j][i] * x[j][i];
            }          
            vjp[j0] = nu[j] * vmP[j];
            for (USI k = 0; k < NC; k++) {
                vjp[j0] += (vmj[j] + tmp[j] + vmx[j][k]) * dnkjdP[k];
            }
            dnkjdP += NC;
            vfP += vjp[j0];
        }

        // Calculate Vfi
        const OCP_DBL* dnkjdN = dnkjdP;
        for (USI i = 0; i < NC; i++) {
            vfi[i] = 0;
            for (USI j = 0; j < NP; j++) {
                const USI j0 = phaseLabel[j];
                vji[j0][i] = 0;
                for (USI k = 0; k < NC; k++) {
                    vji[j0][i] += (vmj[j] + tmp[j] + vmx[j][k]) * dnkjdN[k];
                }
                vfi[i] += vji[j0][i];
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
        OCP_DBL* fugNj = &fugN[j][0];
        for (USI i = 0; i < NC; i++) {
            // for ith components
            Dcopy(NC, bId, fugNj);
            bId += (NP - 1 - j) * NC;
            bId[i] = 1.0;
            bId += (1 + j) * NC;
            fugNj += NC;
        }
        bId += NC;
    }
    // NP - 1 phase
    bId            = &JmatDer[(NP - 1) * (NP * NC * NC)];
    OCP_DBL* fugNj = &fugN[NP - 1][0];
    for (USI i = 0; i < NC; i++) {
        for (USI j = 0; j < NP - 1; j++) {
            Dcopy(NC, bId, fugNj);
            bId += NC;
        }
        Dscalar((NP - 1) * NC, -1.0, bId - (NP - 1) * NC);
        fugNj += NC;
        bId[i] = 1.0;
        bId += NC;
    }
}

void MixtureComp::AssembleRhsVfiVfp_full01()
{
    fill(rhsDer.begin(), rhsDer.end(), 0.0);
    OCP_DBL* rhstmp = &rhsDer[0];
    // d P
    for (USI j = 0; j < NP; j++) {
        for (USI i = 0; i < NC; i++) {
            rhstmp[i] = fugP[NP - 1][i] - fugP[j][i];
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

    vector<OCP_INT> pLabel(numPhase, -1);
    for (USI j = 0; j < NP; j++) {
        pLabel[j] = phaseLabel[j];
    }
    pLabel.back() = WATER;

    // dS / dP, dS / dN
    for (USI j = 0; j < numPhase; j++) {
        const OCP_INT j0 = pLabel[j];
        if (j0 >= 0) {
            OCP_DBL* bId = &dXsdXp[j0 * ncol];
            // dS / dP
            bId[0] = (vjp[j0] * vf - vfP * vj[j0]) / vf2;
            bId++;
            // dS / dN
            for (USI m = 0; m < numCom; m++) {
                bId[m] = (vji[j0][m] * vf - vfi[m] * vj[j0]) / vf2;
            }
        }
    }

    if (NP == 1) {
        // dxij / dNm
        OCP_DBL* bId = &dXsdXp[(numPhase + phaseLabel[0] * numCom) * ncol + 1];
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
        for (USI j = 0; j < NP; j++) {
            const USI j0   = phaseLabel[j];
            OCP_DBL*  bId0 = bId + j0 * numCom * ncol;                      
            OCP_DBL        njDerSum = 0;
            for (USI i = 0; i < NC; i++) njDerSum += dnkjdP[i];
            for (USI i = 0; i < NC; i++) {
                bId0[0] = (dnkjdP[i]  - njDerSum * x[j][i]) / nu[j];
                bId0    += ncol;
            }
            dnkjdP += NC;
        }
        // dxij / dNm
        const OCP_DBL* dnkjdN = dnkjdP;
        for (USI m = 0; m < NC; m++) {
            for (USI j = 0; j < NP; j++) {
                const USI j0       = phaseLabel[j];
                OCP_DBL   njDerSum = 0;
                for (USI i = 0; i < NC; i++) njDerSum += dnkjdN[i];
                OCP_DBL* bId0 = bId + j0 * numCom * ncol + m + 1;
                for (USI i = 0; i < NC; i++) {
                    bId0[0] = (dnkjdN[i] - njDerSum * x[j][i]) / nu[j];
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

    for (USI j = 0; j < NP; j++) {
        vmj[j] = eos.CalVmDer(P, T, x[j], vmP[j], vmx[j]);
    }

    if (NP == 1) {
        const USI j0 = phaseLabel[0];
        vjp[j0] = nu[0] * vmP[0];
        vfP     = vjp[j0];
        for (USI i = 0; i < NC; i++) {
            vji[j0][i] = vmj[0] + vmx[0][i];
            for (USI k = 0; k < NC; k++) {
                vji[j0][i] -= vmx[0][k] * x[0][k];
            }
            vfi[i] = vji[j0][i];
        }
    } 
    else if (NP == 2) {
        // JUST FOR NP = 2
        for (USI j = 0; j < NP; j++) {
            eos.CalLnFugN(P, T, x[j], nu[j], fugN[j]);
            eos.CalLnFugP(P, T, x[j], fugP[j]);
        }

        AssembleMatVfiVfp_full02();
        AssembleRhsVfiVfp_full02();
        LUSolve(NC + 1, NC, &JmatDer[0], &rhsDer[0], &pivot[0]);

        // now d nm0 / dP(dNk) has been available
        OCP_DBL tmp0 = 0;
        OCP_DBL tmp1 = 0;
        for (USI i = 0; i < NC; i++) {
            tmp0 -= vmx[0][i] * x[0][i];
            tmp1 -= vmx[1][i] * x[1][i];
        }
        const USI j0 = phaseLabel[0];
        const USI j1 = phaseLabel[1];
        // vfP
        const OCP_DBL* dnkjdP = &rhsDer[0];
        vjp[j0] = nu[0] * vmP[0];
        vjp[j1] = nu[1] * vmP[1];
        for (USI k = 0; k < NC; k++) {
            vjp[j0] += (vmj[0] + tmp0 + vmx[0][k]) * dnkjdP[k];
            vjp[j1] -= (vmj[1] + tmp1 + vmx[1][k]) * dnkjdP[k];
        }
        vfP = vjp[j0] + vjp[j1];

        // Calculate Vfi
        const OCP_DBL* dnkjdN = dnkjdP + NC;
        for (USI i = 0; i < NC; i++) {
            vji[j0][i] = 0;
            vji[j1][i] = 0;
            for (USI k = 0; k < NC; k++) {
                vji[j0][i] += (vmj[0] + tmp0 + vmx[0][k]) * dnkjdN[k];
                vji[j1][i] += (vmj[1] + tmp1 + vmx[1][k]) * (delta(i, k) - dnkjdN[k]);
            }
            vfi[i] = vji[j0][i] + vji[j1][i];
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
            tmpMat[m] = fugN[0][i * NC + m] + fugN[1][i * NC + m];
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
        rhstmp[i] = fugP[1][i] - fugP[0][i];
    }
    rhstmp += NC;

    // dNk
    for (USI k = 0; k < NC; k++) {
        for (USI i = 0; i < NC; i++) {
            rhstmp[i] = fugN[1][k * NC + i]; // d lnfij / d nkj = d lnfkj / d nij
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
    const USI ncol = numCom + 1;
    const OCP_DBL vf2 = vf * vf;

    vector<OCP_INT> pLabel(numPhase, -1);
    for (USI j = 0; j < NP; j++) {
        pLabel[j] = phaseLabel[j];
    }
    pLabel.back() = WATER;

    for (USI j = 0; j < numPhase; j++) {
        const OCP_INT j0 = pLabel[j];
        if (j0 >= 0) {
            OCP_DBL* bId = &dXsdXp[j0 * ncol];
            // dS / dP
            bId[0] = (vjp[j0] * vf - vfP * vj[j0]) / vf2;
            bId++;
            // dS / dN
            for (USI m = 0; m < numCom; m++) {
                bId[m] = (vji[j0][m] * vf - vfi[m] * vj[j0]) / vf2;
            }
        }
    }

    if (NP == 1) {
        // dxij / dNm
        OCP_DBL* bId = &dXsdXp[(numPhase + phaseLabel[0] * numCom) * ncol + 1];
        for (USI i = 0; i < NC; i++) {
            for (USI m = 0; m < NC; m++) {
                bId[m] = (delta(i, m) * Nh - Ni[i]) / (Nh * Nh);
            }
            bId += ncol;
        }
    } 
    else if (NP == 2)  {
        // NP = 2
        OCP_DBL* bId = &dXsdXp[numPhase * ncol];
        // dxij / dP, dxij / dNm
        const USI      j0       = phaseLabel[0];
        const USI      j1       = phaseLabel[1];
        const OCP_DBL* dnkjdNP  = &rhsDer[0];
        OCP_DBL*       bId1     = bId + j0 * numCom * ncol;
        OCP_DBL*       bId2     = bId + j1 * numCom * ncol;
        OCP_DBL        njDerSum = 0;
        // dxij / dP
        for (USI i = 0; i < NC; i++) njDerSum += dnkjdNP[i];
        for (USI i = 0; i < NC; i++) {
            bId1[0] = (dnkjdNP[i] - njDerSum * x[0][i]) / nu[0];
            bId2[0] = (-dnkjdNP[i] + njDerSum * x[1][i]) / nu[1];
            bId1 += ncol;
            bId2 += ncol;
        }
        dnkjdNP += NC;
        // dxij / dNm
        for (USI m = 0; m < NC; m++) {
            njDerSum = 0;
            for (USI i = 0; i < NC; i++) njDerSum += dnkjdNP[i];
            bId1 = bId + j0 * numCom * ncol + m + 1;
            bId2 = bId + j1 * numCom * ncol + m + 1;
            for (USI i = 0; i < NC; i++) {
                bId1[0] = (dnkjdNP[i] - njDerSum * x[0][i]) / nu[0];
                bId2[0] = ((delta(i, m) - dnkjdNP[i]) - (1 - njDerSum) * x[1][i]) / nu[1];
                bId1 += ncol;
                bId2 += ncol;
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

/////////////////////////////////////////////////////////////////////
// Accelerate PVT
/////////////////////////////////////////////////////////////////////

void MixtureComp::AllocateSkip()
{
    phiN.resize(NC * NC);
    skipMatSTA.resize(NC * NC);
    eigenSkip.resize(NC);
    eigenWork.resize(2 * NC + 1);
}

void MixtureComp::CalPhiNSkip()
{
    eos.CalLnPhiN(P, T, x[0], nu[0], phiN);
}

void MixtureComp::AssembleSkipMatSTA()
{
    // Sysmetric Matrix
    // stored by colum
    vector<OCP_DBL>& xj = x[0];

    for (USI i = 0; i < NC; i++) {
        for (USI j = 0; j <= i; j++) {
            skipMatSTA[i * NC + j] =
                delta(i, j) + nu[0] * sqrt(xj[i] * xj[j]) * phiN[i * NC + j];
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
                stm01Params = SurTenMethod01Params(compParam.parachor.data[tarId], x[0].data(), x[1].data(),
                                               &xiC[0], &xiC[1], &NP);
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
    cout << "SSMSTA:     " << setw(12) << itersSSMSTA << setw(15)
        << itersSSMSTA * 1.0 / countsSSMSTA << endl;
    cout << "NRSTA:      " << setw(12) << itersNRSTA << setw(15)
        << itersNRSTA * 1.0 / countsNRSTA << endl;
    cout << "SSMSP:      " << setw(12) << itersSSMSP << setw(15)
        << itersSSMSP * 1.0 / countsSSMSP << endl;
    cout << "NRSP:       " << setw(12) << itersNRSP << setw(15)
        << itersNRSP * 1.0 / countsNRSP << endl;
    cout << "NRRR:       " << setw(12) << itersRR << setw(15)
        << itersRR * 1.0 / itersSSMSP << endl;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jan/05/2022      Create file                          */
/*----------------------------------------------------------------------------*/