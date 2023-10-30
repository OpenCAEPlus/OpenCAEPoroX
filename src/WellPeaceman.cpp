/*! \file    WellPeaceman.cpp
 *  \brief   Peaceman PeacemanWell class definition
 *  \author  Shizhe Li
 *  \date    Aug/17/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "WellPeaceman.hpp"


void PeacemanWell::InputPerfo(const WellParam& well, const Domain& domain, const USI& wId)
{
    OCP_FUNCNAME;

    numPerf = domain.GetPerfNum(wId);
    perf.resize(numPerf);
    USI pp = 0;
    for (USI p = 0; p < well.GetPerfNum(); p++) {

        const OCP_INT loc = domain.GetPerfLocation(wId, p);       
        if (loc < 0) {
            continue;
        }
        perf[pp].location = loc;
        perf[pp].WI = well.WI[p];
        perf[pp].radius = well.diameter[p] / 2.0;
        perf[pp].kh = well.kh[p];
        perf[pp].skinFactor = well.skinFactor[p];
        if (well.direction[p] == "X" || well.direction[p] == "x") {
            perf[pp].direction = PerfDirection::x;
        }
        else if (well.direction[p] == "Y" || well.direction[p] == "y") {
            perf[pp].direction = PerfDirection::y;
        }
        else if (well.direction[p] == "Z" || well.direction[p] == "z") {
            perf[pp].direction = PerfDirection::z;
        }
        else if (well.direction[p] == "usg") {
            perf[pp].direction = PerfDirection::usg;
        }
        else {
            OCP_ABORT("Wrong direction of perforations!");
        }
        pp++;
    }
    if (pp != numPerf) {
        OCP_ABORT(to_string(numPerf) + "  " + to_string(pp) +  " -- Wrong Perf Setup! Rank : " + to_string(domain.myrank));
    }  
}


void PeacemanWell::Setup(const Bulk& bk, const vector<SolventINJ>& sols)
{
    OCP_FUNCNAME;

    const BulkVarSet& bvs = bk.vs;

    nc = bvs.nc;
    np = bvs.np;
    mixture = bk.PVTm.GetMixture();
    rsTemp = bk.rsTemp;

    qi_lbmol.resize(nc);
    factor.resize(nc);

    SetupUnit();
    SetupOpts(sols);
    // Perf
    for (USI p = 0; p < numPerf; p++) {
        perf[p].state = WellState::open;
        perf[p].depth = bvs.depth[perf[p].location];
        perf[p].multiplier = 1;
        perf[p].qi_lbmol.resize(nc);
        perf[p].transj.resize(np);
        perf[p].qj_ft3.resize(np);
    }
    // dG
    dG.resize(numPerf, 0);

    if (depth < 0) depth = perf[0].depth;

    CalWI(bk);
    // test
    // ShowPerfStatus(bvs);
}



void PeacemanWell::InitWellP(const Bulk& bk)
{
    bhp = bk.vs.P[perf[0].location]; 
    CalPerfP();
}


void PeacemanWell::CheckOptMode(const Bulk& bk)
{
    if (opt.state != WellState::open)  return;

    CalTrans(bk);
    CaldG(bk);

    OCP_FUNCNAME;
    if (opt.initMode == WellOptMode::bhp) {
        if (opt.type == WellType::injector) {
            const OCP_DBL q = CalInjRateMaxBHP(bk);
            if (q > opt.maxRate) {
                opt.mode = WellOptMode::irate;
            }
            else {
                opt.mode = WellOptMode::bhp;
                bhp = opt.tarBHP;
            }
        }
        else {
            opt.mode = WellOptMode::bhp;
            bhp = opt.tarBHP;
        }
    }
    else {
        OCP_DBL q;
        if (opt.type == WellType::injector) {
            q = CalInjRateMaxBHP(bk);
        }
        else {
            q = CalProdRateMinBHP(bk);
        }
        if (q > opt.maxRate) {
            opt.mode = opt.initMode;
        }
        else {
            opt.mode = WellOptMode::bhp;
            bhp = opt.tarBHP;
        }
    }

    CalPerfP();
}


void PeacemanWell::CalFluxInit(const Bulk& bk) 
{
    if (opt.state != WellState::open)  return;

    CalTrans(bk);
    // recalculating seems worse
    // CaldG(bk);
    CalFlux(bk, OCP_TRUE);
}


void PeacemanWell::CalFlux(const Bulk& bk) 
{
    if (opt.state != WellState::open)  return;

    CalTrans(bk);
    CalFlux(bk, OCP_FALSE);
}


ReservoirState PeacemanWell::CheckP(const Bulk& bk)
{
    if (opt.state != WellState::open)  return ReservoirState::well_success;

    OCP_FUNCNAME;
    // 0 : all correct
    // 1 : negative P
    // 2 : outlimited P
    // 3 : crossflow happens

    if (bhp < 0) {
        cout << "### WARNING: Well " << name << " BHP = " << bhp << endl;
        return ReservoirState::well_negative_P;
    }
    for (USI p = 0; p < numPerf; p++) {
        if (perf[p].state == WellState::open && perf[p].P < 0) {
#ifdef DEBUG
            cout << "### WARNING: Well " << name << " Perf[" << p
                << "].P = " << perf[p].P << endl;
#endif // DEBUG
            return ReservoirState::well_negative_P;
        }
    }


    if (opt.mode != WellOptMode::bhp &&
        ((opt.type == WellType::injector && bhp > opt.maxBHP) ||
            (opt.type == WellType::productor && bhp < opt.minBHP))) {
#if _DEBUG
        cout << "### WARNING: Well " << name << " switch to BHPMode" << endl;
#endif
        opt.mode = WellOptMode::bhp;
        bhp      = opt.tarBHP;
        return ReservoirState::well_switch_BHPm;
    }

    return CheckCrossFlow(bk);
}



void PeacemanWell::CalIPRate(const Bulk& bk, const OCP_DBL& dt)
{
    WGIR = 0;
    WWIR = 0;
    WOPR = 0;
    WGPR = 0;
    WWPR = 0;

    if (opt.state != WellState::open)  return;

    if (opt.state == WellState::open) {
        if (opt.type == WellType::productor)  CalProdQj(bk, dt);
        else                                  CalInjQj(bk, dt);
    }
}


OCP_DBL PeacemanWell::CalMaxChangeTime() const
{ 
    if (opt.state != WellState::open)  return 0;
    else                               return (bhp - lbhp); 
}


OCP_DBL PeacemanWell::CalMaxChangeNR()
{
    if (opt.state != WellState::open)  return 0;
    else {
        const OCP_DBL dP = bhp - NRbhp;
        NRbhp = bhp;
        return dP;
    }
}


void PeacemanWell::ResetToLastTimeStep(const Bulk& bk)
{
    if (opt.state != WellState::open)  return;

    bhp   = lbhp;
    NRbhp = lbhp;
    if (opt.mode == WellOptMode::bhp) bhp = opt.tarBHP;
    CalPerfP();
    CalFluxInit(bk);
}


void PeacemanWell::UpdateLastTimeStep()
{ 
    lbhp  = bhp; 
    NRbhp = bhp;
}


void PeacemanWell::CalWI(const Bulk& bk)
{
    OCP_FUNCNAME;

    const BulkVarSet& bvs = bk.vs;

    // this fomular needs to be carefully checked !
    // especially the dz

    for (USI p = 0; p < numPerf; p++) {
        if (perf[p].WI > 0) {
            perf[p].WI *= CONV1;
            break;
        }
        else {
            if (perf[p].direction != PerfDirection::usg) {
                const OCP_USI Idb = perf[p].location;
                const OCP_DBL dx = bvs.dx[Idb];
                const OCP_DBL dy = bvs.dy[Idb];
                const OCP_DBL dz = bvs.dz[Idb] * bvs.ntg[Idb];
                OCP_DBL       ro = 0;
                switch (perf[p].direction) 
                {
                case PerfDirection::x:
                {
                    const OCP_DBL kykz = bvs.rockKy[Idb] * bvs.rockKz[Idb];
                    const OCP_DBL ky_kz = bvs.rockKy[Idb] / bvs.rockKz[Idb];
                    assert(kykz > 0);
                    ro = 0.28 * pow((dy * dy * pow(1 / ky_kz, 0.5) +
                        dz * dz * pow(ky_kz, 0.5)),
                        0.5);
                    ro /= (pow(ky_kz, 0.25) + pow(1 / ky_kz, 0.25));

                    if (perf[p].kh < 0) {
                        perf[p].kh = (dx * pow(kykz, 0.5));
                    }
                    break;
                }

                case PerfDirection::y:
                {
                    const OCP_DBL kzkx = bvs.rockKz[Idb] * bvs.rockKx[Idb];
                    const OCP_DBL kz_kx = bvs.rockKz[Idb] / bvs.rockKx[Idb];
                    assert(kzkx > 0);
                    ro = 0.28 * pow((dz * dz * pow(1 / kz_kx, 0.5) +
                        dx * dx * pow(kz_kx, 0.5)),
                        0.5);
                    ro /= (pow(kz_kx, 0.25) + pow(1 / kz_kx, 0.25));

                    if (perf[p].kh < 0) {
                        perf[p].kh = (dy * pow(kzkx, 0.5));
                    }
                    break;
                }

                case PerfDirection::z:
                {
                    const OCP_DBL kxky = bvs.rockKx[Idb] * bvs.rockKy[Idb];
                    const OCP_DBL kx_ky = bvs.rockKx[Idb] / bvs.rockKy[Idb];
                    assert(kxky > 0);
                    ro = 0.28 * pow((dx * dx * pow(1 / kx_ky, 0.5) +
                        dy * dy * pow(kx_ky, 0.5)),
                        0.5);
                    ro /= (pow(kx_ky, 0.25) + pow(1 / kx_ky, 0.25));

                    if (perf[p].kh < 0) {
                        perf[p].kh = (dz * pow(kxky, 0.5));
                    }
                    break;
                }
                
                default:
                    OCP_ABORT("Wrong direction of perforations!");
                }

                perf[p].WI = (2 * PI) * perf[p].kh /
                    (log(ro / perf[p].radius) + perf[p].skinFactor);
            }
            else {
                perf[p].WI = PI * perf[p].radius * perf[p].radius;
            }

            perf[p].WI *= CONV_DARCY;
        }
    }
}


void PeacemanWell::CalTrans(const Bulk& bk)
{
    OCP_FUNCNAME;

    const BulkVarSet& bvs = bk.vs;

    if (opt.type == WellType::injector) {
        for (USI p = 0; p < numPerf; p++) {
            perf[p].transINJ = 0;
            const OCP_USI k    = perf[p].location;
            const OCP_DBL temp = perf[p].WI * perf[p].multiplier;

            // single phase
            for (USI j = 0; j < np; j++) {
                perf[p].transj[j] = 0;
                const OCP_USI id = k * np + j;
                if (bvs.phaseExist[id]) {
                    perf[p].transj[j] = temp * bvs.kr[id] / bvs.mu[id];
                    perf[p].transINJ += perf[p].transj[j];
                }
            }
            if (ifUseUnweight) {
                perf[p].transINJ = perf[p].WI;
            }
        }
    }
    else {
        for (USI p = 0; p < numPerf; p++) {
            const OCP_USI k    = perf[p].location;
            const OCP_DBL temp = perf[p].WI * perf[p].multiplier;

            // multi phase
            for (USI j = 0; j < np; j++) {
                perf[p].transj[j] = 0;
                OCP_USI id = k * np + j;
                if (bvs.phaseExist[id]) {
                    perf[p].transj[j] = temp * bvs.kr[id] / bvs.mu[id];
                }
            }
        }
    }
}


void PeacemanWell::CalFlux(const Bulk& bk, const OCP_BOOL ReCalXi)
{
    OCP_FUNCNAME;

    const BulkVarSet& bvs = bk.vs;

    fill(qi_lbmol.begin(), qi_lbmol.end(), 0.0);

    if (opt.type == WellType::injector) {

        for (USI p = 0; p < numPerf; p++) {
            const OCP_USI k  = perf[p].location;
            const OCP_DBL dP = bvs.P[k] - perf[p].P;

            perf[p].qt_ft3 = perf[p].transINJ * dP;

            if (ReCalXi) {
                perf[p].xi = bk.PVTm.GetPVT(k)->XiPhase(
                    perf[p].P, opt.injTemp, opt.injZi, opt.injPhase);
            }
            for (USI i = 0; i < nc; i++) {
                perf[p].qi_lbmol[i] = perf[p].qt_ft3 * perf[p].xi * opt.injZi[i];
                qi_lbmol[i] += perf[p].qi_lbmol[i];
            }
        }
    }
    else {

        for (USI p = 0; p < numPerf; p++) {
            const OCP_USI k = perf[p].location;
            perf[p].qt_ft3 = 0;
            fill(perf[p].qi_lbmol.begin(), perf[p].qi_lbmol.end(), 0.0);
            fill(perf[p].qj_ft3.begin(), perf[p].qj_ft3.end(), 0.0);

            for (USI j = 0; j < np; j++) {
                OCP_USI id = k * np + j;
                if (bvs.phaseExist[id]) {
                    OCP_DBL dP = bvs.Pj[id] - perf[p].P;

                    perf[p].qj_ft3[j] = perf[p].transj[j] * dP;
                    perf[p].qt_ft3 += perf[p].qj_ft3[j];

                    OCP_DBL xi = bvs.xi[id];
                    OCP_DBL xij;
                    for (USI i = 0; i < nc; i++) {
                        xij = bvs.xij[id * nc + i];
                        perf[p].qi_lbmol[i] += perf[p].qj_ft3[j] * xi * xij;
                    }
                }
            }
            for (USI i = 0; i < nc; i++) qi_lbmol[i] += perf[p].qi_lbmol[i];
        }
    }
}


/// Pressure in injection well equals maximum ones in injection well,
/// which is input by users. this function is used to check if operation mode of
/// well shoubld be swtched.
OCP_DBL PeacemanWell::CalInjRateMaxBHP(const Bulk& bk)
{
    OCP_FUNCNAME;

    const BulkVarSet& bvs = bk.vs;

    OCP_DBL qj = 0;
    const OCP_DBL Pwell = opt.maxBHP;

    for (USI p = 0; p < numPerf; p++) {

        const OCP_DBL Pperf = Pwell + dG[p];
        const OCP_USI k = perf[p].location;
        const OCP_DBL xi = bk.PVTm.GetPVT(k)->XiPhase(Pperf, opt.injTemp, opt.injZi, opt.injPhase);
        const OCP_DBL dP = Pperf - bvs.P[k];
        qj += perf[p].transINJ * xi * dP;
    }

    const OCP_DBL fac = UnitConvertR2S(opt.injPhase, mixture->CalVmStd(Psurf, Tsurf, &opt.injZi[0], opt.injPhase));
    return qj * fac;
}

/// Pressure in production well equals minial ones in production well,
/// which is input by users. this function is used to check if operation mode of
/// well shoubld be swtched.
OCP_DBL PeacemanWell::CalProdRateMinBHP(const Bulk& bk)
{
    OCP_FUNCNAME;

    const BulkVarSet& bvs = bk.vs;

    OCP_DBL qj = 0;
    const OCP_DBL Pwell = opt.minBHP;

    vector<OCP_DBL> tmpQi_lbmol(nc, 0);
    vector<OCP_DBL> tmpQj(np, 0);

    for (USI p = 0; p < numPerf; p++) {

        OCP_DBL Pperf = Pwell + dG[p];
        OCP_USI k = perf[p].location;

        for (USI j = 0; j < np; j++) {
            OCP_USI id = k * np + j;
            if (bvs.phaseExist[id]) {
                OCP_DBL dP = bvs.Pj[id] - Pperf;
                OCP_DBL temp = perf[p].transj[j] * bvs.xi[id] * dP;
                for (USI i = 0; i < nc; i++) {
                    tmpQi_lbmol[i] += bvs.xij[id * nc + i] * temp;
                }
            }
        }
    }

    qj = 0;
    mixture->CalVStd(Psurf, Tsurf, &tmpQi_lbmol[0]);
    for (USI j = 0; j < np; j++) {
        qj += UnitConvertR2S(j, mixture->GetVarSet().vj[j]) * opt.prodPhaseWeight[j];
    }

    return qj;
}


void PeacemanWell::CalInjQj(const Bulk& bvs, const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    OCP_DBL qj = 0;

    for (USI i = 0; i < nc; i++) {
        qj += qi_lbmol[i];
    }
    const OCP_DBL fac = UnitConvertR2S(opt.injPhase, mixture->CalVmStd(Psurf, Tsurf, &opt.injZi[0], opt.injPhase));

    if (opt.injPhase == PhaseType::wat) {
        WWIR = -qj * fac;
        WWIT += WWIR * dt;
    }
    else {
        WGIR = -qj * fac;
        WGIT += WGIR * dt;
    }
}

void PeacemanWell::CalProdQj(const Bulk& bvs, const OCP_DBL& dt)
{
    mixture->CalVStd(Psurf, Tsurf, &qi_lbmol[0]);
    const auto o = mixture->OilIndex();
    const auto g = mixture->GasIndex();
    const auto w = mixture->WatIndex();

    if (o >= 0) WOPR = UnitConvertR2S(o, mixture->GetVarSet().vj[o]);
    if (g >= 0) WGPR = UnitConvertR2S(g, mixture->GetVarSet().vj[g]);
    if (w >= 0) WWPR = UnitConvertR2S(w, mixture->GetVarSet().vj[w]);

    WOPT += WOPR * dt;
    WGPT += WGPR * dt;
    WWPT += WWPR * dt;
}


ReservoirState PeacemanWell::CheckCrossFlow(const Bulk& bk)
{
    OCP_FUNCNAME;

    OCP_USI  k;
    OCP_BOOL flagC = OCP_TRUE;

    const BulkVarSet& bvs = bk.vs;

    if (opt.type == WellType::productor) {
        for (USI p = 0; p < numPerf; p++) {
            k = perf[p].location;
            OCP_DBL minP = bvs.P[k];
            if (perf[p].state == WellState::open && minP < perf[p].P) {
                cout << std::left << std::setw(12) << name << "  "
                    << "Well P = " << perf[p].P << ", "
                    << "Bulk P = " << minP << endl;
                perf[p].state = WellState::close;
                perf[p].multiplier = 0;
                flagC = OCP_FALSE;
                break;
            }
            else if (perf[p].state == WellState::close && minP > perf[p].P) {
                perf[p].state = WellState::open;
                perf[p].multiplier = 1;
            }
        }
    }
    else {
        for (USI p = 0; p < numPerf; p++) {
            k = perf[p].location;
            if (perf[p].state == WellState::open && bvs.P[k] > perf[p].P) {
                cout << std::left << std::setw(12) << name << "  "
                    << "Well P = " << perf[p].P << ", "
                    << "Bulk P = " << bvs.P[k] << endl;
                perf[p].state = WellState::close;
                perf[p].multiplier = 0;
                flagC = OCP_FALSE;
                break;
            }
            else if (perf[p].state == WellState::close && bvs.P[k] < perf[p].P) {
                perf[p].state = WellState::open;
                perf[p].multiplier = 1;
            }
        }
    }

    OCP_BOOL flag = OCP_FALSE;
    // check well --  if all perf are closed, open the depthest perf
    for (USI p = 0; p < numPerf; p++) {
        if (perf[p].state == WellState::open) {
            flag = OCP_TRUE;
            break;
        }
    }

    if (!flag) {
        // open the deepest perf
        perf.back().state = WellState::open;
        perf.back().multiplier = 1;
        cout << "### WARNING: All perfs of " << name
            << " are closed! Open the last perf!\n";
    }

    if (!flagC) {
        // if crossflow happens, then corresponding perforation will be closed,
        // the multiplier of perforation will be set to zero, so trans of well
        // should be recalculated!
        //
        // dG = ldG;
        CalTrans(bk);
        // CalFlux(bvs);
        // CaldG(bvs);
        // CheckOptMode(bvs);
        return ReservoirState::well_crossflow;
    }

    return ReservoirState::well_success;
}


void PeacemanWell::CalFactor(const Bulk& bk) const
{
    if (opt.mode == WellOptMode::bhp)  return;

    const BulkVarSet& bvs = bk.vs;

    if (opt.type == WellType::productor) {
        if (mixture->IfWellFriend()) {
            // For black oil models -- phase = components
            vector<OCP_DBL> qitmp(nc, 1.0);
            mixture->CalVStd(Psurf, Tsurf, &qitmp[0]);
            for (USI i = 0; i < nc; i++) {
                factor[i] = UnitConvertR2S(i, mixture->GetVarSet().vj[i]) * opt.prodPhaseWeight[i];
            }
        }
        else {
            // For other models
            vector<OCP_DBL> qitmp(nc, 0);
            OCP_DBL         qt = 0;
            OCP_BOOL        flag = OCP_TRUE;
            for (USI i = 0; i < nc; i++) {
                qt += qi_lbmol[i];
                if (qi_lbmol[i] < 0) flag = OCP_FALSE;
            }
            if (qt > TINY && flag) {
                qitmp = qi_lbmol;
            }
            else {
                for (USI p = 0; p < numPerf; p++) {
                    OCP_USI n = perf[p].location;

                    for (USI j = 0; j < np; j++) {
                        const OCP_USI n_np_j = n * np + j;
                        if (!bvs.phaseExist[n_np_j]) continue;
                        for (USI k = 0; k < nc; k++) {
                            qitmp[k] += perf[p].transj[j] * bvs.xi[n_np_j] *
                                bvs.xij[n_np_j * nc + k];
                        }
                    }
                }
            }

            qt = 0;
            for (USI i = 0; i < nc; i++) qt += qitmp[i];
            OCP_DBL qv = 0;
            mixture->CalVStd(Psurf, Tsurf, &qitmp[0]);
            for (USI j = 0; j < np; j++) {
                qv += UnitConvertR2S(j, mixture->GetVarSet().vj[j]) * opt.prodPhaseWeight[j];
            }
            fill(factor.begin(), factor.end(), qv / qt);
            if (factor[0] < 1E-12) {
                OCP_ABORT("Wrong Condition!");
            }
        }
    }
    else if (opt.type == WellType::injector) {
        fill(factor.begin(), factor.end(), UnitConvertR2S(opt.injPhase, mixture->CalVmStd(Psurf, Tsurf, &opt.injZi[0], opt.injPhase)));
    }
    else {
        OCP_ABORT("WRONG Well Type!");
    }
}


/// It calculates pressure difference between perforations iteratively.
/// This function can be used in both black oil model and compositional model.
/// stability of this method shoule be tested.
void PeacemanWell::CaldG(const Bulk& bk)
{
    OCP_FUNCNAME;

    if (opt.state == WellState::open) {
        if (opt.type == WellType::injector) CalInjdG(bk);
        else                                CalProddG01(bk);
        CalPerfP();
    }
}


void PeacemanWell::CalInjdG(const Bulk& bk)
{
    OCP_FUNCNAME;

    const OCP_DBL   maxlen = 10;
    USI             seg_num = 0;
    OCP_DBL         seg_len = 0;
    vector<OCP_DBL> dGperf(numPerf, 0);

    if (depth <= perf.front().depth) {
        // Well is higher
        for (OCP_INT p = numPerf - 1; p >= 0; p--) {
            if (p == 0) {
                seg_num = ceil(fabs((perf[0].depth - depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[0].depth - depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p].depth - perf[p - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p].depth - perf[p - 1].depth) / seg_num;
            }
            OCP_USI n = perf[p].location;
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;
            auto    PVT = bk.PVTm.GetPVT(n);

            for (USI i = 0; i < seg_num; i++) {
                Ptmp -= PVT->RhoPhase(Ptmp, 0, opt.injTemp, opt.injZi, opt.injPhase) * GRAVITY_FACTOR * seg_len;
            }
            dGperf[p] = Pperf - Ptmp;
        }
        dG[0] = dGperf[0];
        for (USI p = 1; p < numPerf; p++) {
            dG[p] = dG[p - 1] + dGperf[p];
        }
    }
    else if (depth >= perf[numPerf - 1].depth) {
        // Well is lower
        for (USI p = 0; p < numPerf; p++) {
            if (p == numPerf - 1) {
                seg_num = ceil(fabs((depth - perf[numPerf - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (depth - perf[numPerf - 1].depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p + 1].depth - perf[p].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p + 1].depth - perf[p].depth) / seg_num;
            }
            OCP_USI n = perf[p].location;
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;

            auto PVT = bk.PVTm.GetPVT(n);
            for (USI i = 0; i < seg_num; i++) {
                Ptmp += PVT->RhoPhase(Ptmp, 0, opt.injTemp, opt.injZi, opt.injPhase) * GRAVITY_FACTOR * seg_len;
            }
            dGperf[p] = Ptmp - Pperf;
        }
        dG[numPerf - 1] = dGperf[numPerf - 1];
        for (OCP_INT p = numPerf - 2; p >= 0; p--) {
            dG[p] = dG[p + 1] + dGperf[p];
        }
    }
}


// Use transj
void PeacemanWell::CalProddG01(const Bulk& bk)
{
    OCP_FUNCNAME;

    const BulkVarSet& bvs = bk.vs;

    const OCP_DBL   maxlen = 10;
    USI             seg_num = 0;
    OCP_DBL         seg_len = 0;
    vector<OCP_DBL> dGperf(numPerf, 0);
    vector<OCP_DBL> tmpNi(nc, 0);
    OCP_DBL         rhotmp, qtacc, rhoacc;

    if (depth <= perf.front().depth) {
        // Well is higher
        for (OCP_INT p = numPerf - 1; p >= 0; p--) {
            if (p == 0) {
                seg_num = ceil(fabs((perf[0].depth - depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[0].depth - depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p].depth - perf[p - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p].depth - perf[p - 1].depth) / seg_num;
            }
            OCP_USI n = perf[p].location;
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;

            // fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (USI j = 0; j < np; j++) {
                const OCP_USI n_np_j = n * np + j;
                if (!bvs.phaseExist[n_np_j]) continue;
                for (USI k = 0; k < nc; k++) {
                    tmpNi[k] += (bvs.P[n] - perf[p].P) * perf[p].transj[j] *
                        bvs.xi[n_np_j] * bvs.xij[n_np_j * nc + k];
                }
            }
            OCP_DBL tmpSum = Dnorm1(nc, &tmpNi[0]);
            if (tmpSum < TINY) {
                for (USI i = 0; i < nc; i++) {
                    tmpNi[i] = bvs.Ni[n * nc + i];
                }
            }

            auto PVT = bk.PVTm.GetPVT(n);
            for (USI i = 0; i < seg_num; i++) {
                qtacc = rhoacc = 0;
                PVT->Flash(Ptmp, bvs.T[n], tmpNi.data());
                for (USI j = 0; j < np; j++) {
                    if (PVT->GetPhaseExist(j)) {
                        rhotmp = PVT->GetRho(j);
                        qtacc += PVT->GetVj(j);
                        rhoacc += PVT->GetVj(j) * rhotmp;
                    }
                }
                Ptmp -= rhoacc / qtacc * seg_len * GRAVITY_FACTOR;
            }
            dGperf[p] = Pperf - Ptmp;
        }
        dG[0] = dGperf[0];
        for (USI p = 1; p < numPerf; p++) {
            dG[p] = dG[p - 1] + dGperf[p];
        }
    }
    else if (depth >= perf[numPerf - 1].depth) {
        // Well is lower
        for (USI p = 0; p < numPerf; p++) {
            if (p == numPerf - 1) {
                seg_num = ceil(fabs((depth - perf[numPerf - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (depth - perf[numPerf - 1].depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p + 1].depth - perf[p].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p + 1].depth - perf[p].depth) / seg_num;
            }
            OCP_USI n = perf[p].location;
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;

            // fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (USI j = 0; j < np; j++) {
                const OCP_USI n_np_j = n * np + j;
                if (!bvs.phaseExist[n_np_j]) continue;
                for (USI k = 0; k < nc; k++) {
                    tmpNi[k] += (bvs.P[n] - perf[p].P) * perf[p].transj[j] *
                        bvs.xi[n_np_j] * bvs.xij[n_np_j * nc + k];
                }
            }
            OCP_DBL tmpSum = Dnorm1(nc, &tmpNi[0]);
            if (tmpSum < TINY) {
                for (USI i = 0; i < nc; i++) {
                    tmpNi[i] = bvs.Ni[n * nc + i];
                }
            }

            auto PVT = bk.PVTm.GetPVT(n);
            for (USI i = 0; i < seg_num; i++) {
                qtacc = rhoacc = 0;
                PVT->Flash(Ptmp, bvs.T[n], tmpNi.data());
                for (USI j = 0; j < np; j++) {
                    if (PVT->GetPhaseExist(j)) {
                        rhotmp = PVT->GetRho(j);
                        qtacc += PVT->GetVj(j);
                        rhoacc += PVT->GetVj(j) * rhotmp;
                    }
                }
                Ptmp += rhoacc / qtacc * seg_len * GRAVITY_FACTOR;
            }
            dGperf[p] = Ptmp - Pperf;
        }
        dG[numPerf - 1] = dGperf[numPerf - 1];
        for (OCP_INT p = numPerf - 2; p >= 0; p--) {
            dG[p] = dG[p + 1] + dGperf[p];
        }
    }
}


// Use bulk
void PeacemanWell::CalProddG02(const Bulk& bk)
{
    OCP_FUNCNAME;

    const BulkVarSet& bvs = bk.vs;

    const OCP_DBL   maxlen = 10;
    USI             seg_num = 0;
    OCP_DBL         seg_len = 0;
    vector<OCP_DBL> dGperf(numPerf, 0);
    vector<OCP_DBL> tmpNi(nc, 0);
    OCP_DBL         rhotmp, qtacc, rhoacc;

    if (depth <= perf.front().depth) {
        // Well is higher
        for (OCP_INT p = numPerf - 1; p >= 0; p--) {
            if (p == 0) {
                seg_num = ceil(fabs((perf[0].depth - depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[0].depth - depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p].depth - perf[p - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p].depth - perf[p - 1].depth) / seg_num;
            }
            OCP_USI n = perf[p].location;
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;

            fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (USI j = 0; j < np; j++) {
                const OCP_USI n_np_j = n * np + j;
                if (!bvs.phaseExist[n_np_j]) continue;
                for (USI k = 0; k < nc; k++) {
                    tmpNi[k] += (perf[p].transj[j] > 0) * bvs.xi[n_np_j] *
                        bvs.xij[n_np_j * nc + k];
                }
            }
            OCP_DBL tmpSum = Dnorm1(nc, &tmpNi[0]);
            if (tmpSum < TINY) {
                for (USI i = 0; i < nc; i++) {
                    tmpNi[i] = bvs.Ni[n * nc + i];
                }
            }

            auto PVT = bk.PVTm.GetPVT(n);
            for (USI i = 0; i < seg_num; i++) {
                qtacc = rhoacc = 0;
                PVT->Flash(Ptmp, bvs.T[n], tmpNi.data());
                for (USI j = 0; j < np; j++) {
                    if (PVT->GetPhaseExist(j)) {
                        rhotmp = PVT->GetRho(j);
                        qtacc += PVT->GetVj(j);
                        rhoacc += PVT->GetVj(j) * rhotmp;
                    }
                }
                Ptmp -= rhoacc / qtacc * seg_len * GRAVITY_FACTOR;
            }
            dGperf[p] = Pperf - Ptmp;
        }
        dG[0] = dGperf[0];
        for (USI p = 1; p < numPerf; p++) {
            dG[p] = dG[p - 1] + dGperf[p];
        }
    }
    else if (depth >= perf[numPerf - 1].depth) {
        // Well is lower
        for (USI p = 0; p < numPerf; p++) {
            if (p == numPerf - 1) {
                seg_num = ceil(fabs((depth - perf[numPerf - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (depth - perf[numPerf - 1].depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p + 1].depth - perf[p].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p + 1].depth - perf[p].depth) / seg_num;
            }
            OCP_USI n = perf[p].location;
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;

            fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (USI j = 0; j < np; j++) {
                const OCP_USI n_np_j = n * np + j;
                if (!bvs.phaseExist[n_np_j]) continue;
                for (USI k = 0; k < nc; k++) {
                    tmpNi[k] += (perf[p].transj[j] > 0) * bvs.xi[n_np_j] *
                        bvs.xij[n_np_j * nc + k];
                }
            }
            OCP_DBL tmpSum = Dnorm1(nc, &tmpNi[0]);
            if (tmpSum < TINY) {
                for (USI i = 0; i < nc; i++) {
                    tmpNi[i] = bvs.Ni[n * nc + i];
                }
            }

            auto PVT = bk.PVTm.GetPVT(n);
            for (USI i = 0; i < seg_num; i++) {
                qtacc = rhoacc = 0;
                PVT->Flash(Ptmp, bvs.T[n], tmpNi.data());
                for (USI j = 0; j < np; j++) {
                    if (PVT->GetPhaseExist(j)) {
                        rhotmp = PVT->GetRho(j);
                        qtacc += PVT->GetVj(j);
                        rhoacc += PVT->GetVj(j) * rhotmp;
                    }
                }
                Ptmp += rhoacc / qtacc * seg_len * GRAVITY_FACTOR;
            }
            dGperf[p] = Ptmp - Pperf;
        }
        dG[numPerf - 1] = dGperf[numPerf - 1];
        for (OCP_INT p = numPerf - 2; p >= 0; p--) {
            dG[p] = dG[p + 1] + dGperf[p];
        }
    }
}


// Use qi_lbmol
void PeacemanWell::CalProddG(const Bulk& bk)
{
    OCP_FUNCNAME;

    const BulkVarSet& bvs = bk.vs;

    const OCP_DBL   maxlen = 5;
    USI             seg_num = 0;
    OCP_DBL         seg_len = 0;
    vector<OCP_DBL> tmpNi(nc, 0);
    vector<OCP_DBL> dGperf(numPerf, 0);
    OCP_DBL         qtacc = 0;
    OCP_DBL         rhoacc = 0;
    OCP_DBL         rhotmp = 0;

    if (depth <= perf.front().depth) {
        // Well is higher

        // check qi_lbmol   ----   test
        if (perf[numPerf - 1].state == WellState::close) {
            for (OCP_INT p = numPerf - 2; p >= 0; p--) {
                if (perf[p].state == WellState::open) {
                    for (USI i = 0; i < nc; i++) {
                        perf[numPerf - 1].qi_lbmol[i] = perf[p].qi_lbmol[i];
                    }
                    break;
                }
            }
        }

        for (OCP_INT p = numPerf - 1; p >= 0; p--) {

            if (p == 0) {
                seg_num = ceil(fabs((perf[0].depth - depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[0].depth - depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p].depth - perf[p - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p].depth - perf[p - 1].depth) / seg_num;
            }

            OCP_USI n = perf[p].location;
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;

            for (USI i = 0; i < nc; i++) {
                tmpNi[i] += perf[p].qi_lbmol[i];
            }

            // check tmpNi
            for (auto& v : tmpNi) {
                v = fabs(v);
            }

            auto PVT = bk.PVTm.GetPVT(n);
            for (USI k = 0; k < seg_num; k++) {
                PVT->Flash(Ptmp, bvs.T[n], tmpNi.data());
                for (USI j = 0; j < np; j++) {
                    if (PVT->GetPhaseExist(j)) {
                        rhotmp = PVT->GetRho(j);
                        qtacc += PVT->GetVj(j) / seg_num;
                        rhoacc += PVT->GetVj(j) * rhotmp / seg_num;
#ifdef DEBUG
                        if (rhotmp <= 0 || !isfinite(rhotmp)) {
                            OCP_ABORT("Wrong rho " + to_string(rhotmp));
                        }
#endif // DEBUG
                    }
                }
                Ptmp -= rhoacc / qtacc * seg_len * GRAVITY_FACTOR;
            }
            dGperf[p] = Pperf - Ptmp;
        }
        dG[0] = dGperf[0];
        for (USI p = 1; p < numPerf; p++) {
            dG[p] = dG[p - 1] + dGperf[p];
        }
    }
    else if (depth >= perf.back().depth) {
        // Well is lower

        // check qi_lbmol   ----   test
        if (perf[0].state == WellState::close) {
            for (USI p = 1; p <= numPerf; p++) {
                if (perf[p].state == WellState::open) {
                    for (USI i = 0; i < nc; i++) {
                        perf[numPerf - 1].qi_lbmol[i] = perf[p].qi_lbmol[i];
                    }
                    break;
                }
            }
        }

        for (USI p = 0; p < numPerf; p++) {
            if (p == numPerf - 1) {
                seg_num = ceil(fabs((depth - perf[numPerf - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (depth - perf[numPerf - 1].depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p + 1].depth - perf[p].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p + 1].depth - perf[p].depth) / seg_num;
            }

            OCP_USI n = perf[p].location;
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;


            fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (OCP_INT p1 = numPerf - 1; p1 - p >= 0; p1--) {
                for (USI i = 0; i < nc; i++) {
                    tmpNi[i] += perf[p1].qi_lbmol[i];
                }
            }

            // check tmpNi
            for (auto& v : tmpNi) {
                v = fabs(v);
            }

            auto PVT = bk.PVTm.GetPVT(n);
            for (USI k = 0; k < seg_num; k++) {
                PVT->Flash(Ptmp, bvs.T[n], tmpNi.data());
                for (USI j = 0; j < np; j++) {
                    if (PVT->GetPhaseExist(j)) {
                        rhotmp = PVT->GetRho(j);
                        qtacc += PVT->GetVj(j) / seg_num;
                        rhoacc += PVT->GetVj(j) * rhotmp / seg_num;
                    }
                }
                Ptmp += rhoacc / qtacc * seg_len * GRAVITY_FACTOR;
            }
            dGperf[p] = Ptmp - Pperf;
        }
        dG[numPerf - 1] = dGperf[numPerf - 1];
        for (OCP_INT p = numPerf - 2; p >= 0; p--) {
            dG[p] = dG[p + 1] + dGperf[p];
        }

    }
    else {
        OCP_ABORT("Wrong well position!");
    }
}


void PeacemanWellIsoT::GetSolutionFIM(const vector<OCP_DBL>& u, OCP_USI& wId)
{
    if (opt.state != WellState::open)  return;
    bhp += u[wId];
    CalPerfP();
    wId += nc + 1;
}


void PeacemanWellIsoT::GetSolutionIMPEC(const vector<OCP_DBL>& u, OCP_USI& wId)
{
    if (opt.state != WellState::open)  return;
    bhp = u[wId];
    CalPerfP();
    wId++;
}


void PeacemanWellIsoT::AssembleMatFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{
    if (opt.state == WellState::open) {
        if (opt.type == WellType::injector)        AssembleMatInjFIM(ls, bk, dt);
        else if (opt.type == WellType::productor)  AssembleMatProdFIM(ls, bk, dt);
        else                                       OCP_ABORT("WRONG Well Type!");
    }
}


void PeacemanWellIsoT::CalResFIM(OCP_USI& wId, OCPNRresidual& res, const Bulk& bk, const OCP_DBL& dt) const
{
    if (opt.state == WellState::open) {

        const USI len = nc + 1;
        // Well to Bulk
        for (USI p = 0; p < numPerf; p++) {
            const OCP_USI k = perf[p].location;
            for (USI i = 0; i < nc; i++) {
                res.resAbs[k * len + 1 + i] += perf[p].qi_lbmol[i] * dt;
            }
        }
        // Well Self
        switch (opt.mode)
        {
        case WellOptMode::bhp:
            res.resAbs[wId] = bhp - opt.tarBHP;
            break;
        case WellOptMode::irate:
        case WellOptMode::orate:
        case WellOptMode::grate:
        case WellOptMode::wrate:
        case WellOptMode::lrate:
            CalFactor(bk);
            if (opt.type == WellType::injector)  res.resAbs[wId] = opt.tarRate;
            else                                 res.resAbs[wId] = -opt.tarRate;
            for (USI i = 0; i < nc; i++) {
                res.resAbs[wId] += qi_lbmol[i] * factor[i];
            }
            res.maxWellRelRes_mol =
                max(res.maxWellRelRes_mol,
                    fabs(res.resAbs[wId] / opt.tarRate));
            break;
        default:
            OCP_ABORT("Wrong well opt mode!");
            break;
        }
        wId += len;
    }
}


void PeacemanWellIsoT::AssembleMatInjFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{
    const BulkVarSet& bvs = bk.vs;

    const USI ncol   = nc + 1;
    const USI ncol2  = np * nc + np;
    const USI bsize  = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
    OCP_USI   n_np_j;

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);

    OCP_DBL mu, muP, dP, transIJ;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, bmat);

    CalFactor(bk);

    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        for (USI j = 0; j < np; j++) {
            n_np_j = n * np + j;
            if (!bvs.phaseExist[n_np_j]) continue;

            mu = bvs.mu[n_np_j];
            muP = bvs.muP[n_np_j];

            dP = bvs.P[n] - bhp - dG[p];

            for (USI i = 0; i < nc; i++) {
                // dQ / dP
                transIJ = perf[p].transj[j] * perf[p].xi * opt.injZi[i];
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu);
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                for (USI k = 0; k < np; k++) {
                    dQdXsB[(i + 1) * ncol2 + k] +=
                        perf[p].WI * perf[p].multiplier * perf[p].xi *
                        opt.injZi[i] * bvs.dKrdS[n_np_j * np + k] * dP / mu;
                }
                // dQ / dxij
                for (USI k = 0; k < nc; k++) {
                    dQdXsB[(i + 1) * ncol2 + np + j * nc + k] +=
                        -transIJ * dP / mu * bvs.mux[n_np_j * nc + k];
                }
            }
        }

        // Bulk to Well
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bvs.dSec_dPri[n * bsize2], 1,
            bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Bulk - Bulk -- add
        ls.AddDiag(n, bmat);

        // Bulk - Well -- insert
        bmat = dQdXpW;
        Dscalar(bsize, dt, bmat.data());
        ls.NewOffDiag(n, wId, bmat);

        // Well
        switch (opt.mode) {
        case WellOptMode::irate:
        case WellOptMode::orate:
        case WellOptMode::grate:
        case WellOptMode::wrate:
        case WellOptMode::lrate:
            // Well - Well -- add
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < nc; i++) {
                bmat[0] += dQdXpW[(i + 1) * ncol] * factor[i];
                bmat[(i + 1) * ncol + i + 1] = 1;
            }
            ls.AddDiag(wId, bmat);

            // Well - Bulk -- insert
            bmat = dQdXpB;
            DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bvs.dSec_dPri[n * bsize2],
                1, bmat.data());
            fill(bmat2.begin(), bmat2.end(), 0.0);
            for (USI i = 0; i < nc; i++) {
                Daxpy(ncol, factor[i], bmat.data() + (i + 1) * ncol, bmat2.data());
            }
            ls.NewOffDiag(wId, n, bmat2);
            break;

        case WellOptMode::bhp:
            // Well - Well -- add
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < ncol; i++) {
                bmat[i * ncol + i] = 1;
            }
            ls.AddDiag(wId, bmat);

            // Well - Bulk -- insert
            fill(bmat.begin(), bmat.end(), 0.0);
            ls.NewOffDiag(wId, n, bmat);
            break;

        default:
            OCP_ABORT("Wrong Well Opt mode!");
            break;
        }
    }
}


void PeacemanWellIsoT::AssembleMatProdFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{
    const BulkVarSet& bvs = bk.vs;

    const USI ncol = nc + 1;
    const USI ncol2 = np * nc + np;
    const USI bsize = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
    OCP_USI   n_np_j;

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);

    OCP_DBL xij, xi, mu, muP, xiP, dP, transIJ, tmp;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, bmat);

    // Set Commponent Weight
    CalFactor(bk);

    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        for (USI j = 0; j < np; j++) {
            n_np_j = n * np + j;
            if (!bvs.phaseExist[n_np_j]) continue;

            dP = bvs.Pj[n_np_j] - bhp - dG[p];
            xi = bvs.xi[n_np_j];
            mu = bvs.mu[n_np_j];
            muP = bvs.muP[n_np_j];
            xiP = bvs.xiP[n_np_j];

            for (USI i = 0; i < nc; i++) {
                xij = bvs.xij[n_np_j * nc + i];
                // dQ / dP
                transIJ = perf[p].transj[j] * xi * xij;
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu) +
                    dP * perf[p].transj[j] * xij * xiP;
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                for (USI k = 0; k < np; k++) {
                    tmp = perf[p].WI * perf[p].multiplier * dP / mu * xi *
                        xij * bvs.dKrdS[n_np_j * np + k];
                    // capillary pressure
                    tmp += transIJ * bvs.dPcdS[n_np_j * np + k];
                    dQdXsB[(i + 1) * ncol2 + k] += tmp;
                }
                // dQ / dCij
                for (USI k = 0; k < nc; k++) {
                    tmp = dP * perf[p].transj[j] * xij *
                        (bvs.xix[n_np_j * nc + k] - xi / mu * bvs.mux[n_np_j * nc + k]);
                    dQdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                }
                dQdXsB[(i + 1) * ncol2 + np + j * nc + i] +=
                    perf[p].transj[j] * xi * dP;
            }
        }

        // Bulk - Bulk -- add
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bvs.dSec_dPri[n * bsize2], 1,
            bmat.data());

        Dscalar(bsize, dt, bmat.data());
        ls.AddDiag(n, bmat);

        // Bulk - Well -- insert
        bmat = dQdXpW;
        Dscalar(bsize, dt, bmat.data());
        ls.NewOffDiag(n, wId, bmat);

        // Well
        switch (opt.mode) {
        case WellOptMode::irate:
        case WellOptMode::orate:
        case WellOptMode::grate:
        case WellOptMode::wrate:
        case WellOptMode::lrate:
            // Well - Well -- add
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < nc; i++) {
                bmat[0] += dQdXpW[(i + 1) * ncol] * factor[i];
                bmat[(i + 1) * ncol + i + 1] = 1;
            }
            ls.AddDiag(wId, bmat);

            // Well - Bulk -- insert
            bmat = dQdXpB;
            DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bvs.dSec_dPri[n * bsize2],
                1, bmat.data());
            fill(bmat2.begin(), bmat2.end(), 0.0);
            for (USI i = 0; i < nc; i++) {
                Daxpy(ncol, factor[i], bmat.data() + (i + 1) * ncol,
                    bmat2.data());
            }
            ls.NewOffDiag(wId, n, bmat2);
            break;

        case WellOptMode::bhp:
            // Well - Well -- add
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < ncol; i++) {
                bmat[i * ncol + i] = 1;
            }
            ls.AddDiag(wId, bmat);

            // Well - Bulk -- insert
            fill(bmat.begin(), bmat.end(), 0.0);
            ls.NewOffDiag(wId, n, bmat);
            break;

        default:
            OCP_ABORT("Wrong Well Opt mode!");
            break;
        }
    }
}


void PeacemanWellIsoT::AssembleMatIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{
    if (opt.state == WellState::open) {
        if (opt.type == WellType::injector)        AssembleMatInjIMPEC(ls, bk, dt);
        else if (opt.type == WellType::productor)  AssembleMatProdIMPEC(ls, bk, dt);
        else                                       OCP_ABORT("WRONG Well Type!");
    }
}


void PeacemanWellIsoT::AssembleMatInjIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{
    const BulkVarSet& bvs = bk.vs;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, 0.0);
    CalFactor(bk);
    OCP_DBL Vfi_zi, valb, valw, bb, bw;

    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;

        Vfi_zi = 0;
        for (USI i = 0; i < nc; i++) {
            Vfi_zi += bvs.vfi[n * nc + i] * opt.injZi[i];
        }

        valw = dt * perf[p].xi * perf[p].transINJ;
        bw = valw * dG[p];
        valb = valw * Vfi_zi;
        bb = valb * dG[p];

        // Bulk to Well
        ls.AddDiag(n, valb);
        ls.NewOffDiag(n, wId, -valb);
        ls.AddRhs(n, bb);

        // Well to Bulk
        switch (opt.mode) {
        case WellOptMode::irate:
        case WellOptMode::orate:
        case WellOptMode::grate:
        case WellOptMode::wrate:
        case WellOptMode::lrate:
            ls.AddDiag(wId, valw * factor[0]);
            ls.NewOffDiag(wId, n, -valw * factor[0]);
            ls.AddRhs(wId, -bw * factor[0]);
            break;
        case WellOptMode::bhp:
            ls.NewOffDiag(wId, n, 0);
            break;
        default:
            OCP_ABORT("Wrong well option mode!");
        }
    }

    // Well Self
    switch (opt.mode) {
    case WellOptMode::irate:
    case WellOptMode::orate:
    case WellOptMode::grate:
    case WellOptMode::wrate:
    case WellOptMode::lrate:
        ls.AddRhs(wId, dt * opt.tarRate);
        break;
    case WellOptMode::bhp:
        ls.AddDiag(wId, dt);
        ls.AddRhs(wId, dt * opt.tarBHP);
        ls.AssignGuess(wId, opt.tarBHP);
        break;
    default:
        OCP_ABORT("Wrong well option mode!");
    }
}


void PeacemanWellIsoT::AssembleMatProdIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{
    const BulkVarSet& bvs = bk.vs;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, 0.0);

    // Set Prod Weight
    CalFactor(bk);

    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;

        OCP_DBL valb = 0;
        OCP_DBL bb = 0;
        OCP_DBL valw = 0;
        OCP_DBL bw = 0;

        for (USI j = 0; j < np; j++) {
            const OCP_USI n_np_j = n * np + j;
            if (!bvs.phaseExist[n_np_j]) continue;

            OCP_DBL tempb = 0;
            OCP_DBL tempw = 0;

            for (USI i = 0; i < nc; i++) {
                tempb += bvs.vfi[n * nc + i] * bvs.xij[n_np_j * nc + i];
                tempw += factor[i] * bvs.xij[n_np_j * nc + i];
            }
            OCP_DBL trans = dt * perf[p].transj[j] * bvs.xi[n_np_j];
            valb += tempb * trans;
            valw += tempw * trans;

            OCP_DBL dP = dG[p] - bvs.Pc[n_np_j];
            bb += tempb * trans * dP;
            bw += tempw * trans * dP;
        }

        // Bulk to Well
        ls.AddDiag(n, valb);
        ls.NewOffDiag(n, wId, -valb);
        ls.AddRhs(n, bb);

        // Well to Bulk
        switch (opt.mode) {
        case WellOptMode::irate:
        case WellOptMode::orate:
        case WellOptMode::grate:
        case WellOptMode::wrate:
        case WellOptMode::lrate:
            ls.AddDiag(wId, -valw);
            ls.NewOffDiag(wId, n, valw);
            ls.AddRhs(wId, bw);
            break;
        case WellOptMode::bhp:
            ls.NewOffDiag(wId, n, 0.0);
            break;
        default:
            OCP_ABORT("Wrong well option mode!");
        }
    }

    // Well Self
    switch (opt.mode) {
    case WellOptMode::irate:
    case WellOptMode::orate:
    case WellOptMode::grate:
    case WellOptMode::wrate:
    case WellOptMode::lrate:
        ls.AddRhs(wId, dt * opt.tarRate);
        break;
    case WellOptMode::bhp:
        ls.AddDiag(wId, dt);
        ls.AddRhs(wId, dt * opt.tarBHP);
        ls.AssignGuess(wId, opt.tarBHP);
        break;
    default:
        OCP_ABORT("Wrong well option mode!");
    }
}


void PeacemanWellT::CalResFIM(OCP_USI& wId, OCPNRresidual& res, const Bulk& bk, const OCP_DBL& dt) const
{
    if (opt.state == WellState::open) {

        const BulkVarSet& bvs = bk.vs;

        const USI len = nc + 2;

        // Well to Bulk
        for (USI p = 0; p < numPerf; p++) {
            const OCP_USI k = perf[p].location;
            // Mass Conservation
            for (USI i = 0; i < nc; i++) {
                res.resAbs[k * len + 1 + i] += perf[p].qi_lbmol[i] * dt;
            }
        }

        // Energy Conservation
        if (opt.type == WellType::injector) {
            const OCP_DBL Hw = bk.PVTm.GetPVT(0)->CalInjWellEnthalpy(opt.injTemp, &opt.injZi[0]);
            for (USI p = 0; p < numPerf; p++) {
                const OCP_USI k = perf[p].location;
                res.resAbs[k * len + 1 + nc] += perf[p].qt_ft3 * perf[p].xi * Hw * dt;
            }
        }
        else {
            for (USI p = 0; p < numPerf; p++) {
                const OCP_USI k = perf[p].location;
                for (USI j = 0; j < np; j++) {
                    res.resAbs[k * len + 1 + nc] += perf[p].qj_ft3[j] *
                        bvs.xi[k * np + j] * bvs.H[k * np + j] * dt;
                }
            }
        }

        switch (opt.mode)
        {
        case WellOptMode::bhp:
            res.resAbs[wId] = bhp - opt.tarBHP;
            break;
        case WellOptMode::irate:
        case WellOptMode::orate:
        case WellOptMode::grate:
        case WellOptMode::wrate:
        case WellOptMode::lrate:
            CalFactor(bk);
            if (opt.type == WellType::injector)  res.resAbs[wId] = opt.tarRate;
            else                                 res.resAbs[wId] = -opt.tarRate;
            for (USI i = 0; i < nc; i++) {
                res.resAbs[wId] += qi_lbmol[i] * factor[i];
            }
            res.maxWellRelRes_mol =
                max(res.maxWellRelRes_mol,
                    fabs(res.resAbs[wId] / opt.tarRate));
            break;
        default:
            OCP_ABORT("Wrong well opt mode!");
            break;
        }
        wId += len;
    }
}


void PeacemanWellT::GetSolutionFIM(const vector<OCP_DBL>& u, OCP_USI& wId)
{
    if (opt.state != WellState::open)  return;
    bhp += u[wId];
    CalPerfP();
    wId += nc + 2;
}


void PeacemanWellT::AssembleMatFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{
    if (opt.state == WellState::open) {
        if (opt.type == WellType::injector)        AssembleMatInjFIM(ls, bk, dt);
        else if (opt.type == WellType::productor)  AssembleMatProdFIM(ls, bk, dt);
        else                                       OCP_ABORT("WRONG Well Type!");
    }
}


void PeacemanWellT::AssembleMatInjFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{

    const BulkVarSet& bvs = bk.vs;

    const USI ncol = nc + 2;
    const USI ncol2 = np * nc + np;
    const USI bsize = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
    OCP_USI   n_np_j;

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);

    OCP_DBL mu, muP, muT, dP, Hw;
    OCP_DBL transJ, transIJ;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, bmat);

    CalFactor(bk);

    Hw = bk.PVTm.GetPVT(0)->CalInjWellEnthalpy(opt.injTemp, &opt.injZi[0]);

    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        dP = bvs.P[n] - bhp - dG[p];

        for (USI j = 0; j < np; j++) {
            n_np_j = n * np + j;
            if (!bvs.phaseExist[n_np_j]) continue;

            mu = bvs.mu[n_np_j];
            muP = bvs.muP[n_np_j];
            muT = bvs.muT[n_np_j];

            for (USI i = 0; i < nc; i++) {

                // Mass Conservation
                if (!ifUseUnweight) {
                    transIJ = perf[p].transj[j] * perf[p].xi * opt.injZi[i];
                    // dQ / dP
                    dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu);
                    dQdXpW[(i + 1) * ncol] += -transIJ;

                    // dQ / dT
                    dQdXpB[(i + 2) * ncol - 1] += transIJ * (-dP * muT / mu);
                    dQdXpW[(i + 2) * ncol - 1] += 0;
                }
                else {
                    // dQ / dP
                    transIJ = perf[p].transINJ * perf[p].xi * opt.injZi[i];
                    dQdXpB[(i + 1) * ncol] += transIJ;
                    dQdXpW[(i + 1) * ncol] += -transIJ;
                }

                if (!ifUseUnweight) {
                    // dQ / dS
                    for (USI k = 0; k < np; k++) {
                        dQdXsB[(i + 1) * ncol2 + k] +=
                            perf[p].WI * perf[p].multiplier * perf[p].xi *
                            opt.injZi[i] * bvs.dKrdS[n_np_j * np + k] * dP / mu;
                    }
                    // dQ / dxij
                    for (USI k = 0; k < nc; k++) {
                        dQdXsB[(i + 1) * ncol2 + np + j * nc + k] +=
                            -transIJ * dP / mu * bvs.mux[n_np_j * nc + k];
                    }
                }
            }

            // Energy Conservation
            if (!ifUseUnweight) {
                transJ = perf[p].transj[j] * perf[p].xi;
                // dQ / dP
                dQdXpB[(nc + 1) * ncol] += transJ * Hw * (1 - dP * muP / mu);
                dQdXpW[(nc + 1) * ncol] += -transJ * Hw;

                // dQ / dT
                dQdXpB[(nc + 2) * ncol - 1] += transJ * Hw * (-dP * muT / mu);
                dQdXpW[(nc + 2) * ncol - 1] += 0;

                // dQ / dS
                for (USI k = 0; k < np; k++) {
                    dQdXsB[(nc + 1) * ncol2 + k] +=
                        perf[p].WI * perf[p].multiplier * perf[p].xi *
                        bvs.dKrdS[n_np_j * np + k] * dP / mu * Hw;
                }
                // dQ / dxij
                for (USI k = 0; k < nc; k++) {
                    dQdXsB[(nc + 1) * ncol2 + np + j * nc + k] +=
                        -transJ * dP / mu * bvs.mux[n_np_j * nc + k] * Hw;
                }
            }
            else {
                transJ = perf[p].transINJ * perf[p].xi;
                // dQ / dP
                dQdXpB[(nc + 1) * ncol] += transJ * Hw;
                dQdXpW[(nc + 1) * ncol] += -transJ * Hw;
            }

            if (ifUseUnweight) break;
        }

        // Bulk to Well
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bvs.dSec_dPri[n * bsize2], 1,
            bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Add
        ls.AddDiag(n, bmat);

        // Insert
        bmat = dQdXpW;
        Dscalar(bsize, dt, bmat.data());
        ls.NewOffDiag(n, wId, bmat);

        // Well
        switch (opt.mode) {
        case WellOptMode::irate:
        case WellOptMode::orate:
        case WellOptMode::grate:
        case WellOptMode::wrate:
        case WellOptMode::lrate:
            // Diag
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < nc; i++) {
                bmat[0] += dQdXpW[(i + 1) * ncol] * factor[i];
                bmat[(i + 1) * ncol + i + 1] = 1;
            }
            bmat[ncol * ncol - 1] = 1;
            ls.AddDiag(wId, bmat);

            // OffDiag
            bmat = dQdXpB;
            DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bvs.dSec_dPri[n * bsize2],
                1, bmat.data());
            fill(bmat2.begin(), bmat2.end(), 0.0);
            for (USI i = 0; i < nc; i++) {
                Daxpy(ncol, factor[i], bmat.data() + (i + 1) * ncol, bmat2.data());
            }
            ls.NewOffDiag(wId, n, bmat2);
            break;

        case WellOptMode::bhp:
            // Diag
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < ncol; i++) {
                bmat[i * ncol + i] = 1;
            }
            // Add
            ls.AddDiag(wId, bmat);
            // OffDiag
            fill(bmat.begin(), bmat.end(), 0.0);
            // Insert
            ls.NewOffDiag(wId, n, bmat);
            break;

        default:
            OCP_ABORT("Wrong Well Opt mode!");
            break;
        }
    }
}


void PeacemanWellT::AssembleMatProdFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{
    const BulkVarSet& bvs = bk.vs;

    const USI ncol = nc + 2;
    const USI ncol2 = np * nc + np;
    const USI bsize = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
    OCP_USI   n_np_j;

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);

    OCP_DBL xij, xi, xiP, xiT, mu, muP, muT, dP, transIJ, transJ, H, HT, Hx, tmp;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, bmat);

    // Set Prod Weight
    CalFactor(bk);

    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        for (USI j = 0; j < np; j++) {
            n_np_j = n * np + j;
            if (!bvs.phaseExist[n_np_j]) continue;

            dP = bvs.Pj[n_np_j] - bhp - dG[p];
            xi = bvs.xi[n_np_j];
            mu = bvs.mu[n_np_j];
            xiP = bvs.xiP[n_np_j];
            xiT = bvs.xiT[n_np_j];
            muP = bvs.muP[n_np_j];
            muT = bvs.muT[n_np_j];
            H = bvs.H[n_np_j];
            HT = bvs.HT[n_np_j];

            // Mass Conservation
            for (USI i = 0; i < nc; i++) {
                xij = bvs.xij[n_np_j * nc + i];
                Hx = bvs.Hx[n_np_j * nc + i];
                // dQ / dP
                transIJ = perf[p].transj[j] * xi * xij;
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu) +
                    dP * perf[p].transj[j] * xij * xiP;
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dT
                dQdXpB[(i + 2) * ncol - 1] +=
                    transIJ * (-dP * muT / mu) + dP * perf[p].transj[j] * xij * xiT;
                dQdXpW[(i + 2) * ncol - 1] += 0;

                // dQ / dS
                for (USI k = 0; k < np; k++) {
                    tmp = perf[p].WI * perf[p].multiplier * dP / mu * xi *
                        xij * bvs.dKrdS[n_np_j * np + k];
                    // capillary pressure
                    tmp += transIJ * bvs.dPcdS[n_np_j * np + k];
                    dQdXsB[(i + 1) * ncol2 + k] += tmp;
                }
                // dQ / dxij
                for (USI k = 0; k < nc; k++) {
                    tmp = dP * perf[p].transj[j] * xij *
                        (bvs.xix[n_np_j * nc + k] - xi / mu * bvs.mux[n_np_j * nc + k]);
                    dQdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                }
                dQdXsB[(i + 1) * ncol2 + np + j * nc + i] +=
                    perf[p].transj[j] * xi * dP;
            }

            // Energy Conservation
            transJ = perf[p].transj[j] * xi;
            // dQ / dP
            dQdXpB[(nc + 1) * ncol] +=
                transJ * (1 - dP * muP / mu) * H + dP * perf[p].transj[j] * xiP * H;
            dQdXpW[(nc + 1) * ncol] += -transJ * H;

            // dQ / dT
            dQdXpB[(nc + 2) * ncol - 1] += transJ * (-dP * muT / mu) * H +
                dP * perf[p].transj[j] * xiT * H +
                transJ * dP * HT;
            dQdXpW[(nc + 2) * ncol - 1] += 0;

            // dQ / dS
            for (USI k = 0; k < np; k++) {
                tmp = perf[p].WI * perf[p].multiplier * dP / mu * xi *
                    bvs.dKrdS[n_np_j * np + k] * H;
                // capillary pressure
                tmp += transJ * bvs.dPcdS[n_np_j * np + k] * H;
                dQdXsB[(nc + 1) * ncol2 + k] += tmp;
            }

            // dQ / dxij
            for (USI k = 0; k < nc; k++) {
                tmp = dP * perf[p].transj[j] *
                    (bvs.xix[n_np_j * nc + k] - xi / mu * bvs.mux[n_np_j * nc + k]) *
                    H +
                    transJ * dP * Hx;
                dQdXsB[(nc + 1) * ncol2 + np + j * nc + k] += tmp;
            }
        }

        // Bulk to Well
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bvs.dSec_dPri[n * bsize2], 1,
            bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Add
        ls.AddDiag(n, bmat);
        // Insert
        bmat = dQdXpW;
        Dscalar(bsize, dt, bmat.data());
        ls.NewOffDiag(n, wId, bmat);

        // Well
        switch (opt.mode) {
        case WellOptMode::irate:
        case WellOptMode::orate:
        case WellOptMode::grate:
        case WellOptMode::wrate:
        case WellOptMode::lrate:
            // Diag
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < nc; i++) {
                bmat[0] += dQdXpW[(i + 1) * ncol] * factor[i];
                bmat[(i + 1) * ncol + i + 1] = 1;
            }
            bmat[ncol * ncol - 1] = 1;
            ls.AddDiag(wId, bmat);

            // OffDiag
            bmat = dQdXpB;
            DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bvs.dSec_dPri[n * bsize2],
                1, bmat.data());
            fill(bmat2.begin(), bmat2.end(), 0.0);
            for (USI i = 0; i < nc; i++) {
                Daxpy(ncol, factor[i], bmat.data() + (i + 1) * ncol,
                    bmat2.data());
            }
            ls.NewOffDiag(wId, n, bmat2);
            break;

        case WellOptMode::bhp:
            // Diag
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < ncol; i++) {
                bmat[i * ncol + i] = 1;
            }
            // Add
            ls.AddDiag(wId, bmat);
            // OffDiag
            fill(bmat.begin(), bmat.end(), 0.0);
            // Insert
            ls.NewOffDiag(wId, n, bmat);
            break;

        default:
            OCP_ABORT("Wrong Well Opt mode!");
            break;
        }
    }
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/17/2023      Create file                          */
/*----------------------------------------------------------------------------*/