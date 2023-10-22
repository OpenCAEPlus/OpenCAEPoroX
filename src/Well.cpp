/*! \file    Well.cpp
 *  \brief   Well class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "Well.hpp"
#include <cmath>


void Well::SetupUnit()
{
    if (mixture->OilIndex() >= 0) unitConvert.push_back(1.0 / CONV1);
    if (mixture->GasIndex() >= 0) unitConvert.push_back(1.0 / CONV2);
    if (mixture->WatIndex() >= 0) unitConvert.push_back(1.0 / CONV1);

    OCP_ASSERT(unitConvert.size() == np, "WRONG unitConvert");
}


OCP_DBL Well::UnitConvertR2S(const PhaseType& pt, const OCP_DBL& val) const
{
    if (opt.injPhase == PhaseType::oil)  return (val * unitConvert[mixture->OilIndex()]);
    if (opt.injPhase == PhaseType::gas)  return (val * unitConvert[mixture->GasIndex()]);
    if (opt.injPhase == PhaseType::wat)  return (val * unitConvert[mixture->WatIndex()]);
    OCP_ABORT("WRONG Phase Type!");
}


OCP_DBL Well::UnitConvertR2S(const USI& j, const OCP_DBL& val) const
{
    OCP_ASSERT(j < np, "WRONG PHASE INDEX!");
    return val * unitConvert[j];
}


OCP_BOOL Well::ApplyOpt(const USI& i)
{
    opt = optSet[i];
    if (i > 0 && opt != optSet[i - 1]) return OCP_TRUE;
    else                               return OCP_FALSE;
}


void Well::SetupOpts(const vector<SolventINJ>& sols)
{
    for (auto& opt : optSet) {
        if (opt.state == WellState::close) continue;
        if (opt.type == WellType::injector) {
            SetupOptsInj(opt, sols);
        }
        else if (opt.type == WellType::productor) {
            SetupOptsProd(opt);
        }
        else {
            OCP_ABORT("Inavailable Well Type!");
        }
    }
}


void Well::SetupOptsInj(WellOpt& opt, const vector<SolventINJ>& sols)
{
    const auto g = mixture->GasIndex();
    const auto w = mixture->WatIndex();
 
    if (opt.injFluidName == "WAT") {     
        if (w < 0) OCP_ABORT("WRONG INJECTED FLUID -- NO WATER!");
        opt.injZi.assign(nc, 0);
        opt.injPhase     = PhaseType::wat;
        opt.injZi.back() =  1.0;
    }
    else if (!sols.empty())
    {
        opt.injPhase  = PhaseType::gas;
        for (auto& s : sols) {
            if (opt.injFluidName == s.name) {
                opt.injZi = s.data;
                opt.injZi.resize(nc);
            }
        }
        if (opt.injZi.size() != nc) {
            OCP_ABORT("Inavailable INJECTED FLUID -- " + opt.injFluidName + "!");
        }
    }
    else if (opt.injFluidName == "GAS") {
        if (g < 0) OCP_ABORT("WRONG INJECTED FLUID -- NO GAS!");
        opt.injZi.assign(nc, 0);
        opt.injPhase  = PhaseType::gas;
        opt.injZi[g]  =  1.0;
    }
    else {
        OCP_ABORT("Inavailable INJECTED FLUID!");
    }

    if (mixture->MixtureType() != OCPMixtureType::THERMALK_OW)
        opt.injTemp = rsTemp;
}


void Well::SetupOptsProd(WellOpt& opt)
{
    const auto o = mixture->OilIndex();
    const auto g = mixture->GasIndex();
    const auto w = mixture->WatIndex();
    const auto l = mixture->LiquidIndex();

    opt.prodPhaseWeight.assign(np, 0);
    switch (opt.mode) {
    case WellOptMode::bhp:
        break;
    case WellOptMode::orate:
        if (o < 0) OCP_ABORT("WRONG PRODUCTED FLUID -- NO OIL!");
        opt.prodPhaseWeight[o] = 1.0;
        break;
    case WellOptMode::grate:
        if (g < 0) OCP_ABORT("WRONG PRODUCTED FLUID -- NO GAS!");
        opt.prodPhaseWeight[g] = 1.0;
        break;
    case WellOptMode::wrate:
        if (w < 0) OCP_ABORT("WRONG PRODUCTED FLUID -- NO WATER!");
        opt.prodPhaseWeight[w] = 1.0;
        break;
    case WellOptMode::lrate:
        if (l.size() == 0) OCP_ABORT("WRONG PRODUCTED FLUID -- NO LIQUID!");
        for (auto& i : l)  opt.prodPhaseWeight[i] = 1.0;
        break;
    default:
        OCP_ABORT("WRONG Opt Mode!");
        break;
    }
}


void Well::ShowPerfStatus(const Bulk& bk) const
{
    OCP_FUNCNAME;

    const BulkVarSet& bvs = bk.vs;

    cout << fixed;
    cout << "----------------------------" << endl;
    cout << name << ":    " << (USI)opt.mode << "   " << setprecision(3) << bhp << endl;
    for (USI p = 0; p < numPerf; p++) {
        vector<OCP_DBL> Qitmp(perf[p].qi_lbmol);
        // OCP_DBL         qt = Dnorm1(nc, &Qitmp[0]);
        OCP_USI n = perf[p].location;
        cout << setw(3) << p << "   " << (USI)perf[p].state << "   " << setw(6)
             << perf[p].location << "  "
             << setprecision(6) << perf[p].WI << "  "               // ccf
             << setprecision(3) << perf[p].radius << "  "           // ccf
             << setw(8) << setprecision(4) << perf[p].kh << "  "    // kh
             << setw(8) << setprecision(2) << perf[p].depth << "  " // depth
             << setprecision(3) << perf[p].P << "  "                // Pp
             << setw(10) << setprecision(3) << bvs.P[n] << "   " // Pb
             << setw(8) << perf[p].qi_lbmol[nc - 1] << "   " << setw(6)
             << setprecision(6) << bvs.S[n * np + 0] << "   " << setw(6)
             << setprecision(6) << bvs.S[n * np + 1] << "   " << setw(6)
             << setprecision(6) << bvs.S[n * np + 2] << endl;
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