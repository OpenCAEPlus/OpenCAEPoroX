/*! \file    AllWells.cpp
 *  \brief   AllWells class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "AllWells.hpp"

/////////////////////////////////////////////////////////////////////
// General
/////////////////////////////////////////////////////////////////////

void AllWells::InputParam(const ParamWell& paramWell, const Domain& domain)
{

    OCP_FUNCNAME;

    for (auto& s : paramWell.solSet) {
        solvents.push_back(s);
    }

    Psurf = paramWell.Psurf;
    Tsurf = paramWell.Tsurf;
 
    const auto my_well = domain.GetWell();
    numWell            = my_well.size();

    if (paramWell.thermal) {
        for (USI w = 0; w < numWell; w++) {
            wells.push_back(new PeacemanWellT());
        }
    }
    else {
        for (USI w = 0; w < numWell; w++) {
            wells.push_back(new PeacemanWellIsoT());
        }
    }

    
    USI         t = paramWell.criticalTime.size();
    vector<USI>         wellOptTime;
    vector<WellOptPair> tmpOptParam;
    for (USI wdst = 0; wdst < numWell; wdst++) {
        const OCP_USI wsrc = my_well[wdst];
        wells[wdst]->name  = paramWell.well[wsrc].name;
        wells[wdst]->group = paramWell.well[wsrc].group;
        wells[wdst]->depth = paramWell.well[wsrc].depth;
        wells[wdst]->InputPerfo(paramWell.well[wsrc], domain, wsrc);
        wells[wdst]->Psurf         = Psurf;
        wells[wdst]->Tsurf         = Tsurf;
        wells[wdst]->ifUseUnweight = paramWell.well[wsrc].ifUseUnweight;
        // opt
        wells[wdst]->optSet.resize(t);

        // If identical optParam.d occurs, use the last one
        tmpOptParam.clear();
        const USI n0 = paramWell.well[wsrc].optParam.size();
        if (n0 == 0) {
            OCP_ABORT("NO Well Opt Param for -- " + paramWell.well[wsrc].name);
        }

        for (USI i = 0; i < n0 - 1; i++) {
            if (paramWell.well[wsrc].optParam[i].d !=
                paramWell.well[wsrc].optParam[i + 1].d) {
                tmpOptParam.push_back(paramWell.well[wsrc].optParam[i]);
            }
        }
        tmpOptParam.push_back(paramWell.well[wsrc].optParam.back());

        const USI n = tmpOptParam.size();
        wellOptTime.clear();
        wellOptTime.resize(n + 1);
        for (USI i = 0; i < n; i++) {
            wellOptTime[i] = tmpOptParam[i].d;
        }
        wellOptTime.back() = t;
        for (USI i = 0; i < n; i++) {
            for (USI d = wellOptTime[i]; d < wellOptTime[i + 1]; d++) {
                wells[wdst]->optSet[d] = WellOpt(tmpOptParam[i].opt);
            }
        }
    }
}


void AllWells::Setup(const Bulk& bk)
{
    OCP_FUNCNAME;
    for (USI w = 0; w < numWell; w++) {
        wells[w]->Setup(bk, solvents);
    }
}


void AllWells::SetupWellGroup(const Bulk& bk)
{
    wellGroup.clear();
    // Field Group, contain all wells
    wellGroup.push_back(WellGroup("Field"));
    for (USI w = 0; w < numWell; w++) {
        wellGroup[0].wId.push_back(w);
        if (wells[w]->WellType() == WellType::injector)
            wellGroup[0].wIdINJ.push_back(w);
        else
            wellGroup[0].wIdPROD.push_back(w);
    }

    // other subgroups
    USI glen = 1;
    USI g    = 1;
    for (USI w = 0; w < numWell; w++) {
        for (g = 1; g < glen; g++) {
            if (wells[w]->group == wellGroup[g].name) {
                // existing group
                wellGroup[g].wId.push_back(w);
                if (wells[w]->WellType() == WellType::injector)
                    wellGroup[g].wIdINJ.push_back(w);
                else
                    wellGroup[g].wIdPROD.push_back(w);
                break;
            }
        }
        if (g == glen && wells[w]->group != "Field") {
            // new group
            wellGroup.push_back(WellGroup(wells[w]->group));
            wellGroup[glen].wId.push_back(w);
            if (wells[w]->WellType() == WellType::injector)
                wellGroup[glen].wIdINJ.push_back(w);
            else
                wellGroup[glen].wIdPROD.push_back(w);
            glen++;
        }
    }

    numGroup = wellGroup.size();
}


void AllWells::SetupWellBulk(Bulk& bk) const
{
    for (auto& w : wells) {
        if (w->IsOpen()) {
            for (auto& p : w->perf) {
                bk.AddWellBulkId(p.Location());
            }
        }
    }
}


void AllWells::ApplyControl(const USI& i)
{
    OCP_FUNCNAME;
    wellOptChange = OCP_FALSE;
    for (USI w = 0; w < numWell; w++) {
        if (wells[w]->ApplyOpt(i)) wellOptChange = OCP_TRUE;
    }
}

void AllWells::InitBHP(const Bulk& bk)
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        wells[w]->InitWellP(bk);
    }
}

void AllWells::PrepareWell(const Bulk& bk)
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        wells[w]->CheckOptMode(bk);
        wells[w]->CalFluxInit(bk);
    }
}


void AllWells::CalFlux(const Bulk& bk)
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        wells[w]->CalFlux(bk);
    }
}


void AllWells::CalIPRT(const Bulk& bk, OCP_DBL dt)
{
    OCP_FUNCNAME;

    FGIR = 0;
    FWIR = 0;
    FOPR = 0;
    FGPR = 0;
    FWPR = 0;
    for (USI w = 0; w < numWell; w++) {
     
        wells[w]->CalIPRate(bk, dt);
        
        FGIR += wells[w]->WGIR;
        FWIR += wells[w]->WWIR;
        FOPR += wells[w]->WOPR;
        FGPR += wells[w]->WGPR;
        FWPR += wells[w]->WWPR;       
    }
    FGIT += FGIR * dt;
    FWIT += FWIR * dt;
    FOPT += FOPR * dt;
    FGPt += FGPR * dt;
    FWPT += FWPR * dt;
}


ReservoirState AllWells::CheckP(const Bulk& bk)
{
    OCP_FUNCNAME;

    OCP_BOOL flagSwitch = OCP_FALSE;
    OCP_BOOL flagCrossf = OCP_FALSE;

    for (USI w = 0; w < numWell; w++) {

        ReservoirState flag = wells[w]->CheckP(bk);

        switch (flag) 
        {
        case ReservoirState::well_negative_P:
            return ReservoirState::well_negative_P;

        case ReservoirState::well_switch_BHPm:
            flagSwitch = OCP_TRUE;
            break;

        case ReservoirState::well_crossflow:
            flagCrossf = OCP_TRUE;
            break;

        case ReservoirState::well_success:
        default:
            break;
        }
    }

    if (flagSwitch || flagCrossf) return ReservoirState::well_switch_BHPm;

    return ReservoirState::well_success;
}

// return the index of Specified well name
USI AllWells::GetIndex(const string& name) const
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w]->name == name) {
            return w;    
        }
    }
    return -1; 
}

USI AllWells::GetWellPerfNum() const
{
    USI numPerf = 0;
    for (USI w = 0; w < numWell; w++) {
        numPerf += wells[w]->numPerf;
    }
    return numPerf;
}

USI AllWells::GetMaxWellPerNum() const
{
    OCP_FUNCNAME;

    USI m = 0;
    for (USI w = 0; w < numWell; w++) {
        m = max(m, wells[w]->numPerf);
    }
    return m;
}


USI AllWells::GetNumOpenWell() const
{
    USI nw = 0;
    for (const auto& w : wells) {
        if (w->IsOpen()) {
            nw++;
        }
    }
    return nw;
}


void AllWells::SetWellVal() const
{
    if (!useVTK) return;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w]->opt.state == WellState::open) {          
            if (wells[w]->opt.mode == WellOptMode::bhp)
                wellVal[w] = wells[w]->opt.tarBHP;
            else
                wellVal[w] = wells[w]->opt.tarRate;
        } else {
            wellVal[w] = 0;
        }
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