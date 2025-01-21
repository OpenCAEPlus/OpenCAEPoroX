/*! \file    OCPOutput.cpp
 *  \brief   OCPOutput class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPOutput.hpp"


OCP_INT  file_myrank;
MPI_Comm file_myComm;


void ItersInfo::Update(const OCPNRsuite& NRs)
{
    numTstep++;
    NRt  += NRs.GetIterNR();
    NRwt += NRs.GetIterNRw();
    LSt  += NRs.GetIterLS();
    LSwt += NRs.GetIterLSw();
}


void Summary::SetupOutputTerm(const OutputSummary& summary_param)
{
    FPR  = summary_param.FPR;
    FTR  = summary_param.FTR;
    FOPR = summary_param.FOPR;
    FOPT = summary_param.FOPT;
    FGPR = summary_param.FGPR;
    FGPt = summary_param.FGPt;
    FWPR = summary_param.FWPR;
    FWPT = summary_param.FWPT;
    FGIR = summary_param.FGIR;
    FGIT = summary_param.FGIT;
    FWIR = summary_param.FWIR;
    FWIT = summary_param.FWIT;

    WOPR = summary_param.WOPR;
    WOPT = summary_param.WOPT;
    WGPR = summary_param.WGPR;
    WGPT = summary_param.WGPT;
    WWPR = summary_param.WWPR;
    WWPT = summary_param.WWPT;
    WGIR = summary_param.WGIR;
    WGIT = summary_param.WGIT;
    WWIR = summary_param.WWIR;
    WWIT = summary_param.WWIT;

    WBHP = summary_param.WBHP;
    DG   = summary_param.DG;

    BPR  = summary_param.BPR;
    SOIL = summary_param.SOIL;
    SGAS = summary_param.SGAS;
    SWAT = summary_param.SWAT;

}

void Summary::Setup(const OutputSummary& summary_param, const Reservoir& rs)
{

    SetupOutputTerm(summary_param);

    const USI maxRowNum = 1000;

    Sumdata.push_back(SumItem("TIME", "-", "DAY", "fixed", maxRowNum));
    Sumdata.push_back(SumItem("TimeStep", "-", "DAY", "fixed", maxRowNum));
    Sumdata.push_back(SumItem("NRiter", "-", "-", "int", maxRowNum));
    Sumdata.push_back(SumItem("NRiterW", "-", "-", "int", maxRowNum));
    Sumdata.push_back(SumItem("NRiter(DDM)", "-", "-", "int", maxRowNum));
    Sumdata.push_back(SumItem("NRiterW(DDM)", "-", "-", "int", maxRowNum));
    Sumdata.push_back(SumItem("LSiter", "-", "-", "int", maxRowNum));
    Sumdata.push_back(SumItem("LS/NR", "-", "-", "float", maxRowNum));
    Sumdata.push_back(SumItem("Runtime", "-", "s", "float", maxRowNum));
    if (FPR) Sumdata.push_back(SumItem("FPR", "-", "PSIA", "float", maxRowNum));
    if (FPR) Sumdata.push_back(SumItem("Volume", "Hydrocarbon", "Ft3", "float", maxRowNum));
    if (FTR) Sumdata.push_back(SumItem("FTR", "-", "F", "float", maxRowNum));
    if (FTR) Sumdata.push_back(SumItem("Volume", "Grid", "Ft3", "float", maxRowNum));
    if (FOPR) Sumdata.push_back(SumItem("FOPR", "-", "STB/DAY", "float", maxRowNum));
    if (FOPT) Sumdata.push_back(SumItem("FOPT", "-", "STB", "float", maxRowNum));
    if (FGPR) Sumdata.push_back(SumItem("FGPR", "-", "MSCF/DAY", "float", maxRowNum));
    if (FGPt) Sumdata.push_back(SumItem("FGPT", "-", "MSCF", "float", maxRowNum));
    if (FWPR) Sumdata.push_back(SumItem("FWPR", "-", "STB/DAY", "float", maxRowNum));
    if (FWPT) Sumdata.push_back(SumItem("FWPT", "-", "STB", "float", maxRowNum));
    if (FGIR) Sumdata.push_back(SumItem("FGIR", "-", "MSCF/DAY", "float", maxRowNum));
    if (FGIT) Sumdata.push_back(SumItem("FGIT", "-", "MSCF", "float", maxRowNum));
    if (FWIR) Sumdata.push_back(SumItem("FWIR", "-", "STB/DAY", "float", maxRowNum));
    if (FWIT) Sumdata.push_back(SumItem("FWIT", "-", "STB", "float", maxRowNum));

    // const Grid&     initGrid = rs.grid;
    // const USI sp = initGrid.GetNumDigitIJK();
    const AllWells& wells    = rs.allWells;
    const USI wellnum = wells.GetWellNum();
    string    wellname;
    USI       num;

    if (WOPR.activity) {
        if (WOPR.obj[0] == "All") {
            WOPR.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WOPR.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WOPR", wellname, "STB/DAY", "float", maxRowNum));
                WOPR.index.push_back(w);
            }
        } else {
            num = WOPR.obj.size();
            for (USI w = 0; w < num; w++) {
                const OCP_INT tmpW = wells.GetIndex(WOPR.obj[w]);
                if (tmpW < 0) continue;
                Sumdata.push_back(SumItem("WOPR", WOPR.obj[w], "STB/DAY", "float", maxRowNum));
                WOPR.index.push_back(tmpW);
            }
        }
    }

    if (WOPT.activity) {
        if (WOPT.obj[0] == "All") {
            WOPT.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WOPT.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WOPT", wellname, "STB", "float", maxRowNum));
                WOPT.index.push_back(w);
            }
        } else {
            num = WOPT.obj.size();
            for (USI w = 0; w < num; w++) {
                const OCP_INT tmpW = wells.GetIndex(WOPT.obj[w]);
                if (tmpW < 0) continue;
                Sumdata.push_back(SumItem("WOPT", WOPT.obj[w], "STB", "float", maxRowNum));
                WOPT.index.push_back(tmpW);
            }
        }
    }

    if (WGPR.activity) {
        if (WGPR.obj[0] == "All") {
            WGPR.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WGPR.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WGPR", wellname, "MSCF/DAY", "float", maxRowNum));
                WGPR.index.push_back(w);
            }
        } else {
            num = WGPR.obj.size();
            for (USI w = 0; w < num; w++) {
                const OCP_INT tmpW = wells.GetIndex(WGPR.obj[w]);
                if (tmpW < 0) continue;
                Sumdata.push_back(SumItem("WGPR", WGPR.obj[w], "MSCF/DAY", "float", maxRowNum));
                WGPR.index.push_back(tmpW);
            }
        }
    }

    if (WGPT.activity) {
        if (WGPT.obj[0] == "All") {
            WGPT.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WGPT.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WGPT", wellname, "MSCF", "float", maxRowNum));
                WGPT.index.push_back(w);
            }
        } else {
            num = WGPT.obj.size();
            for (USI w = 0; w < num; w++) {
                const OCP_INT tmpW = wells.GetIndex(WGPT.obj[w]);
                if (tmpW < 0) continue;
                Sumdata.push_back(SumItem("WGPT", WGPT.obj[w], "MSCF", "float", maxRowNum));
                WGPT.index.push_back(tmpW);
            }
        }
    }

    if (WWPR.activity) {
        if (WWPR.obj[0] == "All") {
            WWPR.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WWPR.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WWPR", wellname, "STB/DAY", "float", maxRowNum));
                WWPR.index.push_back(w);
            }
        } else {
            num = WWPR.obj.size();
            for (USI w = 0; w < num; w++) {
                const OCP_INT tmpW = wells.GetIndex(WWPR.obj[w]);
                if (tmpW < 0) continue;
                Sumdata.push_back(SumItem("WWPR", WWPR.obj[w], "STB/DAY", "float", maxRowNum));
                WWPR.index.push_back(tmpW);
            }
        }
    }

    if (WWPT.activity) {
        if (WWPT.obj[0] == "All") {
            WWPT.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WWPT.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WWPT", wellname, "STB", "float", maxRowNum));
                WWPT.index.push_back(w);
            }
        } else {
            num = WWPT.obj.size();
            for (USI w = 0; w < num; w++) {
                const OCP_INT tmpW = wells.GetIndex(WWPT.obj[w]);
                if (tmpW < 0) continue;
                Sumdata.push_back(SumItem("WWPT", WWPT.obj[w], "STB", "float", maxRowNum));
                WWPT.index.push_back(tmpW);
            }
        }
    }

    if (WGIR.activity) {
        if (WGIR.obj[0] == "All") {
            WGIR.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WGIR.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WGIR", wellname, "MSCF/DAY", "float", maxRowNum));
                WGIR.index.push_back(w);
            }
        } else {
            num = WGIR.obj.size();
            for (USI w = 0; w < num; w++) {
                const OCP_INT tmpW = wells.GetIndex(WGIR.obj[w]);
                if (tmpW < 0) continue;
                Sumdata.push_back(SumItem("WGIR", WGIR.obj[w], "MSCF/DAY", "float", maxRowNum));
                WGIR.index.push_back(tmpW);
            }
        }
    }

    if (WGIT.activity) {
        if (WGIT.obj[0] == "All") {
            WGIT.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WGIT.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WGIT", wellname, "MSCF", "float", maxRowNum));
                WGIT.index.push_back(w);
            }
        } else {
            num = WGIT.obj.size();
            for (USI w = 0; w < num; w++) {
                const OCP_INT tmpW = wells.GetIndex(WGIT.obj[w]);
                if (tmpW < 0) continue;
                Sumdata.push_back(SumItem("WGIT", WGIT.obj[w], "MSCF", "float", maxRowNum));
                WGIT.index.push_back(tmpW);
            }
        }
    }

    if (WWIR.activity) {
        if (WWIR.obj[0] == "All") {
            WWIR.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WWIR.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WWIR", wellname, "STB/DAY", "float", maxRowNum));
                WWIR.index.push_back(w);
            }
        } else {
            num = WWIR.obj.size();
            for (USI w = 0; w < num; w++) {
                const OCP_INT tmpW = wells.GetIndex(WWIR.obj[w]);
                if (tmpW < 0) continue;
                Sumdata.push_back(SumItem("WWIR", WWIR.obj[w], "STB/DAY", "float", maxRowNum));
                WWIR.index.push_back(tmpW);
            }
        }
    }

    if (WWIT.activity) {
        if (WWIT.obj[0] == "All") {
            WWIT.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WWIT.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WWIT", wellname, "STB", "float", maxRowNum));
                WWIT.index.push_back(w);
            }
        } else {
            num = WWIT.obj.size();
            for (USI w = 0; w < num; w++) {
                const OCP_INT tmpW = wells.GetIndex(WWIT.obj[w]);
                if (tmpW < 0) continue;
                Sumdata.push_back(SumItem("WWIT", WWIT.obj[w], "STB", "float", maxRowNum));
                WWIT.index.push_back(tmpW);
            }
        }
    }

    if (WBHP.activity) {
        if (WBHP.obj[0] == "All") {
            WBHP.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WBHP.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WBHP", wellname, "PSIA", "float", maxRowNum));
                WBHP.index.push_back(w);
            }
        } else {
            num = WBHP.obj.size();
            for (USI w = 0; w < num; w++) {
                const OCP_INT tmpW = wells.GetIndex(WBHP.obj[w]);
                if (tmpW < 0) continue;
                Sumdata.push_back(SumItem("WBHP", WBHP.obj[w], "PSIA", "float", maxRowNum));
                WBHP.index.push_back(tmpW);
            }
        }
    }

    //if (DG.activity) {
    //    if (DG.obj[0] == "All") {
    //        DG.obj.clear();
    //        for (USI w = 0; w < wellnum; w++) {
    //            string wellname = wells.GetWellName(w);
    //            DG.obj.push_back(wellname);
    //            USI perfnum = wells.GetWellPerfNum(w);
    //            for (USI p = 0; p < perfnum; p++) {
    //                Sumdata.push_back(SumItem("DG", wellname + " Perf" + to_string(p),
    //                                          "PSIA", "float", maxRowNum));
    //                DG.index.push_back(w);
    //            }
    //        }
    //    } else {
    //        num = DG.obj.size();
    //        for (USI w = 0; w < num; w++) {
    //            USI wId     = wells.GetIndex(DG.obj[w]);
    //            USI perfnum = wells.GetWellPerfNum(wId);
    //            for (USI p = 0; p < perfnum; p++) {
    //                Sumdata.push_back(SumItem("DG", DG.obj[w] + " P" + to_string(p),
    //                                          "PSIA", "float", maxRowNum));
    //                DG.index.push_back(wId);
    //            }
    //        }
    //    }
    //}

    //if (BPR.activity) {
    //    num = BPR.obj.size();
    //    for (USI i = 0; i < num; i++) {
    //        string temp = GetIJKformat(to_string(BPR.obj[i].I), to_string(BPR.obj[i].J),
    //                                   to_string(BPR.obj[i].K), sp);
    //        Sumdata.push_back(SumItem("BPR", temp, "PSIA", "float"));
    //        USI I = BPR.obj[i].I - 1;
    //        USI J = BPR.obj[i].J - 1;
    //        USI K = BPR.obj[i].K - 1;

    //        const OCP_INT tarId = initGrid.GetActIndex(I, J, K);
    //        if (tarId < 0)
    //            OCP_WARNING("Non Fluid Grid: " + GetIJKformat(I + 1, J + 1, K + 1, sp) +
    //                        "in BPR in SUMMARY");
    //        else
    //            BPR.index.push_back(tarId);
    //    }
    //}

    //if (SOIL.activity) {
    //    num = SOIL.obj.size();
    //    for (USI i = 0; i < num; i++) {
    //        string temp =
    //            GetIJKformat(to_string(SOIL.obj[i].I), to_string(SOIL.obj[i].J),
    //                         to_string(SOIL.obj[i].K), sp);
    //        Sumdata.push_back(SumItem("SOIL", temp, "   ", "float"));
    //        USI I = SOIL.obj[i].I - 1;
    //        USI J = SOIL.obj[i].J - 1;
    //        USI K = SOIL.obj[i].K - 1;

    //        const OCP_INT tarId = initGrid.GetActIndex(I, J, K);
    //        if (tarId < 0)
    //            OCP_WARNING("Non Fluid Grid: " + GetIJKformat(I + 1, J + 1, K + 1, sp) +
    //                        "in SOIL in SUMMARY");
    //        else
    //            SOIL.index.push_back(tarId);
    //    }
    //}

    //if (SGAS.activity) {
    //    num = SGAS.obj.size();
    //    for (USI i = 0; i < num; i++) {
    //        string temp =
    //            GetIJKformat(to_string(SGAS.obj[i].I), to_string(SGAS.obj[i].J),
    //                         to_string(SGAS.obj[i].K), sp);
    //        Sumdata.push_back(SumItem("SGAS", temp, "   ", "float"));
    //        USI I = SGAS.obj[i].I - 1;
    //        USI J = SGAS.obj[i].J - 1;
    //        USI K = SGAS.obj[i].K - 1;

    //        const OCP_INT tarId = initGrid.GetActIndex(I, J, K);
    //        if (tarId < 0)
    //            OCP_WARNING("Non Fluid Grid: " + GetIJKformat(I + 1, J + 1, K + 1, sp) +
    //                        "in SOIL in SUMMARY");
    //        else
    //            SGAS.index.push_back(tarId);
    //    }
    //}

    //if (SWAT.activity) {
    //    num = SWAT.obj.size();
    //    for (USI i = 0; i < num; i++) {
    //        string temp =
    //            GetIJKformat(to_string(SWAT.obj[i].I), to_string(SWAT.obj[i].J),
    //                         to_string(SWAT.obj[i].K), sp);
    //        Sumdata.push_back(SumItem("SWAT", temp, "   ", "float"));
    //        USI I = SWAT.obj[i].I - 1;
    //        USI J = SWAT.obj[i].J - 1;
    //        USI K = SWAT.obj[i].K - 1;

    //        const OCP_INT tarId = initGrid.GetActIndex(I, J, K);
    //        if (tarId < 0)
    //            OCP_WARNING("Non Fluid Grid: " + GetIJKformat(I + 1, J + 1, K + 1, sp) +
    //                        "in SWAT in SUMMARY");
    //        else
    //            SWAT.index.push_back(tarId);
    //    }
    //}

}

void Summary::SetVal(const Reservoir& rs, const OCPControl& ctrl, const ItersInfo& iters, GetWallTime& timer)
{
    const Bulk&     bulk  = rs.bulk;
    const AllWells& wells = rs.allWells;

    USI n = 0;

    // TIME
    Sumdata[n++].val.push_back(ctrl.time.GetCurrentTime());
    // TimeStep
    Sumdata[n++].val.push_back(ctrl.time.GetLastDt());
    // NRiter
    Sumdata[n++].val.push_back(iters.GetNRt());
    // NRiter(waste)
    Sumdata[n++].val.push_back(iters.GetNRwt());
    // NRiter(DDM)
    Sumdata[n++].val.push_back(static_cast<OCP_DBL>(OCPITER_NR_DDM));
    // NRiter(DDM)(waste)
    Sumdata[n++].val.push_back(static_cast<OCP_DBL>(OCPITER_NRW_DDM));
    // LSiter
    Sumdata[n++].val.push_back(iters.GetLSt());
    // LSiter / NRiter
    Sumdata[n++].val.push_back(Sumdata[2].val.back() == 0 ? 0 : 1.0 * Sumdata[6].val.back() / Sumdata[2].val.back());
    // runtime
    Sumdata[n++].val.push_back(timer.Stop());

    OCP_DBL  tmpV = 0;
    if (FPR) Sumdata[n++].val.push_back(bulk.CalFPR(tmpV));
    if (FPR) Sumdata[n++].val.push_back(tmpV);
    if (FTR) Sumdata[n++].val.push_back(bulk.CalFTR(tmpV));
    if (FTR) Sumdata[n++].val.push_back(tmpV);
    if (FOPR) Sumdata[n++].val.push_back(wells.GetFOPR());
    if (FOPT) Sumdata[n++].val.push_back(wells.GetFOPT());
    if (FGPR) Sumdata[n++].val.push_back(wells.GetFGPR());
    if (FGPt) Sumdata[n++].val.push_back(wells.GetFGPT());
    if (FWPR) Sumdata[n++].val.push_back(wells.GetFWPR());
    if (FWPT) Sumdata[n++].val.push_back(wells.GetFWPT());
    if (FGIR) Sumdata[n++].val.push_back(wells.GetFGIR());
    if (FGIT) Sumdata[n++].val.push_back(wells.GetFGIT());
    if (FWIR) Sumdata[n++].val.push_back(wells.GetFWIR());
    if (FWIT) Sumdata[n++].val.push_back(wells.GetFWIT());

    USI len = 0;
    // WOPR
    len = WOPR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWOPR(WOPR.index[w]));

    // WOPT
    len = WOPT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWOPT(WOPT.index[w]));

    // WGPR
    len = WGPR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWGPR(WGPR.index[w]));

    // WGPT
    len = WGPT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWGPT(WGPT.index[w]));

    // WWPR
    len = WWPR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWWPR(WWPR.index[w]));

    // WWPT
    len = WWPT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWWPT(WWPT.index[w]));

    // WGIR
    len = WGIR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWGIR(WGIR.index[w]));

    // WGIT
    len = WGIT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWGIT(WGIT.index[w]));

    // WWIR
    len = WWIR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWWIR(WWIR.index[w]));

    // WWIT
    len = WWIT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWWIT(WWIT.index[w]));

    // WBHP
    len = WBHP.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWBHP(WBHP.index[w]));

    //// DG
    //len = DG.obj.size();
    //for (USI w = 0; w < len; w++) {
    //    USI numperf = rs.allWells.GetWellPerfNum(DG.index[w]);
    //    for (USI p = 0; p < numperf; p++) {
    //        Sumdata[n++].val.push_back(wells.GetWellDG(DG.index[w], p));
    //    }
    //}

    //// BPR
    //len = BPR.index.size();
    //for (USI i = 0; i < len; i++) Sumdata[n++].val.push_back(bulk.GetP(BPR.index[i]));

    //// SOIL
    //len = SOIL.index.size();
    //for (USI i = 0; i < len; i++)
    //    Sumdata[n++].val.push_back(bulk.GetSOIL(SOIL.index[i]));

    //// SGAS
    //len = SGAS.index.size();
    //for (USI i = 0; i < len; i++)
    //    Sumdata[n++].val.push_back(bulk.GetSGAS(SGAS.index[i]));

    //// SWAT
    //len = SWAT.index.size();
    //for (USI i = 0; i < len; i++)
    //    Sumdata[n++].val.push_back(bulk.GetSWAT(SWAT.index[i]));
}

/// Write output information in the dir/SUMMARY.out file.
void Summary::PrintInfo(const string& dir, const string& filename, const OCP_INT& rank) const
{
    string   FileOut;
    if (rank >= 0) {
        FileOut = dir + "proc_" + to_string(rank) + "_SUMMARY.out";
    }
    else {
        FileOut = dir + "SUMMARY.out";
    }

    ofstream outF(FileOut);
    if (!outF.is_open()) {
        OCP_ABORT("Can not open " + FileOut);
    }

    const USI ns  = 12;
    const USI col = 10;
    const USI num = Sumdata.size();
    const USI len = Sumdata[0].val.size();

    USI row = 0;
    USI id  = 0;
    USI ID  = 1;

    outF << "SUMMARY OF RUN " + dir + filename << " -- " << Sumdata[0].val.size() << " time step\n";

    while (id != num) {

        outF << "Row " << ++row << "\n";

        // Item
        outF << "\t" << setw(ns) << Sumdata[0].Item;

        id = ID;
        for (USI i = 1; i < col; i++) {
            outF << "\t" << setw(ns) << Sumdata[id++].Item;
            if (id == num) break;
        }
        outF << "\n";

        // Unit
        outF << "\t" << setw(ns) << Sumdata[0].Unit;

        id = ID;
        for (USI i = 1; i < col; i++) {
            outF << "\t" << setw(ns) << Sumdata[id++].Unit;
            if (id == num) break;
        }
        outF << "\n";

        // Obj Name
        outF << "\t" << setw(ns) << Sumdata[0].Obj;

        id = ID;
        for (USI i = 1; i < col; i++) {
            outF << "\t" << setw(ns) << Sumdata[id++].Obj;
            if (id == num) break;
        }
        outF << "\n";

        // Data
        for (USI l = 0; l < len; l++) {

            // Time
            outF << "\t" << setw(ns) << fixed << setprecision(3) << Sumdata[0].val[l];

            id = ID;
            for (USI i = 1; i < col; i++) {
                if (Sumdata[id].Type == "int") {
                    outF << fixed << setprecision(0);
                } else if (Sumdata[id].Type == "fixed") {
                    outF << fixed << setprecision(3);
                } else if (Sumdata[id].Type == "float") {
                    outF << scientific << setprecision(5);
                }
                outF << "\t" << setw(ns) << Sumdata[id++].val[l];
                if (id == num) break;
            }
            outF << "\n";
        }
        ID += (col - 1);
        outF << "\n";
    }
    outF.close();
}


// By files
//void Summary::PostProcess(const string& dir, const string& filename, const OCP_INT& numproc) const
//{
//    vector<SumItem>* sumdata = const_cast<vector<SumItem>*>(&Sumdata);
//    sumdata->clear();
//
//    OCP_INT         bId = 0;
//    USI             varlen = 0;
//    OCP_USI         rowNum;
//
//    vector<SumItem> mySum;
//    vector<string>  buffer;
//    OCP_BOOL        flag = OCP_FALSE;
//
//    for (USI p = 0; p < numproc; p++) {
//        const string myFile = dir + "proc_" + to_string(p) + "_SUMMARY.out";
//        ifstream ifs(myFile, ios::in);
//
//        if (!ifs) {
//            OCP_MESSAGE("Trying to open file: " << (myFile));
//            OCP_ABORT("Failed to open the input file!");
//        }
//
//        // Get time steps
//        ReadLine(ifs, buffer, OCP_FALSE);
//        for (USI i = 0; i < buffer.size(); i++) {
//            if (buffer[i] == "--") {
//                rowNum = stoi(buffer[i + 1]);
//                break;
//            }
//        }
//
//        while (!ifs.eof()) {
//            if (buffer[0] == "TIME") {
//                // begin to input title of data
//                flag = OCP_TRUE;
//                bId += varlen;
//                varlen = buffer.size();
//                for (const auto& v : buffer) {
//                    mySum.push_back(SumItem(v, rowNum));
//                }
//                ReadLine(ifs, buffer, OCP_FALSE);
//                OCP_ASSERT(varlen == buffer.size(), "Mismatch Value");
//                for (USI i = 0; i < varlen; i++) {
//                    mySum[bId + i].Unit = buffer[i];
//                }
//                ReadLine(ifs, buffer, OCP_FALSE);
//                OCP_ASSERT(varlen == buffer.size(), "Mismatch Value");
//                for (USI i = 0; i < varlen; i++) {
//                    mySum[bId + i].Obj = buffer[i];
//                }
//            }
//            ReadLine(ifs, buffer, OCP_FALSE);
//            if (flag && buffer.size() == varlen && buffer[0] != "Row") {
//                // begin to input value of data               
//                for (USI i = 0; i < varlen; i++) {
//                    mySum[bId + i].val.push_back(stod(buffer[i]));
//                }
//            }
//            else {
//                flag = OCP_FALSE;
//            }
//        }
//
//        ifs.close();
//        if (remove(myFile.c_str()) != 0) {
//            OCP_WARNING("Failed to delete " + myFile);
//        }
//    }
//
//    // post process
//    // combine some values
//    vector<SumItem> tmpdata;
//    for (const auto& s : mySum) {
//        USI n = 0;
//        for (; n < tmpdata.size(); n++) {
//            if (tmpdata[n] == s) {
//                if (s.Item == "FOPR" || s.Item == "FGPR" || s.Item == "FWPR" ||
//                    s.Item == "FOPT" || s.Item == "FGPT" || s.Item == "FWPT" ||
//                    s.Item == "FGIR" || s.Item == "FGIT" || s.Item == "FWIR" ||
//                    s.Item == "FWIT") {
//                    for (USI i = 0; i < rowNum; i++) {
//                        sumdata->at(n).val[i] += s.val[i];
//                    }
//                }
//                else if (s.Item == "FPR" || s.Item == "FTR") {
//                    for (USI i = 0; i < rowNum; i++) {
//                        sumdata->at(n).val[i] += s.val[i] * (&s + 1)->val[i];
//                        sumdata->at(n + 1).val[i] += (&s + 1)->val[i];
//                    }
//                }
//                break;
//            }
//        }
//        if (n == sumdata->size()) {
//            sumdata->push_back(s);
//            tmpdata.push_back(SumItem(s.Item, s.Obj));
//            // prepare for calculating for FPR or FPT
//            if (s.Item == "Volume") {
//                for (USI i = 0; i < rowNum; i++) {
//                    sumdata->at(n - 1).val[i] *= sumdata->at(n).val[i];
//                }
//            }
//        }
//    }
//    // set output precision
//    for (USI n = 0; n < sumdata->size(); n++) {
//        if (sumdata->at(n).Item == "FPR" || sumdata->at(n).Item == "FTR") {
//            for (USI i = 0; i < rowNum; i++) {
//                sumdata->at(n).val[i] /= sumdata->at(n + 1).val[i];
//            }
//            sumdata->at(n).Type = "float";
//            continue;
//        }
//        if (sumdata->at(n).Item == "TIME") {
//            sumdata->at(n).Type = "fixed";
//            continue;
//        }
//        if (sumdata->at(n).Item == "NRiter" || sumdata->at(n).Item == "LSiter") {
//            sumdata->at(n).Type = "int";
//            continue;
//        }
//        else {
//            sumdata->at(n).Type = "float";
//            continue;
//        }
//    }
//
//    PrintInfo(dir, filename, -1);
//}


// By communication
void Summary::PostProcess(const string& dir, const string& filename, const OCP_INT& numproc) const
{
    vector<SumItem>* sumdata = const_cast<vector<SumItem>*>(&Sumdata);
    OCP_USI         rowNum = Sumdata[0].val.size();
    double* buffer = new double[rowNum];

    // post process
    // combine some values
    for (USI n = 0; n < sumdata->size(); n++) {
        auto& s = sumdata->at(n);
        if (s.Item == "FOPR" || s.Item == "FGPR" || s.Item == "FWPR" ||
            s.Item == "FOPT" || s.Item == "FGPT" || s.Item == "FWPT" ||
            s.Item == "FGIR" || s.Item == "FGIT" || s.Item == "FWIR" ||
            s.Item == "FWIT") {
            MPI_Reduce(s.val.data(), buffer, rowNum, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, file_myComm);
            memcpy(s.val.data(), buffer, sizeof(double) * rowNum);

        }
        else if (s.Item == "FPR" || s.Item == "FTR") {
            for (USI i = 0; i < rowNum; i++) {
                s.val[i] = s.val[i] * (&s + 1)->val[i];
            }
            MPI_Reduce((&s + 1)->val.data(), buffer, rowNum, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, file_myComm);
            memcpy((&s + 1)->val.data(), buffer, sizeof(double) * rowNum);
            MPI_Reduce(s.val.data(), buffer, rowNum, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, file_myComm);
            for (USI i = 0; i < rowNum; i++) {
                sumdata->at(n).val[i] = buffer[i] / sumdata->at(n + 1).val[i];
            }
        }
    }
    // gather extra column
    int* extra_column_count = new int[numproc];
    int local_extra_column_count = 0;
    int ItemStringLenPresum[4];

    for (USI n = 0; n < sumdata->size(); n++) {
        auto& s = sumdata->at(n);
        if (s.Item == "FOPR" || s.Item == "FGPR" || s.Item == "FWPR" ||
            s.Item == "FOPT" || s.Item == "FGPT" || s.Item == "FWPT" ||
            s.Item == "FGIR" || s.Item == "FGIT" || s.Item == "FWIR" ||
            s.Item == "FWIT" || s.Item == "FPR" || s.Item == "FTR" ||
            s.Item == "Volume" || s.Item == "TIME" || s.Item == "TimeStep" ||
            s.Item == "NRiter" || s.Item == "NRiterW" || s.Item == "NRiter(DDM)" ||
            s.Item == "NRiterW(DDM)" || s.Item == "LSiter" || s.Item == "LS/NR" ||
            s.Item == "Runtime") {
            continue;
        }
        local_extra_column_count++;
    }
    MPI_Gather(&local_extra_column_count, 1, MPI_INT, extra_column_count, 1, MPI_INT, MASTER_PROCESS, file_myComm);

    if (file_myrank == MASTER_PROCESS) {
        char ItemStringBuffer[256];
        for (USI proc = 0; proc < numproc; proc++) {
            if (proc != MASTER_PROCESS) {
                for (USI i = 0; i < extra_column_count[proc]; i++) {
                    MPI_Recv(ItemStringLenPresum, 4, MPI_INT, proc, 0, file_myComm, MPI_STATUS_IGNORE);
                    MPI_Recv(ItemStringBuffer, ItemStringLenPresum[3], MPI_CHAR, proc, 1, file_myComm, MPI_STATUS_IGNORE);
                    MPI_Recv(buffer, rowNum, MPI_DOUBLE, proc, 2, file_myComm, MPI_STATUS_IGNORE);
                    string Item(ItemStringBuffer, ItemStringLenPresum[0]);
                    string Obj(ItemStringBuffer + ItemStringLenPresum[0], ItemStringLenPresum[1] - ItemStringLenPresum[0]);
                    string Unit(ItemStringBuffer + ItemStringLenPresum[1], ItemStringLenPresum[2] - ItemStringLenPresum[1]);
                    string Type(ItemStringBuffer + ItemStringLenPresum[2], ItemStringLenPresum[3] - ItemStringLenPresum[2]);
                    SumItem column(Item, Obj);
                    bool flag = false;
                    for (USI n = 0; n < sumdata->size(); n++) {
                        if (sumdata->at(n) == column) {
                            flag = true;
                        }
                    }
                    if (flag) {
                        continue;
                    }
                    sumdata->emplace_back(Item, Obj, Unit, Type, rowNum);
                    for (USI i = 0; i < rowNum; i++)
                        sumdata->rbegin()->val.push_back(buffer[i]);
                }
            }
        }
    }
    else {
        string ItemStringBuffer;
        for (USI n = 0; n < sumdata->size(); n++) {
            auto& s = sumdata->at(n);
            if (s.Item == "FOPR" || s.Item == "FGPR" || s.Item == "FWPR" ||
                s.Item == "FOPT" || s.Item == "FGPT" || s.Item == "FWPT" ||
                s.Item == "FGIR" || s.Item == "FGIT" || s.Item == "FWIR" ||
                s.Item == "FWIT" || s.Item == "FPR" || s.Item == "FTR" ||
                s.Item == "Volume" || s.Item == "TIME" || s.Item == "TimeStep" ||
                s.Item == "NRiter" || s.Item == "NRiterW" || s.Item == "NRiter(DDM)" ||
                s.Item == "NRiterW(DDM)" || s.Item == "LSiter" || s.Item == "LS/NR" ||
                s.Item == "Runtime") {
                continue;
            }
            ItemStringLenPresum[0] = s.Item.length();
            ItemStringLenPresum[1] = ItemStringLenPresum[0] + s.Obj.length();
            ItemStringLenPresum[2] = ItemStringLenPresum[1] + s.Unit.length();
            ItemStringLenPresum[3] = ItemStringLenPresum[2] + s.Type.length();
            ItemStringBuffer = s.Item + s.Obj + s.Unit + s.Type;
            MPI_Send(ItemStringLenPresum, 4, MPI_INT, MASTER_PROCESS, 0, file_myComm);
            MPI_Send(ItemStringBuffer.c_str(), ItemStringLenPresum[3], MPI_CHAR, MASTER_PROCESS, 1, file_myComm);
            MPI_Send(s.val.data(), rowNum, MPI_DOUBLE, MASTER_PROCESS, 2, file_myComm);
        }
    }

    delete[]extra_column_count;
    delete[]buffer;
    // set output precision
    if (file_myrank == MASTER_PROCESS) {
        for (USI n = 0; n < sumdata->size(); n++) {
            if (sumdata->at(n).Item == "TIME") {
                sumdata->at(n).Type = "fixed";
                continue;
            }
            if (sumdata->at(n).Item == "NRiter" || sumdata->at(n).Item == "LSiter") {
                sumdata->at(n).Type = "int";
                continue;
            }
            else {
                sumdata->at(n).Type = "float";
                continue;
            }
        }

        PrintInfo(dir, filename, -1);
    }
}



void CriticalInfo::Setup()
{
    // Allocate memory
    const USI maxRowNum = 1000;

    Sumdata.push_back(SumItem("TIME", "-", "DAY", "fixed", maxRowNum));
    Sumdata.push_back(SumItem("dt/NRiter", "-", "DAY/-", "fixed", maxRowNum));
    Sumdata.push_back(SumItem("LSiter", "-", "-", "int", maxRowNum));
    Sumdata.push_back(SumItem("dPBmax", "-", "-", "float", maxRowNum));
    Sumdata.push_back(SumItem("dPWmax", "-", "-", "float", maxRowNum));
    Sumdata.push_back(SumItem("dTmax", "-", "-", "float", maxRowNum));
    Sumdata.push_back(SumItem("eVmax", "-", "-", "float", maxRowNum));
    Sumdata.push_back(SumItem("dSmax", "-", "-", "float", maxRowNum));
    Sumdata.push_back(SumItem("dNmax", "-", "-", "float", maxRowNum));
    Sumdata.push_back(SumItem("CFL", "-", "-", "float", maxRowNum));

}

void CriticalInfo::SetVal(const OCPControl& ctrl, const OCPNRsuite& NRs)
{

    // for NR step
    const auto& LS  = NRs.GetIterNRLS();
    const auto& dPB = NRs.DPBmaxNR();
    const auto& dPW = NRs.DPWmaxNR();
    const auto& dT  = NRs.DTmaxNR();
    const auto& eV  = NRs.EVmaxNR();
    const auto& dS  = NRs.DSmaxNR();
    const auto& dN  = NRs.DNmaxNR();

    USI n = 0;
    for (INT i = 0; i < dPB.size(); i++) {
        // Time
        Sumdata[n++].val.push_back(-1);
        // Time step size
        Sumdata[n++].val.push_back(-(i + 1));
        // LS
        Sumdata[n++].val.push_back(LS[i]);
        // dPBmax
        Sumdata[n++].val.push_back(dPB[i]);
        // dPWmax
        Sumdata[n++].val.push_back(dPW[i]);
        // dTmax
        Sumdata[n++].val.push_back(dT[i]);
        // eVmax
        Sumdata[n++].val.push_back(eV[i]);
        // dSmax
        Sumdata[n++].val.push_back(dS[i]);
        // dNmax
        Sumdata[n++].val.push_back(dN[i]);
        // CFL
        Sumdata[n++].val.push_back(0.0);

        n = 0;
    }

    // for time step
    // Time
    Sumdata[n++].val.push_back(ctrl.time.GetCurrentTime());
    // Time step
    Sumdata[n++].val.push_back(ctrl.time.GetLastDt());
    // LS
    Sumdata[n++].val.push_back(NRs.GetIterLS());
    // dPBmax
    Sumdata[n++].val.push_back(NRs.DPBmaxT());
    // dPWmax
    Sumdata[n++].val.push_back(NRs.DPWmaxT());
    // dTmax
    Sumdata[n++].val.push_back(NRs.DTmaxT());
    // eVmax
    Sumdata[n++].val.push_back(NRs.EVmaxT());
    // dSmax
    Sumdata[n++].val.push_back(NRs.DSmaxT());
    // dNmax
    Sumdata[n++].val.push_back(NRs.DNmaxT());
    // CFL
    Sumdata[n++].val.push_back(NRs.GetMaxCFL());
}

void CriticalInfo::PrintFastReview(const string& dir, const string& filename, const OCP_INT& rank, const ItersInfo& iters) const
{
    string   FileOut;
    if (rank >= 0) {
        FileOut = dir + "proc_" + to_string(rank) + "_FastReview.out";
    }
    else {
        FileOut = dir + "FastReview.out";
    }

    const USI ns      = 12;
    ofstream outF(FileOut);
    if (!outF.is_open()) {
        OCP_ABORT("Can not open " + FileOut);
    }

    outF << "FastReview OF RUN " + dir + filename
        << " -- " << iters.GetNumTimeStep() << " time step "
        << " -- " << iters.GetNRt() << " (+" << iters.GetNRwt() << ") NR step"
        << " : " << iters.GetLSt() << " (+" << iters.GetLSwt() << ") LS step\n";

    // Item
    for (const auto& v : Sumdata) {
        outF << setw(ns) << v.Item;
    }
    outF << "\n";
    // Unit
    for (const auto& v : Sumdata) {
        outF << setw(ns) << v.Unit;
    }
    outF << "\n";
    // Value
    const OCP_USI num = Sumdata[0].val.size();
    for (OCP_USI n = 0; n < num; n++) {
        for (const auto& v : Sumdata) {

            outF.unsetf(ios_base::fixed);
            outF.unsetf(ios_base::scientific);

            if (rank < 0 && v.Item == "TIME" && v.val[n] < 0) {
                outF << setw(ns) << "-";
                continue;
            }

            // for NR STEP
            if (v.Item == "dt/NRiter" && v.val[n] < 0) {
                outF << fixed << setprecision(0) << setw(ns) << -v.val[n];
                continue;
            }

            if (v.Type == "fixed") {
                outF << fixed << setprecision(3);
            }
            else if (v.Type == "int") {
                outF << fixed << setprecision(0);
            }
            else if (v.Type == "float") {
                outF << scientific << setprecision(3);
            }

            outF << setw(ns) << v.val[n];
        }
        outF << "\n";
    }
    outF.close();
}


/// Combine all files into 1 by Master process
//void CriticalInfo::PostProcess(const string& dir, const string& filename, const OCP_INT& numproc, const ItersInfo& iters) const
//{
//
//    vector<SumItem>* sumdata = const_cast<vector<SumItem>*>(&Sumdata);
//    sumdata->clear();
//
//    OCP_USI         rowNum = 0;
//    vector<string>  buffer;
//
//    for (USI p = 0; p < numproc; p++) {
//        const string myFile = dir + "proc_" + to_string(p) + "_FastReview.out";
//        ifstream ifs(myFile, ios::in);
//
//        if (!ifs) {
//            OCP_MESSAGE("Trying to open file: " << (myFile));
//            OCP_ABORT("Failed to open the input file!");
//        }
//
//        if (p == 0) {
//            // Get time steps
//            ReadLine(ifs, buffer, OCP_FALSE);
//            for (USI i = 0; i < buffer.size(); i++) {
//                if (buffer[i] == "--") {
//                    rowNum += stoi(buffer[i + 1]);
//                    continue;
//                }
//            }
//            // Get Item
//            ReadLine(ifs, buffer, OCP_FALSE);
//            for (USI i = 0; i < buffer.size(); i++) {
//                sumdata->push_back(SumItem(buffer[i], rowNum));
//            }
//            // Get Unit
//            ReadLine(ifs, buffer, OCP_FALSE);
//            for (USI i = 0; i < buffer.size(); i++) {
//                sumdata->at(i).Unit = buffer[i];
//            }
//            // Get val
//            for (OCP_USI n = 0; n < rowNum; n++) {
//                ReadLine(ifs, buffer, OCP_FALSE);
//                for (USI i = 0; i < buffer.size(); i++) {
//                    sumdata->at(i).val.push_back(stod(buffer[i]));
//                }
//            }
//        }
//        else {
//            // skip first three lines
//            ReadLine(ifs, buffer, OCP_FALSE);
//            ReadLine(ifs, buffer, OCP_FALSE);
//            ReadLine(ifs, buffer, OCP_FALSE);
//            // Get val
//            for (OCP_USI n = 0; n < rowNum; n++) {
//                ReadLine(ifs, buffer, OCP_FALSE);
//                for (USI i = 2; i < buffer.size(); i++) {
//                    sumdata->at(i).val[n] = max(sumdata->at(i).val[n], static_cast<OCP_DBL>(stod(buffer[i])));
//                }
//            }
//        }
//        ifs.close();
//        if (remove(myFile.c_str()) != 0) {
//            OCP_WARNING("Failed to delete " + myFile);
//        }
//    }
//    for (auto& v : *sumdata) {
//        if (v.Item == "TIME" || v.Item == "dt") {
//            v.Type = "fixed";
//        }
//        else if (v.Item == "LS") {
//            v.Type = "int";
//        }
//        else {
//            v.Type = "float";
//        }
//    }
//    PrintFastReview(dir, filename, -1, iters);
//}


// By Communication
void CriticalInfo::PostProcess(const string& dir, const string& filename, const OCP_INT& numproc, const ItersInfo& iters) const
{

    vector<SumItem>* sumdata = const_cast<vector<SumItem>*>(&Sumdata);

    OCP_USI         rowNum = Sumdata[0].val.size();
    double* buffer = new double[rowNum];
    for (OCP_USI i = 2; i < sumdata->size(); i++) {
        MPI_Reduce(sumdata->at(i).val.data(), buffer, rowNum, MPI_DOUBLE, MPI_MAX, MASTER_PROCESS, file_myComm);
        if (file_myrank == MASTER_PROCESS) {
            memcpy(sumdata->at(i).val.data(), buffer, sizeof(double) * rowNum);
        }
    }
    delete[]buffer;

    if (file_myrank == MASTER_PROCESS) {
        for (auto& v : *sumdata) {
            if (v.Item == "TIME" || v.Item == "dt") {
                v.Type = "fixed";
            }
            else if (v.Item == "LS") {
                v.Type = "int";
            }
            else {
                v.Type = "float";
            }
        }
        PrintFastReview(dir, filename, -1, iters);
    }
}



void OutGridVarSet::Setup(const OutGridParam& param, const Bulk& bk)
{
    DEPTH     = param.DEPTH;
    PRE       = param.PRE;
    PHASEP    = param.PHASEP;
    COMPM     = param.COMPM;
    SOIL      = param.SOIL;
    SGAS      = param.SGAS;
    SWAT      = param.SWAT;
    DENO      = param.DENO;
    DENG      = param.DENG;
    DENW      = param.DENW;
    KRO       = param.KRO;
    KRG       = param.KRG;
    KRW       = param.KRW;
    BOIL      = param.BOIL;
    BGAS      = param.BGAS;
    BWAT      = param.BWAT;
    VOIL      = param.VOIL;
    VGAS      = param.VGAS;
    VWAT      = param.VWAT;
    XMF       = param.XMF;
    YMF       = param.YMF;
    PCW       = param.PCW;
    CO2       = param.CO2;
    SATNUM    = param.SATNUM;
    PERMX     = param.PERMX;
    PERMY     = param.PERMY;
    PERMZ     = param.PERMZ;
    DSAT      = param.DSAT;
    DP        = param.DP;
    CSFLAG    = param.CSFLAG;
    ITERNRDDM = param.ITERNRDDM;
    ITERLSDDM = param.ITERLSDDM;
    TIMELSDDM = param.TIMELSDDM;

    nc = bk.GetComNum();
    np = bk.GetPhaseNum();

    // correct wrong output request
    if (!bk.IfOilExist()) {
        SOIL = OCP_FALSE;
        DENO = OCP_FALSE;
        BOIL = OCP_FALSE;
        KRO = OCP_FALSE;
        VOIL = OCP_FALSE;
        XMF = OCP_FALSE;
    }
    if (!bk.IfGasExist()) {
        SGAS = OCP_FALSE;
        DENG = OCP_FALSE;
        BGAS = OCP_FALSE;
        KRG = OCP_FALSE;
        VGAS = OCP_FALSE;
        YMF = OCP_FALSE;
    }
    if (!bk.IfWatExist()) {
        SWAT = OCP_FALSE;
        DENW = OCP_FALSE;
        BWAT = OCP_FALSE;
        KRW = OCP_FALSE;
        VWAT = OCP_FALSE;
    }

    bgpnum = 0;
    if (DEPTH)    bgpnum++;
    if (PRE)      bgpnum++;
    if (COMPM)    bgpnum += nc;
    if (PHASEP)   bgpnum += np;
    if (SOIL)     bgpnum++;
    if (SGAS)     bgpnum++;
    if (SWAT)     bgpnum++;
    //if (DENO)     bgpnum++;
    //if (DENG)     bgpnum++;
    //if (DENW)     bgpnum++;
    //if (KRO)      bgpnum++;
    //if (KRG)      bgpnum++;
    //if (KRW)      bgpnum++;
    //if (BOIL)     bgpnum++;
    //if (BGAS)     bgpnum++;
    //if (BWAT)     bgpnum++;
    //if (VOIL)     bgpnum++;
    //if (VGAS)     bgpnum++;
    //if (VWAT)     bgpnum++;
    //if (XMF)      bgpnum++;
    //if (YMF)      bgpnum++;
    //if (PCW)      bgpnum++;
    if (CO2)       bgpnum++;
    if (SATNUM)    bgpnum++;
    if (PERMX)     bgpnum++;
    if (PERMY)     bgpnum++;
    if (DSAT)      bgpnum += np;
    if (DP)        bgpnum += np;
    if (CSFLAG)    bgpnum++;
    if (ITERNRDDM) bgpnum++;
    if (ITERLSDDM) bgpnum++;
    if (TIMELSDDM) bgpnum++;
    //if (PERMY)    bgpnum++;
    //if (PERMZ)    bgpnum++;
}



//void Out4RPT::InputParam(const OutputRPTParam& RPTparam)
//{
//    useRPT = RPTparam.useRPT;
//    if (!useRPT) return;
//
//    bgp.Setup(RPTparam.bgp, rs.bulk);
//}

//void Out4RPT::Setup(const string& dir, const Reservoir& rs)
//{
//    if (!useRPT) return;
//
//    string   FileOut = dir + "RPT.out";
//    ofstream outF(FileOut);
//    if (!outF.is_open()) {
//        OCP_ABORT("Can not open " + FileOut);
//    }
//    outF.close();
//
//    const Grid& initGrid = rs.grid;
//
//    nx       = initGrid.nx;
//    ny       = initGrid.ny;
//    numGrid  = initGrid.numGrid;
//    IJKspace = initGrid.numDigutIJK;
//}

//void Out4RPT::PrintRPT(const string&    dir,
//                       const Reservoir& rs,
//                       const OCP_DBL&   days) const
//{
//
//    if (!useRPT) return;
//
//    string   FileOut = dir + "RPT.out";
//    ofstream outRPT;
//    outRPT.open(FileOut, ios::app);
//    if (!outRPT.is_open()) {
//        OCP_ABORT("Can not open " + FileOut);
//    }
//
//    const Grid& initGrid = rs.grid;
//    const Bulk& bulk     = rs.bulk;
//
//    const USI              np     = bulk.numPhase;
//    const USI              nc     = bulk.numCom;
//    const USI              OIndex = bulk.phase2Index[OIL];
//    const USI              GIndex = bulk.phase2Index[GAS];
//    const USI              WIndex = bulk.phase2Index[WATER];
//    const vector<GB_Pair>& g2bp   = initGrid.map_All2Act;
//
//    outRPT << OCP_SEP02(50) << "\n";
//
//    // Well Info
//    USI numWell = rs.allWells.GetWellNum();
//    outRPT << "Well Information"
//           << "                    ";
//    outRPT << fixed << setprecision(3) << days << "  DAYS"
//           << "\n";
//    // INJ
//    for (USI w = 0; w < numWell; w++) {
//        if (rs.allWells.wells[w].opt.type == INJ) {
//            outRPT << "-------------------------------------"
//                   << "\n";
//            outRPT << rs.allWells.wells[w].name << "   " << w << "   "
//                   << rs.allWells.wells[w].depth << " (feet)     ";
//            outRPT << rs.allWells.wells[w].I << "   " << rs.allWells.wells[w].J << "\n";
//
//            if (rs.allWells.wells[w].opt.state == OPEN) {
//                outRPT << "OPEN\t" << rs.allWells.wells[w].WGIR << " (MSCF/DAY)\t"
//                       << rs.allWells.wells[w].WWIR << " (STB/DAY)"
//                       << "\n";
//            } else {
//                outRPT << "SHUTIN"
//                       << "\n";
//            }
//            // perf
//            for (USI p = 0; p < rs.allWells.wells[w].numPerf; p++) {
//                outRPT << "perf" << p << "   " << rs.allWells.wells[w].perf[p].I
//                       << "   " << rs.allWells.wells[w].perf[p].J << "   "
//                       << rs.allWells.wells[w].perf[p].K << "   "
//                       << rs.allWells.wells[w].perf[p].depth << "   ";
//                if (rs.allWells.wells[w].perf[p].state == OPEN) {
//                    outRPT << "OPEN";
//                } else {
//                    outRPT << "SHUTIN";
//                }
//                outRPT << "   " << rs.allWells.wells[w].perf[p].location << "\n";
//            }
//        }
//    }
//    // PROD
//    for (USI w = 0; w < numWell; w++) {
//        if (rs.allWells.wells[w].opt.type == PROD) {
//            outRPT << "-------------------------------------"
//                   << "\n";
//            outRPT << rs.allWells.wells[w].name << "   " << w << "   "
//                   << rs.allWells.wells[w].depth << " (feet)     ";
//            outRPT << rs.allWells.wells[w].I << "   " << rs.allWells.wells[w].J << "\n";
//
//            if (rs.allWells.wells[w].opt.state == OPEN) {
//                outRPT << "OPEN\t" << rs.allWells.wells[w].WOPR << " (STB/DAY)\t"
//                       << rs.allWells.wells[w].WGPR << " (MSCF/DAY)\t"
//                       << rs.allWells.wells[w].WWPR << " (STB/DAY)"
//                       << "\n";
//            } else {
//                outRPT << "SHUTIN"
//                       << "\n";
//            }
//            // perf
//            for (USI p = 0; p < rs.allWells.wells[w].numPerf; p++) {
//                outRPT << "perf" << p << "   " << rs.allWells.wells[w].perf[p].I
//                       << "   " << rs.allWells.wells[w].perf[p].J << "   "
//                       << rs.allWells.wells[w].perf[p].K << "   "
//                       << rs.allWells.wells[w].perf[p].depth << "   ";
//                if (rs.allWells.wells[w].perf[p].state == OPEN) {
//                    outRPT << "OPEN";
//                } else {
//                    outRPT << "SHUTIN";
//                }
//                outRPT << "   " << rs.allWells.wells[w].perf[p].location << "\n";
//            }
//        }
//    }
//
//    outRPT << "\n\n";
//
//    static OCP_BOOL flag = OCP_FALSE;
//    // Print once
//    if (flag) {
//        PrintRPT_Scalar(outRPT, "DX : feet", days, &initGrid.dx[0], 1, g2bp, OCP_FALSE);
//        PrintRPT_Scalar(outRPT, "DY : feet", days, &initGrid.dy[0], 1, g2bp, OCP_FALSE);
//        PrintRPT_Scalar(outRPT, "DZ : feet", days, &initGrid.dz[0], 1, g2bp, OCP_FALSE);
//        PrintRPT_Scalar(outRPT, "Depth : feet", days, &initGrid.depth[0], 1, g2bp,
//                        OCP_FALSE);
//        PrintRPT_Scalar(outRPT, "PERMX : MDarcy", days, &initGrid.kx[0], 1, g2bp,
//                        OCP_FALSE);
//        PrintRPT_Scalar(outRPT, "PERMY : MDarcy", days, &initGrid.ky[0], 1, g2bp,
//                        OCP_FALSE);
//        PrintRPT_Scalar(outRPT, "PERMZ : MDarcy", days, &initGrid.kz[0], 1, g2bp,
//                        OCP_FALSE);
//        flag = OCP_FALSE;
//    }
//
//    // PRESSURE
//    if (bgp.PRE) {
//        PrintRPT_Scalar(outRPT, "PRESSURE : psia", days, &bulk.P[0], 1, g2bp, OCP_TRUE);
//    }
//
//    // DENSITY of OIL
//    if (bgp.DENO && bulk.oil) {
//        PrintRPT_Scalar(outRPT, "DENO : lb/ft3", days, &bulk.rho[OIndex], np, g2bp,
//                        OCP_TRUE);
//    }
//    outRPT << endl;
//
//    // DENSITY of GAS
//    if (bgp.DENG && bulk.gas) {
//        PrintRPT_Scalar(outRPT, "DENG : lb/ft3", days, &bulk.rho[GIndex], np, g2bp,
//                        OCP_TRUE);
//    }
//
//    // DENSITY of WATER
//    if (bgp.DENW && bulk.water) {
//        PrintRPT_Scalar(outRPT, "DENW : lb/ft3", days, &bulk.rho[WIndex], np, g2bp,
//                        OCP_TRUE);
//    }
//
//    // SATURATION of OIL
//    if (bgp.SOIL && bulk.oil) {
//        PrintRPT_Scalar(outRPT, "SOIL         ", days, &bulk.S[OIndex], np, g2bp,
//                        OCP_TRUE);
//    }
//
//    // SATURATION of GAS
//    if (bgp.SGAS && bulk.gas) {
//        PrintRPT_Scalar(outRPT, "SGAS         ", days, &bulk.S[GIndex], np, g2bp,
//                        OCP_TRUE);
//    }
//
//    // SATURATION of WATER
//    if (bgp.SWAT && bulk.water) {
//        PrintRPT_Scalar(outRPT, "SWAT         ", days, &bulk.S[WIndex], np, g2bp,
//                        OCP_TRUE);
//    }
//
//    // Relative Permeability of OIL
//    if (bgp.KRO && bulk.oil) {
//        PrintRPT_Scalar(outRPT, "KRO          ", days, &bulk.kr[OIndex], np, g2bp,
//                        OCP_TRUE);
//    }
//
//    // Relative Permeability of GAS
//    if (bgp.KRG && bulk.gas) {
//        PrintRPT_Scalar(outRPT, "KRG          ", days, &bulk.kr[GIndex], np, g2bp,
//                        OCP_TRUE);
//    }
//
//    // Relative Permeability of WATER
//    if (bgp.KRW && bulk.water) {
//        PrintRPT_Scalar(outRPT, "KRW          ", days, &bulk.kr[WIndex], np, g2bp,
//                        OCP_TRUE);
//    }
//
//    // Molar Density of OIL
//    if (bgp.BOIL && bulk.oil && bulk.IfUseEoS()) {
//        PrintRPT_Scalar(outRPT, "BOIL : lb-M/rb", days, &bulk.xi[OIndex], np, g2bp,
//                        OCP_TRUE);
//    }
//
//    // Molar Density of GAS
//    if (bgp.BGAS && bulk.gas && bulk.IfUseEoS()) {
//        PrintRPT_Scalar(outRPT, "BGAS : lb-M/rb", days, &bulk.xi[GIndex], np, g2bp,
//                        OCP_TRUE);
//    }
//
//    // Molar Density of WATER
//    if (bgp.BWAT && bulk.water) {
//        PrintRPT_Scalar(outRPT, "BWAT : lb-M/rb", days, &bulk.xi[WIndex], np, g2bp,
//                        OCP_TRUE, (CONV1 * 19.437216));
//    }
//
//    // Viscosity of OIL
//    if (bgp.VOIL && bulk.oil) {
//        PrintRPT_Scalar(outRPT, "VOIL : lb-M/rb", days, &bulk.mu[OIndex], np, g2bp,
//                        OCP_TRUE);
//    }
//
//    // Viscosity of GAS
//    if (bgp.VGAS && bulk.gas) {
//        PrintRPT_Scalar(outRPT, "VGAS : lb-M/rb", days, &bulk.mu[GIndex], np, g2bp,
//                        OCP_TRUE);
//    }
//
//    // Viscosity of WATER
//    if (bgp.VWAT && bulk.water) {
//        PrintRPT_Scalar(outRPT, "VWAT : lb-M/rb", days, &bulk.mu[WIndex], np, g2bp,
//                        OCP_TRUE);
//    }
//
//    // liquid component mole fractions.
//    if (bgp.XMF && bulk.IfUseEoS()) {
//        for (USI i = 0; i < nc - 1; i++) {
//            PrintRPT_Scalar(outRPT, "XMF : Oil  " + to_string(i + 1) + "th Component",
//                            days, &bulk.xij[OIndex * nc + i], np * nc, g2bp, OCP_TRUE);
//        }
//    }
//
//    // gas component mole fractions.
//    if (bgp.YMF && bulk.IfUseEoS()) {
//        for (USI i = 0; i < nc - 1; i++) {
//            PrintRPT_Scalar(outRPT, "YMF : Gas  " + to_string(i + 1) + "th Component",
//                            days, &bulk.xij[GIndex * nc + i], np * nc, g2bp, OCP_TRUE);
//        }
//    }
//
//    // Po - Pw
//    if (bgp.PCW) {
//        PrintRPT_Scalar(outRPT, "PCW : psia  ", days, &bulk.Pc[WIndex], np, g2bp,
//                        OCP_TRUE);
//
//        // PrintRPT_Scalar(outRPT, "PPCW : psia ", days,
//        // &rs.optFeatures.scalePcow.scaleVal[0],
//        //                 1, g2bp, OCP_TRUE, 103.053 - 0.116);
//    }
//
//    outRPT.close();
//}

//void Out4RPT::GetIJKGrid(USI& i, USI& j, USI& k, const OCP_ULL& n) const
//{
//    // i,j,k begin from 1
//    // n must be the index of grids instead bulks
//    k = n / (nx * ny) + 1;
//    j = (n - (k - 1) * nx * ny) / nx + 1;
//    i = n - (k - 1) * nx * ny - (j - 1) * nx + 1;
//}


static const INT timeInfoLen = 256;

void Out4VTK::Setup(const OutputVTKParam& VTKParam, const string& dir, const Reservoir& rs)
{
    useVTK = VTKParam.useVTK;
    if (!useVTK) return;

    bgp.Setup(VTKParam.bgp, rs.bulk);
    out4vtk.Setup(VTKParam.bgp.ASCII, VTKParam.bgp.DOUBLE);

    Initialize(dir, rs);
}


void Out4VTK::Initialize(const string& dir, const Reservoir& rs)
{
    if (!useVTK) return;

    timeInfo = new OCP_CHAR[timeInfoLen];

    // output the gloabl index of grids belonging to current domain
    const Domain& doman = rs.domain;
    if (doman.global_numproc > 1) {
        myFile = dir + "proc" + to_string(doman.global_rank) + "_vtktmp.out";
        ofstream outF(myFile, ios::out | ios::binary);
        if (!outF.is_open()) {
            OCP_WARNING("Can not open " + myFile);
        }

        const OCP_USI numInteriorGrid = doman.GetNumGridInterior();
        outF.write((const OCP_CHAR*)&numInteriorGrid, sizeof(numInteriorGrid));
        outF.write((const OCP_CHAR*)&doman.GetGrid()[0], numInteriorGrid * sizeof(doman.GetGrid()[0]));
        outF.close();
    }
    else {
        myFile = dir + "main_vtktmp.out";
        ofstream outF(myFile, ios::out | ios::binary);
        if (!outF.is_open()) {
            OCP_WARNING("Can not open " + myFile);
        }
        outF.close();
    }    
}


void Out4VTK::PrintVTK(const Reservoir& rs, const OCPControl& ctrl) const
{
    if (!useVTK)         return;

    countPrint++;

    const BulkVarSet& bvs = rs.bulk.vs;

    const auto nb     = bvs.nbI;
    const auto np     = bvs.np;
    const auto nc     = bvs.nc;
    const auto OIndex = bvs.o;
    const auto GIndex = bvs.g;
    const auto WIndex = bvs.w;

    ofstream outF(myFile, ios::app | ios::binary);
    vector<OCP_DBL> tmpV(nb);
    // output    
    // output time info
    {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3) << ctrl.GetOCPFile() << ctrl.time.GetCurrentTime();
        std::string str = oss.str() + TIMEUNIT + "_";
        strncpy(timeInfo, str.c_str(), timeInfoLen);
        outF.write(timeInfo, timeInfoLen);
    }

    if (bgp.bgpnum == 0) {
        outF.close();
        return;
    }

    // output depth 
    if (bgp.DEPTH) {
        for (OCP_USI n = 0; n < nb; n++)
            tmpV[n] = bvs.depth[n];
        outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
    }
    // output physical variables
    if (bgp.PRE) {
        for (OCP_USI n = 0; n < nb; n++)
            tmpV[n] = bvs.P[n];
        outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
    }
    if (bgp.COMPM) {
        for (USI i = 0; i < nc; i++) {
            for (OCP_USI n = 0; n < nb; n++)
                tmpV[n] = bvs.Ni[n * nc + i];
            outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
        }
    }
    if (bgp.PHASEP) {
        for (USI j = 0; j < np; j++) {
            for (OCP_USI n = 0; n < nb; n++)
                tmpV[n] = bvs.Pj[n * np + j];
            outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
        }
    }
    if (bgp.SOIL) {
        for (OCP_USI n = 0; n < nb; n++)
            tmpV[n] = bvs.S[n * np + OIndex];
        outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
    }
    if (bgp.SGAS) {
        for (OCP_USI n = 0; n < nb; n++)
            tmpV[n] = bvs.S[n * np + GIndex];
        outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
    }
    if (bgp.SWAT) {
        for (OCP_USI n = 0; n < nb; n++)
            tmpV[n] = bvs.S[n * np + WIndex];
        outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
    }

    // add CO2 concentration for SPE11
    if (bgp.CO2) {
        for (OCP_USI n = 0; n < nb; n++)
            tmpV[n] = bvs.rho[n * np + WIndex] * bvs.xij[(n * np + WIndex) * nc + GIndex];
        outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
    }

    // add SAT region for SPE11
    if (bgp.SATNUM) {
        for (OCP_USI n = 0; n < nb; n++)
            tmpV[n] = static_cast<OCP_DBL>(rs.bulk.SATm.GetSATNUM(n));
        outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
    }

    if (bgp.PERMX) {
        for (OCP_USI n = 0; n < nb; n++)
            tmpV[n] = bvs.rockKx[n];
        outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
    }

    if (bgp.PERMY) {
        for (OCP_USI n = 0; n < nb; n++)
            tmpV[n] = bvs.rockKy[n];
        outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
    }

    OCP_DBL dt = ctrl.time.GetLastDt();
    dt = dt < 1E-8 ? 1 : dt;

    if (bgp.DSAT) {
        for (USI j = 0; j < np; j++) {
            for (OCP_USI n = 0; n < nb; n++)
                tmpV[n] = fabs(bvs.S[n * np + j] - bvs.lS[n * np + j]);
            outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
        }
    }
    if (bgp.DP) {
        for (USI j = 0; j < np; j++) {
            for (OCP_USI n = 0; n < nb; n++)
                tmpV[n] = fabs(bvs.Pj[n * np + j] - bvs.lPj[n * np + j]);
            outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
        }
    }
    if (bgp.CSFLAG) {
        OCP_DBL flag = -1.0;
        if (rs.domain.cs_group_global_rank_for_output.size() > 1) {
            flag = static_cast<OCP_DBL>(*rs.domain.cs_group_global_rank_for_output.begin());
        }
        fill(tmpV.begin(), tmpV.end(), flag);
        outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
    }
    if (bgp.ITERNRDDM) {
        fill(tmpV.begin(), tmpV.end(), static_cast<OCP_DBL>(OCPITER_NR_DDM));
        outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
    }
    if (bgp.ITERLSDDM) {
        fill(tmpV.begin(), tmpV.end(), static_cast<OCP_DBL>(OCPITER_LS_DDM));
        outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
    }
    if (bgp.TIMELSDDM) {
        fill(tmpV.begin(), tmpV.end(), OCPTIME_LSOLVER_DDM);
        outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
    }
         
    outF.close();
}


void Out4VTK::PostProcess(const string& dir, const string& filename, const OCP_INT& numproc) const
{
    if (numproc > 1)  PostProcessP(dir, filename, numproc);
    else              PostProcessS(dir, filename);
}


void Out4VTK::PostProcessP(const string& dir, const string& filename, const OCP_INT& numproc) const
{
    if (!useVTK) return;

    // Input Points
    const string srcFile = dir + "TSTEP.vtk";
    numGrid = out4vtk.Init(dir, srcFile, "RUN of " + dir + filename);

    // Input cell values    
    vector<vector<OCP_DBL>> gridVal;        // each row is on a node of TSTEP in turn
    vector<OCP_ULL>         global_index;
    vector<OCP_USI>         mypart(numGrid);
    OCP_DBL*                workPtr;
    vector<OCP_DBL>         tmpVal;
    

    vector<string>          timeSeries;

    for (USI p = 0; p < numproc; p++) {
        const string myfile = dir + "proc" + to_string(p) + "_vtktmp.out";
        ifstream inV(myfile, ios::in | ios::binary);
        if (!inV.is_open()) {
            OCP_WARNING("Can not open " + myfile);
        }
        
        OCP_USI numGridLoc = 0;
        inV.read((OCP_CHAR*)(&numGridLoc), sizeof(numGridLoc));
        global_index.resize(numGridLoc);
        inV.read((OCP_CHAR*)(&global_index[0]), sizeof(global_index[0]) * numGridLoc);

        for (OCP_USI n = 0; n < numGridLoc; n++)
            mypart[global_index[n]] = p;


		USI index = 0;
		while (index < countPrint) {
			// input time info
			inV.read(timeInfo, timeInfoLen);
			if (p == 0) {
				timeSeries.push_back(string(timeInfo));
			}

			// input grid info
            if (bgp.bgpnum == 0) {
                index++;
                continue;
            }

			tmpVal.resize(numGridLoc * bgp.bgpnum);

			if (index >= gridVal.size()) {
				gridVal.push_back(vector<OCP_DBL>());
				gridVal.back().resize(bgp.bgpnum * numGrid);
			}
			workPtr = &gridVal[index][0];
			inV.read((OCP_CHAR*)(&tmpVal[0]), sizeof(tmpVal[0]) * tmpVal.size());
			const OCP_DBL* tmpVal_ptr = &tmpVal[0];

            if (bgp.DEPTH) {
				for (OCP_USI n = 0; n < numGridLoc; n++) {
					workPtr[global_index[n]] = tmpVal_ptr[n];
				}
				workPtr += numGrid;
				tmpVal_ptr += numGridLoc;
			}
			if (bgp.PRE) {
				for (OCP_USI n = 0; n < numGridLoc; n++) {
					workPtr[global_index[n]] = tmpVal_ptr[n];
				}
				workPtr += numGrid;
				tmpVal_ptr += numGridLoc;
			}
			if (bgp.COMPM) {
				for (USI i = 0; i < bgp.nc; i++) {
					for (OCP_USI n = 0; n < numGridLoc; n++) {
						workPtr[global_index[n]] = tmpVal_ptr[n];
					}
					workPtr += numGrid;
					tmpVal_ptr += numGridLoc;
				}
			}
			if (bgp.PHASEP) {
				for (USI j = 0; j < bgp.np; j++) {
					for (OCP_USI n = 0; n < numGridLoc; n++) {
						workPtr[global_index[n]] = tmpVal_ptr[n];
					}
					workPtr += numGrid;
					tmpVal_ptr += numGridLoc;
				}
			}
			if (bgp.SOIL) {
				for (OCP_USI n = 0; n < numGridLoc; n++) {
					workPtr[global_index[n]] = tmpVal_ptr[n];
				}
				workPtr += numGrid;
				tmpVal_ptr += numGridLoc;
			}
			if (bgp.SGAS) {
				for (OCP_USI n = 0; n < numGridLoc; n++) {
					workPtr[global_index[n]] = tmpVal_ptr[n];
				}
				workPtr += numGrid;
				tmpVal_ptr += numGridLoc;
			}
			if (bgp.SWAT) {
				for (OCP_USI n = 0; n < numGridLoc; n++) {
					workPtr[global_index[n]] = tmpVal_ptr[n];
				}
				workPtr += numGrid;
				tmpVal_ptr += numGridLoc;
			}
			if (bgp.CO2) {
				for (OCP_USI n = 0; n < numGridLoc; n++) {
					workPtr[global_index[n]] = tmpVal_ptr[n];
				}
				workPtr += numGrid;
				tmpVal_ptr += numGridLoc;
			}
			if (bgp.SATNUM) {
				for (OCP_USI n = 0; n < numGridLoc; n++) {
					workPtr[global_index[n]] = tmpVal_ptr[n];
				}
				workPtr += numGrid;
				tmpVal_ptr += numGridLoc;
			}
			if (bgp.PERMX) {
				for (OCP_USI n = 0; n < numGridLoc; n++) {
					workPtr[global_index[n]] = tmpVal_ptr[n];
				}
				workPtr += numGrid;
				tmpVal_ptr += numGridLoc;
			}
            if (bgp.PERMY) {
                for (OCP_USI n = 0; n < numGridLoc; n++) {
                    workPtr[global_index[n]] = tmpVal_ptr[n];
                }
                workPtr += numGrid;
                tmpVal_ptr += numGridLoc;
            }
			if (bgp.DSAT) {
				for (USI j = 0; j < bgp.np; j++) {
					for (OCP_USI n = 0; n < numGridLoc; n++) {
						workPtr[global_index[n]] = tmpVal_ptr[n];
					}
					workPtr += numGrid;
					tmpVal_ptr += numGridLoc;
				}
			}
            if (bgp.DP) {
                for (USI j = 0; j < bgp.np; j++) {
                    for (OCP_USI n = 0; n < numGridLoc; n++) {
                        workPtr[global_index[n]] = tmpVal_ptr[n];
                    }
                    workPtr += numGrid;
                    tmpVal_ptr += numGridLoc;
                }
            }
			if (bgp.CSFLAG) {
				for (OCP_USI n = 0; n < numGridLoc; n++) {
					workPtr[global_index[n]] = tmpVal_ptr[n];
				}
				workPtr += numGrid;
				tmpVal_ptr += numGridLoc;
			}
            if (bgp.ITERNRDDM) {
                for (OCP_USI n = 0; n < numGridLoc; n++) {
                    workPtr[global_index[n]] = tmpVal_ptr[n];
                }
                workPtr += numGrid;
                tmpVal_ptr += numGridLoc;
            }
            if (bgp.ITERLSDDM) {
                for (OCP_USI n = 0; n < numGridLoc; n++) {
                    workPtr[global_index[n]] = tmpVal_ptr[n];
                }
                workPtr += numGrid;
                tmpVal_ptr += numGridLoc;
            }
            if (bgp.TIMELSDDM) {
                for (OCP_USI n = 0; n < numGridLoc; n++) {
                    workPtr[global_index[n]] = tmpVal_ptr[n];
                }
                workPtr += numGrid;
                tmpVal_ptr += numGridLoc;
            }

			index++;
		}
        inV.close();
		if (remove(myfile.c_str()) != 0) {
			OCP_WARNING("Failed to delete " + myfile);
		}
    }

    // OutPut

    USI numTstep = timeSeries.size();
    if (numTstep != gridVal.size()) {
        cout << numTstep << "   " << gridVal.size() << endl;
        OCP_ABORT("Something Wrong in the temporary files");
    }

    if (numTstep == 0) {
        ofstream source(srcFile, ios::app);
        source << "\n" << VTK_CELL_DATA << " " << numGrid;
        // Ouput partition
        out4vtk.OutputCELL_DATA_SCALARS(source, "PARTITION", VTK_UNSIGNED_INT, mypart, 0, numGrid, 0);
        source.close();
    }
    else {
        for (USI t = 0; t < numTstep; t++) {
            const string dstFile = dir + timeSeries[t] + to_string(t) + ".vtk";
            ifstream source(srcFile, ios::binary);
            ofstream dest(dstFile, ios::binary);
            dest << source.rdbuf();
            source.close();

            dest << "\n" << VTK_CELL_DATA << " " << numGrid;
            // Ouput partition
            out4vtk.OutputCELL_DATA_SCALARS(dest, "PARTITION", VTK_UNSIGNED_INT, mypart, 0, numGrid, 0);

            OCP_ULL bId = 0;
            if (bgp.DEPTH) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "DEPTH", VTK_FLOAT, gridVal[t], bId, numGrid, 3);
                bId += numGrid;
            }
            if (bgp.PRE) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "PRESSURE", VTK_FLOAT, gridVal[t], bId, numGrid, 3);
                bId += numGrid;
            }
            if (bgp.COMPM) {
                for (USI i = 0; i < bgp.nc; i++) {
                    out4vtk.OutputCELL_DATA_SCALARS(dest, "COMPM-" + to_string(i), VTK_FLOAT, gridVal[t], bId, numGrid, 3);
                    bId += numGrid;
                }
            }
            if (bgp.PHASEP) {
                for (USI j = 0; j < bgp.np; j++) {
                    out4vtk.OutputCELL_DATA_SCALARS(dest, "PHASEP-" + to_string(j), VTK_FLOAT, gridVal[t], bId, numGrid, 3);
                    bId += numGrid;
                }
            }
            if (bgp.SOIL) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "SOIL", VTK_FLOAT, gridVal[t], bId, numGrid, 6);
                bId += numGrid;
            }
            if (bgp.SGAS) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "SGAS", VTK_FLOAT, gridVal[t], bId, numGrid, 6);
                bId += numGrid;
            }
            if (bgp.SWAT) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "SWAT", VTK_FLOAT, gridVal[t], bId, numGrid, 6);
                bId += numGrid;
            }
            if (bgp.CO2) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "CO2", VTK_FLOAT, gridVal[t], bId, numGrid, 6);
                bId += numGrid;
            }
            if (bgp.SATNUM) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "SATNUM", VTK_UNSIGNED_INT, gridVal[t], bId, numGrid, 0);
                bId += numGrid;
            }
            if (bgp.PERMX) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "PERMX", VTK_FLOAT, gridVal[t], bId, numGrid, 0);
                bId += numGrid;
            }
            if (bgp.PERMY) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "PERMY", VTK_FLOAT, gridVal[t], bId, numGrid, 0);
                bId += numGrid;
            }
            if (bgp.DSAT) {
                for (USI j = 0; j < bgp.np; j++) {
                    out4vtk.OutputCELL_DATA_SCALARS(dest, "DS-" + to_string(j), VTK_FLOAT, gridVal[t], bId, numGrid, 3);
                    bId += numGrid;
                }
            }
            if (bgp.DP) {
                for (USI j = 0; j < bgp.np; j++) {
                    out4vtk.OutputCELL_DATA_SCALARS(dest, "DP-" + to_string(j), VTK_FLOAT, gridVal[t], bId, numGrid, 3);
                    bId += numGrid;
                }
            }
            if (bgp.CSFLAG) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "CSFLAG", VTK_INT, gridVal[t], bId, numGrid, 0);
                bId += numGrid;
            }
            if (bgp.ITERNRDDM) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "ITERNRDDM", VTK_INT, gridVal[t], bId, numGrid, 0);
                bId += numGrid;
            }
            if (bgp.ITERLSDDM) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "ITERLSDDM", VTK_INT, gridVal[t], bId, numGrid, 0);
                bId += numGrid;
            }
            if (bgp.TIMELSDDM) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "TIMELSDDM", VTK_FLOAT, gridVal[t], bId, numGrid, 0);
                bId += numGrid;
            }

            dest.close();
        }
        if (remove(srcFile.c_str()) != 0) {
            OCP_WARNING("Failed to delete " + srcFile);
        }
    }
}


void Out4VTK::PostProcessS(const string& dir, const string& filename) const
{
    if (!useVTK) return;

    const string srcFile = dir + "TSTEP.vtk";
    numGrid = out4vtk.Init(dir, srcFile, "RUN of " + dir + filename);
    
    ifstream inV(myFile, ios::in | ios::binary);
    if (!inV.is_open()) {
        OCP_WARNING("Can not open " + myFile);
    }   

    vector<OCP_DBL> tmpVal;
    USI             index = 0; 
    while (index < countPrint) {
        // input time info
        inV.read(timeInfo, timeInfoLen);
        string tmpS(timeInfo);

        // input grid info
        if (bgp.bgpnum == 0) {
            index++;
            continue;
        }

        tmpVal.resize(bgp.bgpnum * numGrid);

        inV.read((OCP_CHAR*)(&tmpVal[0]), sizeof(tmpVal[0]) * tmpVal.size());

        const string dstFile = dir + tmpS + to_string(index++) + ".vtk";
        ifstream source(srcFile, ios::binary);
        ofstream dest(dstFile, ios::binary);
        dest << source.rdbuf();
        source.close();

        dest << "\n" << VTK_CELL_DATA << " " << numGrid;
        OCP_ULL bId = 0;
        if (bgp.DEPTH) {
            out4vtk.OutputCELL_DATA_SCALARS(dest, "DEPTH", VTK_FLOAT, tmpVal, bId, numGrid, 3);
            bId += numGrid;
        }
        if (bgp.PRE) {
            out4vtk.OutputCELL_DATA_SCALARS(dest, "PRESSURE", VTK_FLOAT, tmpVal, bId, numGrid, 3);
            bId += numGrid;
        }
        if (bgp.COMPM) {
            for (USI i = 0; i < bgp.nc; i++) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "COMPM-" + to_string(i), VTK_FLOAT, tmpVal, bId, numGrid, 3);
                bId += numGrid;
            }
        }
        if (bgp.PHASEP) {
            for (USI j = 0; j < bgp.np; j++) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "PHASEP-" + to_string(j), VTK_FLOAT, tmpVal, bId, numGrid, 3);
                bId += numGrid;
            }
        }
        if (bgp.SOIL) {
            out4vtk.OutputCELL_DATA_SCALARS(dest, "SOIL", VTK_FLOAT, tmpVal, bId, numGrid, 6);
            bId += numGrid;
        }
        if (bgp.SGAS) {
            out4vtk.OutputCELL_DATA_SCALARS(dest, "SGAS", VTK_FLOAT, tmpVal, bId, numGrid, 6);
            bId += numGrid;
        }
        if (bgp.SWAT) {
            out4vtk.OutputCELL_DATA_SCALARS(dest, "SWAT", VTK_FLOAT, tmpVal, bId, numGrid, 6);
            bId += numGrid;
        }
        if (bgp.CO2) {
            out4vtk.OutputCELL_DATA_SCALARS(dest, "CO2", VTK_FLOAT, tmpVal, bId, numGrid, 6);
            bId += numGrid;
        }
        if (bgp.SATNUM) {
            out4vtk.OutputCELL_DATA_SCALARS(dest, "SATNUM", VTK_UNSIGNED_INT, tmpVal, bId, numGrid, 0);
            bId += numGrid;
        }
        if (bgp.PERMX) {
            out4vtk.OutputCELL_DATA_SCALARS(dest, "PERMX", VTK_FLOAT, tmpVal, bId, numGrid, 0);
            bId += numGrid;
        }
        if (bgp.PERMY) {
            out4vtk.OutputCELL_DATA_SCALARS(dest, "PERMY", VTK_FLOAT, tmpVal, bId, numGrid, 0);
            bId += numGrid;
        }
        if (bgp.DSAT) {
            for (USI j = 0; j < bgp.np; j++) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "DS-" + to_string(j), VTK_FLOAT, tmpVal, bId, numGrid, 3);
                bId += numGrid;
            }
        }
        if (bgp.DP) {
            for (USI j = 0; j < bgp.np; j++) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "DP-" + to_string(j), VTK_FLOAT, tmpVal, bId, numGrid, 3);
                bId += numGrid;
            }
        }
        if (bgp.CSFLAG) {
            out4vtk.OutputCELL_DATA_SCALARS(dest, "CSFLAG", VTK_INT, tmpVal, bId, numGrid, 0);
            bId += numGrid;
        }
        if (bgp.ITERNRDDM) {
            out4vtk.OutputCELL_DATA_SCALARS(dest, "ITERNRDDM", VTK_INT, tmpVal, bId, numGrid, 0);
            bId += numGrid;
        }
        if (bgp.ITERLSDDM) {
            out4vtk.OutputCELL_DATA_SCALARS(dest, "ITERLSDDM", VTK_INT, tmpVal, bId, numGrid, 0);
            bId += numGrid;
        }
        if (bgp.TIMELSDDM) {
            out4vtk.OutputCELL_DATA_SCALARS(dest, "TIMELSDDM", VTK_FLOAT, tmpVal, bId, numGrid, 0);
            bId += numGrid;
        }

        dest.close();
    }

    inV.close();
    
    if (remove(srcFile.c_str()) != 0) {
        OCP_WARNING("Failed to delete " + srcFile);
    }
    if (remove(myFile.c_str()) != 0) {
        OCP_WARNING("Failed to delete " + myFile);
    }

    free(timeInfo);
}


void OCPOutput::Setup(const ParamOutput& paramOutput, const OCPControl& ctrl, const Reservoir& rs)
{
    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Input Output Params -- begin");
    }

    workDir  = ctrl.GetWorkDir();
    fileName = ctrl.GetOCPFile();
    SetupComm(rs.GetDomain());

    summary.Setup(paramOutput.summary, rs);
    out4VTK.Setup(paramOutput.outVTKParam, workDir, rs);
    crtInfo.Setup();
    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Input Output Params -- end");
    }
}


void OCPOutput::SetupComm(const Domain& domain)
{
    myComm  = domain.global_comm;
    numproc = domain.global_numproc;
    myrank  = domain.global_rank;
}


void OCPOutput::SetValAtTimeStep(const Reservoir& rs, const OCPControl& ctrl, const OCPNRsuite& NR, GetWallTime& timer_total)
{
    GetWallTime timer;
    timer.Start();

    iters.Update(NR);
    summary.SetVal(rs, ctrl, iters, timer_total);
    crtInfo.SetVal(ctrl, NR);

    OCPTIME_OUTPUT += timer.Stop();
}

void OCPOutput::PrintAtTimeStep() const
{
    GetWallTime timer;
    timer.Start();

    summary.PrintInfo(workDir, fileName, (numproc > 1 ? myrank : -1));
    crtInfo.PrintFastReview(workDir, fileName, (numproc > 1 ? myrank : -1), iters);

    OCPTIME_OUTPUT += timer.Stop();
}

void OCPOutput::PrintInfoSched(const Reservoir&  rs,
                               const OCPControl& ctrl,
                               const OCP_DBL&    runtime) const
{
    OCP_DBL time = ctrl.time.GetCurrentTime();

    // print timing info on the screen
    if (ctrl.printLevel >= PRINT_MIN && myrank == MASTER_PROCESS) {
        cout << "Timestep " << setw(6) << left << iters.GetNumTimeStep() << ": " << fixed
             << setw(10) << setprecision(3) << right << time << " " << TIMEUNIT
             << "    Wall time: " << runtime << " Sec" << endl;
    }

    // print to output file
    GetWallTime timer;
    timer.Start();
    //out4RPT.PrintRPT(workDir, rs, days);
    out4VTK.PrintVTK(rs, ctrl);
    OCPTIME_OUTPUT += timer.Stop();
}


// By file
//void OCPOutput::PostProcess() const
//{
//    MPI_Barrier(myComm);
//    GetWallTime timer;
//    timer.Start();
//    if (numproc > 1 && myrank == MASTER_PROCESS) {
//        summary.PostProcess(workDir, fileName, numproc);
//        crtInfo.PostProcess(workDir, fileName, numproc, iters);           
//    }
//    if (myrank == MASTER_PROCESS) {
//        out4VTK.PostProcess(workDir, fileName, numproc);
//    }
//    
//    
//    OCPTIME_OUTPUT += timer.Stop();
//}


// By communication
void OCPOutput::PostProcess() const
{
    GetWallTime timer;
    timer.Start();
    file_myComm = myComm;
    file_myrank = myrank;
    summary.PostProcess(workDir, fileName, numproc);
    crtInfo.PostProcess(workDir, fileName, numproc, iters);
    if (myrank == MASTER_PROCESS) {
        out4VTK.PostProcess(workDir, fileName, numproc);
    }


    OCPTIME_OUTPUT += timer.Stop();
}


void OCPOutput::PrintCurrentTimeIter(const OCPControl& ctrl) const
{
    if (ctrl.printLevel >= PRINT_SOME && CURRENT_RANK == MASTER_PROCESS) {
        cout << "### DEBUG: " << setprecision(3) << scientific << ctrl.time.GetCurrentTime()
            << " " << TIMEUNIT;
        cout << ",  NR: " << iters.GetNRt() << ",  LS: " << iters.GetLSt()
            << ",  Last dt: " << ctrl.time.GetLastDt() << " " << TIMEUNIT << endl;
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/09/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/