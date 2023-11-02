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


void ItersInfo::Update(const OCPNRsuite& NRs)
{
    numTstep++;
    NRt  += NRs.GetIterNR();
    NRwt += NRs.GetIterNRw();
    LSt  += NRs.GetIterLS();
    LSwt += NRs.GetIterLSw();
}


void Summary::InputParam(const OutputSummary& summary_param)
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

    // cout << "Summary::InputParam" << endl;
}

void Summary::Setup(const Reservoir& rs)
{
    const USI maxRowNum = 1000;

    Sumdata.push_back(SumItem("TIME", "-", "DAY", "fixed", maxRowNum));
    Sumdata.push_back(SumItem("NRiter", "-", "-", "int", maxRowNum));
    Sumdata.push_back(SumItem("LSiter", "-", "-", "int", maxRowNum));
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

void Summary::SetVal(const Reservoir& rs, const OCPControl& ctrl, const ItersInfo& iters)
{
    const Bulk&     bulk  = rs.bulk;
    const AllWells& wells = rs.allWells;

    USI n = 0;

    // TIME
    Sumdata[n++].val.push_back(ctrl.time.GetCurrentTime());
    // NRiter
    Sumdata[n++].val.push_back(iters.GetNRt());
    // LSiter
    Sumdata[n++].val.push_back(iters.GetLSt());

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


void Summary::PostProcess(const string& dir, const string& filename, const OCP_INT& numproc) const
{
    vector<SumItem>* sumdata = const_cast<vector<SumItem>*>(&Sumdata);
    sumdata->clear();

    OCP_INT         bId = 0;
    USI             varlen = 0;
    OCP_USI         rowNum;

    vector<SumItem> mySum;
    vector<string>  buffer;
    OCP_BOOL        flag = OCP_FALSE;

    for (USI p = 0; p < numproc; p++) {
        const string myFile = dir + "proc_" + to_string(p) + "_SUMMARY.out";
        ifstream ifs(myFile, ios::in);

        if (!ifs) {
            OCP_MESSAGE("Trying to open file: " << (myFile));
            OCP_ABORT("Failed to open the input file!");
        }

        // Get time steps
        ReadLine(ifs, buffer, OCP_FALSE);
        for (USI i = 0; i < buffer.size(); i++) {
            if (buffer[i] == "--") {
                rowNum = stoi(buffer[i + 1]);
                break;
            }
        }

        while (!ifs.eof()) {
            if (buffer[0] == "TIME") {
                // begin to input title of data
                flag = OCP_TRUE;
                bId += varlen;
                varlen = buffer.size();
                for (const auto& v : buffer) {
                    mySum.push_back(SumItem(v, rowNum));
                }
                ReadLine(ifs, buffer, OCP_FALSE);
                OCP_ASSERT(varlen == buffer.size(), "Mismatch Value");
                for (USI i = 0; i < varlen; i++) {
                    mySum[bId + i].Unit = buffer[i];
                }
                ReadLine(ifs, buffer, OCP_FALSE);
                OCP_ASSERT(varlen == buffer.size(), "Mismatch Value");
                for (USI i = 0; i < varlen; i++) {
                    mySum[bId + i].Obj = buffer[i];
                }
            }
            ReadLine(ifs, buffer, OCP_FALSE);
            if (flag && buffer.size() == varlen && buffer[0] != "Row") {
                // begin to input value of data               
                for (USI i = 0; i < varlen; i++) {
                    mySum[bId + i].val.push_back(stod(buffer[i]));
                }
            }
            else {
                flag = OCP_FALSE;
            }
        }

        ifs.close();
        if (remove(myFile.c_str()) != 0) {
            OCP_WARNING("Failed to delete " + myFile);
        }
    }

    // post process
    // combine some values
    vector<SumItem> tmpdata;
    for (const auto& s : mySum) {
        USI n = 0;
        for (; n < tmpdata.size(); n++) {
            if (tmpdata[n] == s) {
                if (s.Item == "FOPR" || s.Item == "FGPR" || s.Item == "FWPR" ||
                    s.Item == "FOPT" || s.Item == "FGPT" || s.Item == "FWPT" ||
                    s.Item == "FGIR" || s.Item == "FGIT" || s.Item == "FWIR" ||
                    s.Item == "FWIT") {
                    for (USI i = 0; i < rowNum; i++) {
                        sumdata->at(n).val[i] += s.val[i];
                    }
                }
                else if (s.Item == "FPR" || s.Item == "FTR") {
                    for (USI i = 0; i < rowNum; i++) {
                        sumdata->at(n).val[i] += s.val[i] * (&s + 1)->val[i];
                        sumdata->at(n + 1).val[i] += (&s + 1)->val[i];
                    }
                }
                break;
            }
        }
        if (n == sumdata->size()) {
            sumdata->push_back(s);
            tmpdata.push_back(SumItem(s.Item, s.Obj));
            // prepare for calculating for FPR or FPT
            if (s.Item == "Volume") {
                for (USI i = 0; i < rowNum; i++) {
                    sumdata->at(n - 1).val[i] *= sumdata->at(n).val[i];
                }
            }
        }
    }
    // set output precision
    for (USI n = 0; n < sumdata->size(); n++) {
        if (sumdata->at(n).Item == "FPR" || sumdata->at(n).Item == "FTR") {
            for (USI i = 0; i < rowNum; i++) {
                sumdata->at(n).val[i] /= sumdata->at(n + 1).val[i];
            }
            sumdata->at(n).Type = "float";
            continue;
        }
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
void CriticalInfo::PostProcess(const string& dir, const string& filename, const OCP_INT& numproc, const ItersInfo& iters) const
{

    vector<SumItem>* sumdata = const_cast<vector<SumItem>*>(&Sumdata);
    sumdata->clear();

    OCP_USI         rowNum = 0;
    vector<string>  buffer;

    for (USI p = 0; p < numproc; p++) {
        const string myFile = dir + "proc_" + to_string(p) + "_FastReview.out";
        ifstream ifs(myFile, ios::in);

        if (!ifs) {
            OCP_MESSAGE("Trying to open file: " << (myFile));
            OCP_ABORT("Failed to open the input file!");
        }

        if (p == 0) {
            // Get time steps
            ReadLine(ifs, buffer, OCP_FALSE);
            for (USI i = 0; i < buffer.size(); i++) {
                if (buffer[i] == "--") {
                    rowNum += stoi(buffer[i + 1]);
                    continue;
                }
            }
            // Get Item
            ReadLine(ifs, buffer, OCP_FALSE);
            for (USI i = 0; i < buffer.size(); i++) {
                sumdata->push_back(SumItem(buffer[i], rowNum));
            }
            // Get Unit
            ReadLine(ifs, buffer, OCP_FALSE);
            for (USI i = 0; i < buffer.size(); i++) {
                sumdata->at(i).Unit = buffer[i];
            }
            // Get val
            for (OCP_USI n = 0; n < rowNum; n++) {
                ReadLine(ifs, buffer, OCP_FALSE);
                for (USI i = 0; i < buffer.size(); i++) {
                    sumdata->at(i).val.push_back(stod(buffer[i]));
                }
            }
        }
        else {
            // skip first three lines
            ReadLine(ifs, buffer, OCP_FALSE);
            ReadLine(ifs, buffer, OCP_FALSE);
            ReadLine(ifs, buffer, OCP_FALSE);
            // Get val
            for (OCP_USI n = 0; n < rowNum; n++) {
                ReadLine(ifs, buffer, OCP_FALSE);
                for (USI i = 2; i < buffer.size(); i++) {
                    sumdata->at(i).val[n] = max(sumdata->at(i).val[n], static_cast<OCP_DBL>(stod(buffer[i])));
                }
            }
        }
        ifs.close();
        if (remove(myFile.c_str()) != 0) {
            OCP_WARNING("Failed to delete " + myFile);
        }
    }
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


void BasicGridProperty::SetBasicGridProperty(const BasicGridPropertyParam& param)
{
    PRE  = param.PRE;
    PGAS = param.PGAS;
    PWAT = param.PWAT;
    SOIL = param.SOIL;
    SGAS = param.SGAS;
    SWAT = param.SWAT;
    DENO = param.DENO;
    DENG = param.DENG;
    DENW = param.DENW;
    KRO  = param.KRO;
    KRG  = param.KRG;
    KRW  = param.KRW;
    BOIL = param.BOIL;
    BGAS = param.BGAS;
    BWAT = param.BWAT;
    VOIL = param.VOIL;
    VGAS = param.VGAS;
    VWAT = param.VWAT;
    XMF  = param.XMF;
    YMF  = param.YMF;
    PCW  = param.PCW;
}


void BasicGridProperty::Check(const Reservoir& rs)
{
    // correct wrong output request
    if (!rs.IfOilExist()) {
        SOIL = OCP_FALSE;
        DENO = OCP_FALSE;
        BOIL = OCP_FALSE;
        KRO  = OCP_FALSE;
        VOIL = OCP_FALSE;
        XMF  = OCP_FALSE;
    }
    if (!rs.IfGasExist()) {
        SGAS = OCP_FALSE;
        DENG = OCP_FALSE;
        BGAS = OCP_FALSE;
        KRG  = OCP_FALSE;
        VGAS = OCP_FALSE;
        YMF  = OCP_FALSE;
    }
    if (!rs.IfWatExist()) {
        SWAT = OCP_FALSE;
        DENW = OCP_FALSE;
        BWAT = OCP_FALSE;
        KRW  = OCP_FALSE;
        VWAT = OCP_FALSE;
    }

    bgpnum = 0;
    if (PRE)      bgpnum++;
    if (PGAS)     bgpnum++;
    if (PWAT)     bgpnum++;
    if (SOIL)     bgpnum++;
    if (SGAS)     bgpnum++;
    if (SWAT)     bgpnum++;
    if (DENO)     bgpnum++;
    if (DENG)     bgpnum++;
    if (DENW)     bgpnum++;
    if (KRO)      bgpnum++;
    if (KRG)      bgpnum++;
    if (KRW)      bgpnum++;
    if (BOIL)     bgpnum++;
    if (BGAS)     bgpnum++;
    if (BWAT)     bgpnum++;
    if (VOIL)     bgpnum++;
    if (VGAS)     bgpnum++;
    if (VWAT)     bgpnum++;
    if (XMF)      bgpnum++;
    if (YMF)      bgpnum++;
    if (PCW)      bgpnum++;
}


void Out4RPT::InputParam(const OutputRPTParam& RPTparam)
{
    useRPT = RPTparam.useRPT;
    if (!useRPT) return;

    bgp.SetBasicGridProperty(RPTparam.bgp);
}

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

void Out4RPT::GetIJKGrid(USI& i, USI& j, USI& k, const OCP_ULL& n) const
{
    // i,j,k begin from 1
    // n must be the index of grids instead bulks
    k = n / (nx * ny) + 1;
    j = (n - (k - 1) * nx * ny) / nx + 1;
    i = n - (k - 1) * nx * ny - (j - 1) * nx + 1;
}

void Out4VTK::InputParam(const OutputVTKParam& VTKParam)
{
    useVTK = VTKParam.useVTK;
    if (!useVTK) return;

    bgp.SetBasicGridProperty(VTKParam.bgp);
}

void Out4VTK::Setup(const string& dir, const Reservoir& rs)
{
    if (!useVTK) return;

    bgp.Check(rs);

    // output the gloabl index of grids belonging to current domain
    const Domain& doman = rs.domain;
    if (doman.numproc > 1) {
        myFile = dir + "proc" + to_string(doman.myrank) + "_vtktmp.out";
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
        if (bgp.bgpnum == 0) {
            if (remove(myFile.c_str()) != 0) {
                OCP_WARNING("Fail to delete " + myFile);
            }
        }
    }
}

void Out4VTK::PrintVTK(const Reservoir& rs) const
{
    if (!useVTK)         return;
    if (bgp.bgpnum == 0) return;

    const BulkVarSet& bcv = rs.bulk.vs;

    const auto nb     = bcv.nbI;
    const auto np     = bcv.np;
    const auto OIndex = bcv.o;
    const auto GIndex = bcv.g;
    const auto WIndex = bcv.w;

    ofstream outF(myFile, ios::app | ios::binary);
    vector<OCP_DBL> tmpV(nb);
    // output    
    if (bgp.PRE)
        outF.write((const OCP_CHAR*)&bcv.P[0], nb * sizeof(bcv.P[0]));
    if (bgp.SOIL) {
        for (OCP_USI n = 0; n < nb; n++)
            tmpV[n] = bcv.S[n * np + OIndex];
        outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
    }
    if (bgp.SGAS) {
        for (OCP_USI n = 0; n < nb; n++)
            tmpV[n] = bcv.S[n * np + GIndex];
        outF.write((const OCP_CHAR*)&tmpV[0], nb * sizeof(tmpV[0]));
    }
    if (bgp.SWAT) {
        for (OCP_USI n = 0; n < nb; n++)
            tmpV[n] = bcv.S[n * np + WIndex];
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
    numGrid = out4vtk.InitASCII(dir, srcFile, "RUN of " + dir + filename);

    // Input cell values    
    vector<vector<OCP_DBL>> gridVal;        // each row is on a node of TSTEP in turn
    vector<OCP_ULL>         global_index;
    vector<OCP_USI>         mypart(numGrid);
    OCP_DBL*                workPtr;
    vector<OCP_DBL>         tmpVal;
    

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

        if (bgp.bgpnum > 0) {
            tmpVal.resize(numGridLoc * bgp.bgpnum);
            USI index = 0;
            while (OCP_TRUE) {
                if (index >= gridVal.size()) {
                    gridVal.push_back(vector<OCP_DBL>(bgp.bgpnum * numGrid));

                }
                workPtr = &gridVal[index][0];
                inV.read((OCP_CHAR*)(&tmpVal[0]), sizeof(tmpVal[0]) * tmpVal.size());
                const OCP_DBL* tamVap_ptr = &tmpVal[0];

                if (bgp.PRE) {
                    for (OCP_USI n = 0; n < numGridLoc; n++) {
                        workPtr[global_index[n]] = tamVap_ptr[n];
                    }
                    workPtr    += numGrid;
                    tamVap_ptr += numGridLoc;
                }
                if (bgp.SOIL) {
                    for (OCP_USI n = 0; n < numGridLoc; n++) {
                        workPtr[global_index[n]] = tamVap_ptr[n];
                    }
                    workPtr    += numGrid;
                    tamVap_ptr += numGridLoc;
                }
                if (bgp.SGAS) {
                    for (OCP_USI n = 0; n < numGridLoc; n++) {
                        workPtr[global_index[n]] = tamVap_ptr[n];
                    }
                    workPtr    += numGrid;
                    tamVap_ptr += numGridLoc;
                }
                if (bgp.SWAT) {
                    for (OCP_USI n = 0; n < numGridLoc; n++) {
                        workPtr[global_index[n]] = tamVap_ptr[n];
                    }
                    workPtr    += numGrid;
                    tamVap_ptr += numGridLoc;
                }

                index++;
                if (inV.eof()) {
                    break;
                }
            }
        }
        inV.close();
		if (remove(myfile.c_str()) != 0) {
			OCP_WARNING("Failed to delete " + myfile);
		}
    }

    // OutPut
    USI numTstep = gridVal.size();
    if (numTstep == 0) {
        ofstream source(srcFile, ios::app);
        source << "\n" << VTK_CELL_DATA << " " << numGrid;
        // Ouput partition
        out4vtk.OutputCELL_DATA_SCALARS(source, "PARTITION", VTK_UNSIGNED_INT, mypart, 0, numGrid, 0);
        source.close();
    }
    else {
        for (USI i = 0; i < numTstep; i++) {
            const string dstFile = dir + "TSTEP" + to_string(i) + ".vtk";
            ifstream source(srcFile, ios::binary);
            ofstream dest(dstFile, ios::binary);
            dest << source.rdbuf();
            source.close();

            dest << "\n" << VTK_CELL_DATA << " " << numGrid;
            // Ouput partition
            out4vtk.OutputCELL_DATA_SCALARS(dest, "PARTITION", VTK_UNSIGNED_INT, mypart, 0, numGrid, 0);

            OCP_ULL bId = 0;
            if (bgp.PRE) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "PRESSURE", VTK_FLOAT, gridVal[i], bId, numGrid, 3);
                bId += numGrid;
            }
            if (bgp.SOIL) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "SOIL", VTK_FLOAT, gridVal[i], bId, numGrid, 3);
                bId += numGrid;
            }
            if (bgp.SGAS) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "SGAS", VTK_FLOAT, gridVal[i], bId, numGrid, 3);
                bId += numGrid;
            }
            if (bgp.SWAT) {
                out4vtk.OutputCELL_DATA_SCALARS(dest, "SWAT", VTK_FLOAT, gridVal[i], bId, numGrid, 3);
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
    numGrid = out4vtk.InitASCII(dir, srcFile, "RUN of " + dir + filename);
 
    // Input cell values
    if (bgp.bgpnum == 0) return;
    vector<OCP_DBL> tmpVal(bgp.bgpnum * numGrid);
    ifstream inV(myFile, ios::in | ios::binary);
    if (!inV.is_open()) {
        OCP_WARNING("Can not open " + myFile);
    }   
    USI index = 0;
    while (OCP_TRUE) {

        inV.read((OCP_CHAR*)(&tmpVal[0]), sizeof(tmpVal[0]) * tmpVal.size());

        const string dstFile = dir + "TSTEP" + to_string(index++) + ".vtk";
        ifstream source(srcFile, ios::binary);
        ofstream dest(dstFile, ios::binary);
        dest << source.rdbuf();
        source.close();

        dest << "\n" << VTK_CELL_DATA << " " << numGrid;
        OCP_ULL bId = 0;
        if (bgp.PRE) {
            out4vtk.OutputCELL_DATA_SCALARS(dest, "PRESSURE", VTK_FLOAT, tmpVal, bId, numGrid, 3);
            bId += numGrid;
        }
        if (bgp.SOIL) {
            out4vtk.OutputCELL_DATA_SCALARS(dest, "SOIL", VTK_FLOAT, tmpVal, bId, numGrid, 3);
            bId += numGrid;
        }
        if (bgp.SGAS) {
            out4vtk.OutputCELL_DATA_SCALARS(dest, "SGAS", VTK_FLOAT, tmpVal, bId, numGrid, 3);
            bId += numGrid;
        }
        if (bgp.SWAT) {
            out4vtk.OutputCELL_DATA_SCALARS(dest, "SWAT", VTK_FLOAT, tmpVal, bId, numGrid, 3);
            bId += numGrid;
        }

        dest.close();
        if (inV.eof()) {
            inV.close();
            break;
        }
    }
    
    if (remove(srcFile.c_str()) != 0) {
        OCP_WARNING("Failed to delete " + srcFile);
    }
    if (remove(myFile.c_str()) != 0) {
        OCP_WARNING("Failed to delete " + myFile);
    }
}


void OCPOutput::InputParam(const ParamOutput& paramOutput)
{
    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Input Output Params -- begin");
    }

    summary.InputParam(paramOutput.summary);
    out4RPT.InputParam(paramOutput.outRPTParam);
    out4VTK.InputParam(paramOutput.outVTKParam);

    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Input Output Params -- end");
    }
}

void OCPOutput::Setup(const Reservoir& reservoir, const OCPControl& ctrl, const Domain& domain)
{
    myComm  = domain.myComm;
    numproc = domain.numproc;
    myrank  = domain.myrank;

    workDir  = ctrl.GetWorkDir();
    fileName = ctrl.GetOCPFile();
    summary.Setup(reservoir);
    crtInfo.Setup();
    // out4RPT.Setup(workDir, reservoir);
    out4VTK.Setup(workDir, reservoir);
}

void OCPOutput::SetVal(const Reservoir& reservoir, const OCPControl& ctrl, const OCPNRsuite& NR)
{
    GetWallTime timer;
    timer.Start();

    iters.Update(NR);
    summary.SetVal(reservoir, ctrl, iters);
    crtInfo.SetVal(ctrl, NR);

    OCPTIME_OUTPUT += timer.Stop() / TIME_S2MS;
}

void OCPOutput::PrintInfo() const
{
    GetWallTime timer;
    timer.Start();

    summary.PrintInfo(workDir, fileName, (numproc > 1 ? myrank : -1));
    crtInfo.PrintFastReview(workDir, fileName, (numproc > 1 ? myrank : -1), iters);

    OCPTIME_OUTPUT += timer.Stop() / TIME_S2MS;
}

void OCPOutput::PrintInfoSched(const Reservoir&  rs,
                               const OCPControl& ctrl,
                               const OCP_DBL&    time) const
{
    OCP_DBL days = ctrl.time.GetCurrentTime();

    // print timing info on the screen
    if (ctrl.printLevel >= PRINT_MIN && myrank == MASTER_PROCESS) {
        cout << "Timestep " << setw(6) << left << iters.GetNumTimeStep() << ": " << fixed
             << setw(10) << setprecision(3) << right << days << TIMEUNIT
             << "    Wall time: " << time / TIME_S2MS << " Sec" << endl;
    }

    // print to output file
    GetWallTime timer;
    timer.Start();
    //out4RPT.PrintRPT(workDir, rs, days);
    out4VTK.PrintVTK(rs);
    OCPTIME_OUTPUT += timer.Stop() / TIME_S2MS;
}


void OCPOutput::PostProcess() const
{
    MPI_Barrier(myComm);
    GetWallTime timer;
    timer.Start();
    if (numproc > 1 && myrank == MASTER_PROCESS) {
        summary.PostProcess(workDir, fileName, numproc);
        crtInfo.PostProcess(workDir, fileName, numproc, iters);           
    }
    if (myrank == MASTER_PROCESS) {
        out4VTK.PostProcess(workDir, fileName, numproc);
    }
    
    
    OCPTIME_OUTPUT += timer.Stop() / TIME_S2MS;
}


void OCPOutput::PrintCurrentTimeIter(const OCPControl& ctrl) const
{
    if (ctrl.printLevel >= PRINT_SOME && CURRENT_RANK == MASTER_PROCESS) {
        cout << "### DEBUG: " << setprecision(3) << scientific << ctrl.time.GetCurrentTime()
            << TIMEUNIT;
        cout << ",  NR: " << iters.GetNRt() << ",  LS: " << iters.GetLSt()
            << ",  Last dt: " << ctrl.time.GetLastDt() << TIMEUNIT << endl;
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