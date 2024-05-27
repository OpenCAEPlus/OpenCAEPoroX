/*! \file    ParamRead.cpp
 *  \brief   ParamRead class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>

#include "ParamRead.hpp"


void ParamRead::Print(std::ostream &out)
{
    char indent = ' ';

    //
    out << "unitType: " << paramRs.unitType << '\n';

    //
    out << "blockOil? " << paramRs.blackOil << '\n';

    // InputCOMPS
    out << "\nInputCOMPS\n";
    out << "comps: " << paramRs.comps << '\n'
        << "numCom: " << paramRs.numCom << '\n'
        << "comsParam.numCom: " << paramRs.comsParam.numCom << '\n';

    out << "\nLBCcoef\n";
    for (auto val: paramRs.comsParam.LBCcoef)
        out << val << indent;

    out << "\nSSMparamSTA\n";
    for (auto val: paramRs.comsParam.SSMparamSTA)
        out << val << indent;

    out << "\nNRparamSTA\n";
    for (auto val: paramRs.comsParam.NRparamSTA)
        out << val << indent;

    out << "\nSSMparamSP\n";
    for (auto val: paramRs.comsParam.SSMparamSP)
        out << val << indent;

    out << "\nNRparamSP\n";
    for (auto val: paramRs.comsParam.NRparamSP)
        out << val << indent;

    out << "\nRRparam\n";
    for (auto val: paramRs.comsParam.RRparam)
        out << val << indent;

    // THERMAL
    out << "\nTHERMAL\n";
    out << paramRs.thermal << indent
        << paramWell.thermal << indent
        << indent;
    //
    out << "\noil, gas, water, disGas\n";
    out << paramRs.oil << indent
        << paramRs.gas <<indent
        << paramRs.water << indent
        << paramRs.disGas << indent
        << paramRs.GRAVDR << indent
        << paramRs.initType << indent
        << '\n';

    // RTEMP
    out << paramRs.rsTemp << '\n';

    // InputTABLE
    out << "\nInputTABLE\n";
    out << "SWFN_T: " << paramRs.SWFN_T.name << '\n';
    for (auto first: paramRs.SWFN_T.data)
    {
        for (auto second: first)
        {
            for (auto third: second)
                out << third << indent;
            out << '\n';
        }
        out << '\n';
    }
    //
    out << "SWOF_T: " << paramRs.SWOF_T.name << '\n';
    for (auto first: paramRs.SWOF_T.data)
    {
        for (auto second: first)
        {
            for (auto third: second)
                out << third << indent;
            out << '\n';
        }
        out << '\n';
    }
    //
    out << "SGFN_T: " << paramRs.SGFN_T.name << '\n';
    for (auto first: paramRs.SGFN_T.data)
    {
        for (auto second: first)
        {
            for (auto third: second)
                out << third << indent;
            out << '\n';
        }
        out << '\n';
    }
    //
    out << "SGOF_T " << paramRs.SGOF_T.name << '\n';
    for (auto first: paramRs.SGOF_T.data)
    {
        for (auto second: first)
        {
            for (auto third: second)
                out << third << indent;
            out << '\n';
        }
        out << '\n';
    }
    //
    out << "SOF3_T " << paramRs.SOF3_T.name << '\n';
    for (auto first: paramRs.SOF3_T.data)
    {
        for (auto second: first)
        {
            for (auto third: second)
                out << third << indent;
            out << '\n';
        }
        out << '\n';
    }
    //
    out << "PBVD_T " << paramRs.PBVD_T.name << '\n';
    for (auto first: paramRs.PBVD_T.data)
    {
        for (auto second: first)
        {
            for (auto third: second)
                out << third << indent;
            out << '\n';
        }
        out << '\n';
    }
    //
    out << "PVCO_T " << paramRs.PVCO_T.name << '\n';
    for (auto first: paramRs.PVCO_T.data)
    {
        for (auto second: first)
        {
            for (auto third: second)
                out << third << indent;
            out << '\n';
        }
        out << '\n';
    }
    //
    out << paramRs.PVDO_T.name << '\n';
    for (auto first: paramRs.PVDO_T.data)
    {
        for (auto second: first)
        {
            for (auto third: second)
                out << third << indent;
            out << '\n';
        }
        out << '\n';
    }
    //
    out << paramRs.PVCDO_T.name << '\n';
    for (auto first: paramRs.PVCDO_T.data)
    {
        for (auto second: first)
        {
            for (auto third: second)
                out << third << indent;
            out << '\n';
        }
        out << '\n';
    }
    //
    out << paramRs.PVDG_T.name << '\n';
    for (auto first: paramRs.PVDG_T.data)
    {
        for (auto second: first)
        {
            for (auto third: second)
                out << third << indent;
            out << '\n';
        }
        out << '\n';
    }
    //
    out << paramRs.PVTW_T.name << '\n';
    for (auto first: paramRs.PVTW_T.data)
    {
        for (auto second: first)
        {
            for (auto third: second)
                out << third << indent;
            out << '\n';
        }
        out << '\n';
    }
    //
    out << paramRs.ZMFVD_T.name << '\n';
    for (auto first: paramRs.ZMFVD_T.data)
    {
        for (auto second: first)
        {
            for (auto third: second)
                out << third << indent;
            out << '\n';
        }
        out << '\n';
    }
    //
    out << paramRs.TEMPVD_T.name << '\n';
    for (auto first: paramRs.TEMPVD_T.data)
    {
        for (auto second: first)
        {
            for (auto third: second)
                out << third << indent;
            out << '\n';
        }
        out << '\n';
    }

    // InputTABLE2
    out << "\nInputTABLE2\n";
    out << paramRs.comsParam.viscTab.name << '\n';
    for (auto first: paramRs.comsParam.viscTab.data)
    {
        for (auto second: first.data)
        {
            for (auto third: second)
            {
                for (auto forth: third)
                    out << forth << indent;
                out << '\n';
            }
            out << '\n';
        }
        out << '\n';
    }
    //
    out << paramRs.PVTH2O.name << '\n';
    for (auto first: paramRs.PVTH2O.data)
    {
        for (auto second: first.data)
        {
            for (auto third: second)
            {
                for (auto forth: third)
                    out << forth << indent;
                out << '\n';
            }
            out << '\n';
        }
        out << '\n';
    }
    //
    out << paramRs.PVTCO2.name << '\n';
    for (auto first: paramRs.PVTCO2.data)
    {
        for (auto second: first.data)
        {
            for (auto third: second)
            {
                for (auto forth: third)
                    out << forth << indent;
                out << '\n';
            }
            out << '\n';
        }
        out << '\n';
    }


    // GARCIAW...
    out << "GARCIAW\n";
    out << paramRs.GARCIAW << indent << "\n";

    // InputROCK
    out << "\nInputROCK";
    for (auto rock: paramRs.rockSet)
        out << "\ntype: " << rock.type << ", Pref: " << rock.Pref << ", cp1: " << rock.cp1 << ", cp2: " << rock.cp2 << indent;

    // InputROCKT
    out << "\nInputROCKT";
    for (auto rock: paramRs.rockSet)
        out << "\ntype: " << rock.type << ", Pref: " << rock.Pref << ", Tref: "
            << rock.Tref << ", cp1: " << rock.cp1 << ' ' << rock.ct << ' ' << rock.cpt << ' '
            << rock.HCP1 << ' ' << rock.HCP2 << indent;

    // InputBrooksCorey
    out << "\nInputBrooksCorey\n";
    for (auto bc: paramRs.BCparam)
    {
        out << bc.sw_imm << ' ' << bc.sn_imm << ' ' << bc.Pentry << ' ' << bc.Pcmax << ' ' << bc.Cw_kr << ' '
            << bc.Cn_kr << ' ' << bc.C_pc << indent;
    }

    // InputHLOSS
    out << "\nInputHLOSS\n";
    out << paramRs.hLoss.ifHLoss << ' '
        << paramRs.hLoss.obUse << ' '
        << paramRs.hLoss.obC << ' '
        << paramRs.hLoss.obK << ' '
        << indent;

    // InputMISCSTR
    out << "\nInputMISCSTR\n";
    out << paramRs.miscstr.ifMiscible << ' '
        << paramRs.miscstr.surTenRef << ' '
        << paramRs.miscstr.surTenEpt << ' '
        << paramRs.miscstr.surTenPc << ' '
        << paramRs.miscstr.surTenExp << indent;

    // InputGRAVITY
    out << "\nInputGRAVITY\n";
    out << paramRs.gravity.activity << ' '
        << paramRs.gravity.data[0] << ' ' << paramRs.gravity.data[1] << ' ' << paramRs.gravity.data[2] << indent;

    // InputDENSITY
    out << "\nInputDENSITY\n";
    out << paramRs.density.activity << ' '
        << paramRs.density.data[0] << ' ' << paramRs.density.data[1] << ' ' << paramRs.density.data[2] << indent;

    // InputTHCON
    out << "\nInputTHCON\n";
    out << paramRs.ifThcon << ' '
        << paramRs.thcono << ' ' << paramRs.thcong << ' ' << paramRs.thconw << ' ' << paramRs.thconr << indent;

    // InputEQUIL
    out << "\nInputEQUIL\n";
    for (auto equil: paramRs.EQUIL)
    {
        for (auto val: equil.data)
            out << val << ' ';
    }

    // InputTABDIMS
    out << "InputTABDIMS\n";
    out << paramRs.NTSFUN << ' ' << paramRs.NTPVT << ' ' << paramRs.NTROOC << ' '
        << paramRs.comsParam.NTPVT << indent;

    // InputMETHOD
    out << "\nInputMETHOD\n";
    for (auto val: paramControl.method)
        out << val << indent;
    for (auto val: paramControl.lsFile)
        out << val << indent;

    // InputTUNING
    out << "\nInputTUNING\n";
    for (auto pair: paramControl.tuning_T)
    {
        out << pair.d << indent;
        for (auto first: pair.Tuning)
            for (auto second: first)
                out << second << indent;
        out << '\n';
    }

    // InputWELSPECS, InputCOMPDAT
    out << "\nInputWELSPECS, InputCOMPDAT\n";
    for (auto wl: paramWell.well)
    {
        if (wl.gridType == GridType::unstructured)
        {
            out << "well grid type: unstructured\n" << wl.name << ' ' << wl.group << ' '
                << wl.X << ' ' << wl.Y << ' ' << wl.Z << indent;
            //
            out << "\nX_perf:\n";
            for (auto val: wl.X_perf)
                out << val << indent;
            out << "\nY_perf:\n";
            for (auto val: wl.Y_perf)
                out << val << indent;
            out << "\nZ_perf:\n";
            for (auto val: wl.Z_perf)
                out << val << indent;
            out << "\nWI:\n";
            for (auto val: wl.WI)
                out << val << indent;
            out << "\ndiameter:\n";
            for (auto val: wl.diameter)
                out << val << indent;
            out << "\nkh:\n";
            for (auto val: wl.kh)
                out << val << indent;
            out << "\nskinFactor:\n";
            for (auto val: wl.skinFactor)
                out << val << indent;
            out << "\ndirection:\n";
            for (auto val: wl.direction)
                out << val << indent;
            //
        }
        else if (wl.gridType == GridType::structured)
        {
            out << "well grid type: structured\nwell name: " << wl.name << ", group: " << wl.group << ", I: "
                << wl.I << ", J: " << wl.J << ", depth: " << wl.depth << indent;
            //
            out << "\nI_perf:\n";
            for (auto val: wl.I_perf)
                out << val << indent;
            out << "\nJ_perf:\n";
            for (auto val: wl.J_perf)
                out << val << indent;
            out << "\nK_perf:\n";
            for (auto val: wl.K_perf)
                out << val << indent;
            out << "\nWI:\n";
            for (auto val: wl.WI)
                out << val << indent;
            out << "\ndiameter:\n";
            for (auto val: wl.diameter)
                out << val << indent;
            out << "\nkh:\n";
            for (auto val: wl.kh)
                out << val << indent;
            out << "\nskinFactor:\n";
            for (auto val: wl.skinFactor)
                out << val << indent;
            out << "\ndirection:\n";
            for (auto val: wl.direction)
                out << val << indent;
            //
        }
        else
            OCP_ABORT("Wrong well grid type!");

        out << '\n';
    }

    // InputWCONINJE, InputWCONPROD, InputUNWEIGHT
    out << "\nInputWCONINJE, InputWCONPROD, InputUNWEIGHT\n";
    for (auto wl: paramWell.well)
    {
        out << "ifUseUnweight: " << wl.ifUseUnweight << indent;
        for (auto opt: wl.optParam)
            out << "\nwell name: " << wl.name << ", well operations, time: " << opt.d
                << ", type: " << opt.opt.type << ' ' << opt.opt.fluidType << ", state: " << opt.opt.state << ", mode: "
                << opt.opt.mode << ", maxRate: " << opt.opt.maxRate << ", maxBHP: " << opt.opt.maxBHP << ", minBHP: "
                << opt.opt.minBHP << ", injTemp: " << opt.opt.injTemp
                << indent;

        out << '\n';
    }

    // InputTSTEP
    out << "\nInputTSTEP\n";
    out << "paramWell.criticalTime:\n";
    int i=0;
    for (auto val: paramWell.criticalTime)
    {
        out << "idx: " << i << ", " << val << '\n';
        ++i;
    }
    out << "paramControl.criticalTime:\n";
    for (auto val: paramControl.criticalTime)
        out << val << ' ';

    // InputWELTARG
    // InputWTEMP

    // InputWELLSTRE
    out << "\nInputWELLSTRE\n";
    for (auto sol: paramWell.solSet)
    {
        out << "\nsolSet name: " << sol.name << indent;
        for (auto val: sol.comRatio)
            out << val << indent;
    }

    // InputPSURF
    out << "\nPsurf: " << paramWell.Psurf << ' ' << paramRs.Psurf << indent;

    // InputTSURF
    out << "\nTsurf: " << paramWell.Tsurf << ' ' << paramRs.Tsurf << indent;

    // InputWELINITP
    out << "\nInputWELINITP:\n";
    for (auto wl: paramWell.well)
        out << wl.name << ' ' << wl.initP << "\n";

    // InputSUMMARY fff
    out << "\nInputSUMMARY:\n";
    out << paramOutput.summary.FPR << ' '
        << paramOutput.summary.FTR << ' '
        << paramOutput.summary.FOPR << ' '
        << paramOutput.summary.FOPT << ' '
        << paramOutput.summary.FGPR << ' '
        << paramOutput.summary.FGPt << ' '
        << paramOutput.summary.FWPR << ' '
        << paramOutput.summary.FWPT << ' '
        << paramOutput.summary.FGIR << ' '
        << paramOutput.summary.FGIT << ' '
        << paramOutput.summary.FWIR << ' '
        << paramOutput.summary.FWIT << ' '
        << indent;

    // InputRPTSCHED fff

    // InputNCOMPS
    out << paramRs.comsParam.numCom << ' ' << paramRs.numCom << indent;

    // InputCNAMES
    out << "\nInputCNAMES:\n";
    for (auto val: paramRs.comsParam.Cname)
        out << val << ' ';

    // InputCOMPONENTS
    out << "\nInputCOMPONENTS\n";
    out << paramRs.comsParam.Tc.activity << indent;
    for (auto comp: paramRs.comsParam.Tc.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.Pc.activity << indent;
    for (auto comp: paramRs.comsParam.Pc.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.Vc.activity << indent;
    for (auto comp: paramRs.comsParam.Vc.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.Zc.activity << indent;
    for (auto comp: paramRs.comsParam.Zc.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.MW.activity << indent;
    for (auto comp: paramRs.comsParam.MW.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.Acf.activity << indent;
    for (auto comp: paramRs.comsParam.Acf.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.OmegaA.activity << indent;
    for (auto comp: paramRs.comsParam.OmegaA.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.OmegaB.activity << indent;
    for (auto comp: paramRs.comsParam.OmegaB.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.Vshift.activity << indent;
    for (auto comp: paramRs.comsParam.Vshift.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.parachor.activity << indent;
    for (auto comp: paramRs.comsParam.parachor.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.Vcvis.activity << indent;
    for (auto comp: paramRs.comsParam.Vcvis.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.Zcvis.activity << indent;
    for (auto comp: paramRs.comsParam.Zcvis.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.molden.activity << indent;
    for (auto comp: paramRs.comsParam.molden.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.cp.activity << indent;
    for (auto comp: paramRs.comsParam.cp.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.ct1.activity << indent;
    for (auto comp: paramRs.comsParam.ct1.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.ct2.activity << indent;
    for (auto comp: paramRs.comsParam.ct2.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.cpt.activity << indent;
    for (auto comp: paramRs.comsParam.cpt.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.cpl1.activity << indent;
    for (auto comp: paramRs.comsParam.cpl1.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.cpl2.activity << indent;
    for (auto comp: paramRs.comsParam.cpl2.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.cpl3.activity << indent;
    for (auto comp: paramRs.comsParam.cpl3.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.cpl4.activity << indent;
    for (auto comp: paramRs.comsParam.cpl4.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.cpg1.activity << indent;
    for (auto comp: paramRs.comsParam.cpg1.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.cpg2.activity << indent;
    for (auto comp: paramRs.comsParam.cpg2.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.cpg3.activity << indent;
    for (auto comp: paramRs.comsParam.cpg3.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.cpg4.activity << indent;
    for (auto comp: paramRs.comsParam.cpg4.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.hvapr.activity << indent;
    for (auto comp: paramRs.comsParam.hvapr.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.hvr.activity << indent;
    for (auto comp: paramRs.comsParam.hvr.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.ev.activity << indent;
    for (auto comp: paramRs.comsParam.ev.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.avisc.activity << indent;
    for (auto comp: paramRs.comsParam.avisc.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.bvisc.activity << indent;
    for (auto comp: paramRs.comsParam.bvisc.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.avg.activity << indent;
    for (auto comp: paramRs.comsParam.avg.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }
    out << paramRs.comsParam.bvg.activity << indent;
    for (auto comp: paramRs.comsParam.bvg.data)
    {
        for (auto& v: comp)
            out << v << ' ';
        out << indent;
    }

    // InputRefPR
    for (auto val: paramRs.comsParam.Pref)
        out << val << ' ';
    for (auto val: paramRs.comsParam.Tref)
        out << val << ' ';

    // InputLBCCOEF
    cout << paramRs.comsParam.NTPVT << indent;
    for (auto val: paramRs.comsParam.LBCcoef)
        out << val << ' ';

    // InputBIC
    for (auto first: paramRs.comsParam.BIC)
    {
        for (auto val: first)
            out << val << ' ';
        out << indent;
    }

    // InputSSMSTA
    for (auto val: paramRs.comsParam.SSMparamSTA)
        out << val << ' ';
    // InputSSMSP
    for (auto val: paramRs.comsParam.SSMparamSP)
        out << val << ' ';
    // InputNRSTA
    for (auto val: paramRs.comsParam.NRparamSTA)
        out << val << ' ';
    // InputNRSP
    for (auto val: paramRs.comsParam.NRparamSP)
        out << val << ' ';
    // InputRR
    for (auto val: paramRs.comsParam.RRparam)
        out << val << ' ';

    // InputBoundary
    for (auto bd: paramRs.BDparam)
        out << bd.name << ' ' << bd.constP << ' ' << bd.P << indent;

    // InpuCurTime
    out << paramControl.curTime << indent;

    // InpuMaxSimTime
    out << paramControl.MaxSimTime << indent;
}


/// Initialize paramRs, paramWell, and paramControl.
void ParamRead::Init()
{
    paramRs.Init();
    paramWell.Init();
    paramControl.Init(workDir, fileName);

    // For HiSim input
    recurrent_type = 1;
}

/// Get workDir and fileName from inputFile.
void ParamRead::GetDirAndName()
{
#if defined(_CONSOLE) || defined(_WIN32) || defined(_WIN64)
    // for Window file system
    OCP_INT pos = inputFile.find_last_of('\\') + 1;
    workDir     = inputFile.substr(0, pos);
    fileName    = inputFile.substr(pos, inputFile.size() - pos);
#else
    // for Linux and Mac OSX file system
    OCP_INT pos = inputFile.find_last_of('/') + 1;
    workDir     = inputFile.substr(0, pos);
    fileName    = inputFile.substr(pos, inputFile.size() - pos);
#endif
}

/// This is the general interface for reading input files.
void ParamRead::ReadInputFile(const string& filename, int type)
{
    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Input Reservoir Param -- begin");
    }

    inputFile = filename;
    GetDirAndName();
    Init();
    if (type == 0)
        ReadFile(inputFile);
    else
    {
        ReadFileHiSim(inputFile);
        ReadRestFile(inputFile);
    }
    {
        ofstream inputfile("./inputfile2");
        Print(inputfile);
        inputfile.close();
    }
    CheckParam();
    SetUnit(paramRs.unitType);

    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Input Reservoir Param -- end");
    }
}

/// Read parameters from a file, which is called in ReadInputFile.
void ParamRead::ReadFile(const string& filename)
{
    ifstream ifs(filename, ios::in);
    if (!ifs) {
        OCP_MESSAGE("Trying to open file: " << (filename));
        OCP_ABORT("Failed to open the input file!");
    }

    while (!ifs.eof()) {
        vector<string> vbuf;
        if (!ReadLine(ifs, vbuf)) break;
        string keyword = vbuf[0];

        switch (Map_Str2Int(&keyword[0], keyword.size())) {

            case Map_Str2Int("FIELD", 5):
                if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT)
                    cout << "FIELD" << endl;
                paramRs.unitType = "FIELD";
                break;

            case Map_Str2Int("METRIC", 6):
                if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT)
                    cout << "METRIC" << endl;
                paramRs.unitType = "METRIC";
                break;

            case Map_Str2Int("SPE11A", 6):
                if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT)
                    cout << "SPE11A" << endl;
                paramRs.unitType = "SPE11A";
                break;

            case Map_Str2Int("SPE11Amg", 8):
                if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT)
                    cout << "SPE11Amg" << endl;
                paramRs.unitType = "SPE11Amg";
                break;

            case Map_Str2Int("BLACKOIL", 8):
                paramRs.blackOil = OCP_TRUE;
                break;

            case Map_Str2Int("COMPS", 5):
                paramRs.InputCOMPS(ifs);
                break;

            case Map_Str2Int("THERMAL", 7):
                paramRs.thermal    = OCP_TRUE;
                paramWell.thermal  = OCP_TRUE;
                paramControl.model = OCPModel::thermal;
                break;

            case Map_Str2Int("OIL", 3):
                paramRs.oil = OCP_TRUE;
                break;

            case Map_Str2Int("GAS", 3):
                paramRs.gas = OCP_TRUE;
                break;

            case Map_Str2Int("WATER", 5):
                paramRs.water = OCP_TRUE;
                break;

            case Map_Str2Int("DISGAS", 6):
                paramRs.disGas = OCP_TRUE;
                break;

            case Map_Str2Int("GRAVDR", 6):
                paramRs.GRAVDR = OCP_TRUE;
                break;

            case Map_Str2Int("INITPTN0", 8):
            case Map_Str2Int("INITPTN1", 8):
                paramRs.initType = "INITPTN";
                break;

            case Map_Str2Int("EQUILWAT", 8):
                paramRs.initType = "EQUILWAT";
                break;

            case Map_Str2Int("RTEMP", 5):
                paramRs.InputRTEMP(ifs);
                break;
                
            case Map_Str2Int("SWFN", 4):
            case Map_Str2Int("SWOF", 4):
            case Map_Str2Int("SGFN", 4):
            case Map_Str2Int("SGOF", 4):
            case Map_Str2Int("SOF3", 4):
            case Map_Str2Int("PVCO", 4):
            case Map_Str2Int("PVDO", 4):
            case Map_Str2Int("PVCDO", 5):
            case Map_Str2Int("PVDG", 4):
            case Map_Str2Int("PVTW", 4):
            case Map_Str2Int("PBVD", 4):
            case Map_Str2Int("ZMFVD", 5):
            case Map_Str2Int("TEMPVD", 6):
                paramRs.InputTABLE(ifs, keyword);
                break;

            case Map_Str2Int("VISCTAB", 7):
            case Map_Str2Int("PVTH2O", 6):
            case Map_Str2Int("PVTCO2", 6):
                paramRs.InputTABLE2(ifs, keyword);
                break;

            case Map_Str2Int("GARCIAW", 7):
                paramRs.GARCIAW = OCP_TRUE;
                break;

            case Map_Str2Int("ROCK", 4):
                paramRs.InputROCK(ifs);
                break;

            case Map_Str2Int("ROCKT", 5):
                paramRs.InputROCKT(ifs);
                break;

            case Map_Str2Int("BCPERM", 6):
                paramRs.InputBrooksCorey(ifs);
                break;

            case Map_Str2Int("HLOSS", 5):
                paramRs.InputHLOSS(ifs);
                break;

            case Map_Str2Int("MISCSTR", 7):
                paramRs.InputMISCSTR(ifs);
                break;

            case Map_Str2Int("GRAVITY", 7):
                paramRs.InputGRAVITY(ifs);
                break;

            case Map_Str2Int("DENSITY", 7):
                paramRs.InputDENSITY(ifs);
                break;

            case Map_Str2Int("THCONO", 6):
            case Map_Str2Int("THCONG", 6):
            case Map_Str2Int("THCONW", 6):
            case Map_Str2Int("THCONR", 6):
                paramRs.InputTHCON(ifs, keyword);
                break;

            case Map_Str2Int("EQUIL", 5):
                paramRs.InputEQUIL(ifs);
                break;

            case Map_Str2Int("TABDIMS", 7):
                paramRs.InputTABDIMS(ifs);
                break;

            case Map_Str2Int("INCLUDE", 7):
                ReadINCLUDE(ifs);
                break;

            case Map_Str2Int("METHOD", 6):
                paramControl.InputMETHOD(ifs);
                break;

            case Map_Str2Int("TUNING", 6):
                paramControl.InputTUNING(ifs);
                break;

            case Map_Str2Int("WELSPECS", 8):
                paramWell.InputWELSPECS(ifs);
                break;

            case Map_Str2Int("COMPDAT", 7):
                paramWell.InputCOMPDAT(ifs);
                break;

            case Map_Str2Int("WCONINJE", 8):
                paramWell.InputWCONINJE(ifs);
                break;

            case Map_Str2Int("WCONPROD", 8):
                paramWell.InputWCONPROD(ifs);
                break;

            case Map_Str2Int("UNWEIGHT", 8):
                paramWell.InputUNWEIGHT(ifs);
                break;

            case Map_Str2Int("TSTEP", 5):
                paramWell.InputTSTEP(ifs);
                paramControl.criticalTime = paramWell.criticalTime;
                break;

            case Map_Str2Int("WELTARG", 7):
            case Map_Str2Int("WELLTARG", 8):
                paramWell.InputWELTARG(ifs);
                break;

            case Map_Str2Int("WTEMP", 5):
                paramWell.InputWTEMP(ifs);
                break;

            case Map_Str2Int("WELLSTRE", 8):
                paramWell.InputWELLSTRE(ifs);
                break;

            case Map_Str2Int("PSURF", 5):
                paramWell.InputPSURF(ifs);
                paramRs.Psurf = paramWell.Psurf;
                break;

            case Map_Str2Int("TSURF", 5):
                paramWell.InputTSURF(ifs);
                paramRs.Tsurf = paramWell.Tsurf;
                break;

            case Map_Str2Int("WELINITP", 8):
                paramWell.InputWELINITP(ifs);
                break;

            case Map_Str2Int("SUMMARY", 7):
                paramOutput.InputSUMMARY(ifs);
                break;

            case Map_Str2Int("RPTSCHED", 8):
            case Map_Str2Int("VTKSCHED", 8):
                paramOutput.InputRPTSCHED(ifs, keyword);
                break;

            case Map_Str2Int("NCOMPS", 6):
                paramRs.InputNCOMPS(ifs);
                break;

            case Map_Str2Int("CNAMES", 6):
                paramRs.InputCNAMES(ifs);
                break;

            case Map_Str2Int("TCRIT", 5):
            case Map_Str2Int("PCRIT", 5):
            case Map_Str2Int("VCRIT", 5):
            case Map_Str2Int("ZCRIT", 5):
            case Map_Str2Int("MW", 2):
            case Map_Str2Int("ACF", 3):
            case Map_Str2Int("OMEGAA", 6):
            case Map_Str2Int("OMEGAB", 6):
            case Map_Str2Int("SSHIFT", 6):
            case Map_Str2Int("PARACHOR", 8):
            case Map_Str2Int("VCRITVIS", 8):
            case Map_Str2Int("MOLDEN", 6):
            case Map_Str2Int("CP", 2):
            case Map_Str2Int("CT1", 3):
            case Map_Str2Int("CT2", 3):
            case Map_Str2Int("CPT", 3):
            case Map_Str2Int("CPL1", 4):
            case Map_Str2Int("CPL2", 4):
            case Map_Str2Int("CPL3", 4):
            case Map_Str2Int("CPL4", 4):
            case Map_Str2Int("CPG1", 4):
            case Map_Str2Int("CPG2", 4):
            case Map_Str2Int("CPG3", 4):
            case Map_Str2Int("CPG4", 4):
            case Map_Str2Int("HVAPR", 5):
            case Map_Str2Int("HVR", 3):
            case Map_Str2Int("EV", 2):
            case Map_Str2Int("AVSIC", 5):
            case Map_Str2Int("BVSIC", 5):
            case Map_Str2Int("AVG", 3):
            case Map_Str2Int("BVG", 3):
                paramRs.InputCOMPONENTS(ifs, keyword);
                break;

            case Map_Str2Int("PRSR", 4):
            case Map_Str2Int("TEMR", 4):
                paramRs.InputRefPR(ifs, keyword);
                break;

            case Map_Str2Int("LBCCOEF", 7):
                paramRs.InputLBCCOEF(ifs);
                break;

            case Map_Str2Int("BIC", 3):
                paramRs.InputBIC(ifs);
                break;

            case Map_Str2Int("SSMSTA", 6):
                paramRs.InputSSMSTA(ifs);
                break;

            case Map_Str2Int("SSMSP", 5):
                paramRs.InputSSMSP(ifs);
                break;

            case Map_Str2Int("NRSTA", 5):
                paramRs.InputNRSTA(ifs);
                break;

            case Map_Str2Int("NRSP", 4):
                paramRs.InputNRSP(ifs);
                break;

            case Map_Str2Int("RR", 2):
                paramRs.InputRR(ifs);
                break;

            case Map_Str2Int("BOUNDARY", 8):
                paramRs.InputBoundary(ifs);
                break;

            case Map_Str2Int("CURTIME", 7):
                paramControl.InpuCurTime(ifs);
                break;

            case Map_Str2Int("MAXSTIME", 8):
                paramControl.InpuMaxSimTime(ifs);
                break;

            default: // skip non-keywords
                break;
        }
    }

    ifs.close();
}

void ParamRead::ReadFileHiSim(const string& filename)
{
    ifstream input(filename, ios::in);
    if (!input) {
        OCP_MESSAGE("Trying to open file: " << (filename));
        OCP_ABORT("Failed to open the input file!");
    }

    /// Read input data
    std::vector<std::string> words;
    std::string buff;
    std::vector<std::vector<std::string>> main_data;
    vector<int> include_file_pos;
    while (GetLineSkipComments(input, buff))
    {
        words = strip_split(buff);
        if (words[0] == "INCLUDE")
        {
            include_file_pos.push_back(main_data.size());
        }
        else
        {
            DealDefault(words);
            if (!words.empty())
                main_data.push_back(words);
        }
    }
    input.close();
    /// Read include files
    vector<vector<vector<string>>> include_data;
    for (int i=0; i<include_file_pos.size(); ++i)
    {
        ifstream include(workDir + main_data[include_file_pos[i]][0] + ".dat", ios::in);
        if (!include)
        {
            OCP_MESSAGE("Trying to open file: " << (workDir + main_data[include_file_pos[i]][0] + ".dat"));
            OCP_ABORT("Failed to open the input file!");
        }

        vector<vector<string>> each_include;
        while (GetLineSkipComments(include, buff))
        {
            words = strip_split(buff);
            DealDefault(words); // assume no include in include files
            if (!words.empty())
                each_include.push_back(words);
        }
        include.close();

        include_data.push_back(each_include);
    }
    /// Insert include data into main input data
    vector<vector<string>> input_data;
    int pos_idx = 0;
    for (int i=0; i<main_data.size(); ++i)
    {
        if (include_file_pos.size() > pos_idx && i == include_file_pos[pos_idx])
        {
            for (auto itm: include_data[pos_idx])
                input_data.push_back(itm);
            pos_idx++;
        }
        else
        {
            input_data.push_back(main_data[i]);
        }
    }

    /// Compute section starts
    int GRID_start = -1;
    int WELL_start = -1;
    int PROPS_start = -1;
    int SOLUTION_start = -1;
    int SCHEDULE_start = -1;
    int TUNE_start = -1;
    int MODEL_start = 0;
    for (int i=1; i<input_data.size(); ++i)
    {
        if (input_data[i][0] == "GRID")
            GRID_start = i;
        else if (input_data[i][0] == "WELL" && input_data[i].size() == 1)
            WELL_start = i;
        else if (input_data[i][0] == "PROPS")
            PROPS_start = i;
        else if (input_data[i][0] == "SCHEDULE")
            SCHEDULE_start = i;
        else if (input_data[i][0] == "SOLUTION")
            SOLUTION_start = i;
        else if (input_data[i][0] == "TUNE")
            TUNE_start = i;
    }

    /// Compute section ends
    std::vector<int> all_starts;
    all_starts.push_back(MODEL_start);
    if (GRID_start > -1) all_starts.push_back(GRID_start);
    if (WELL_start > -1) all_starts.push_back(WELL_start);
    if (PROPS_start > -1) all_starts.push_back(PROPS_start);
    if (SOLUTION_start > -1) all_starts.push_back(SOLUTION_start);
    if (SCHEDULE_start > -1) all_starts.push_back(SCHEDULE_start);
    if (TUNE_start > -1) all_starts.push_back(TUNE_start);
    //
    std::sort(all_starts.begin(), all_starts.end());
    //
    int MODEL_end = all_starts[1];
    //
    int GRID_end = input_data.size();
    for (auto itm: all_starts)
        if (itm > GRID_start)
        {
            GRID_end = itm;
            break;
        }
    //
    int WELL_end = input_data.size();
    for (auto itm: all_starts)
        if (itm > WELL_start)
        {
            WELL_end = itm;
            break;
        }
    //
    int PROPS_end = input_data.size();
    for (auto itm: all_starts)
        if (itm > PROPS_start)
        {
            PROPS_end = itm;
            break;
        }
    //
    int SOLUTION_end = input_data.size();
    for (auto itm: all_starts)
        if (itm > SOLUTION_start)
        {
            SOLUTION_end = itm;
            break;
        }
    //
    int SCHEDULE_end = input_data.size();
    for (auto itm: all_starts)
        if (itm > SCHEDULE_start)
        {
            SCHEDULE_end = itm;
            break;
        }
    //
    int TUNE_end = input_data.size();
    for (auto itm: all_starts)
        if (itm > TUNE_start)
        {
            TUNE_end = itm;
            break;
        }


    /***************************************************************************
     *
     *                    Read and Set parameters.
     *
     ***************************************************************************/

    /// MODEL section
    for (int i=MODEL_start; i<MODEL_end; ++i)
    {
        words = input_data[i];
        if (words[0] == "MODELTYPE")
        {
            if (words[1] == "BlackOil" || words[1] == "BLACKOIL")
                paramRs.blackOil = OCP_TRUE;
        }
        else if (words[0] == "SOLNT")
        {
        }
        else if (words[0] == "FIELD")
        {
            paramRs.unitType = words[0];
        }
        else
        {
            OCP_MESSAGE("Found not supported keywords in MODEL section: " << buff);
            OCP_ABORT("Wrong keywords!");
        }
    }


    /// WELL section
    std::vector<std::string> markers;
    std::vector<std::vector<std::string>> wells_data;
    for (int i=WELL_start+1; i<WELL_end; ++i)
    {
        words = input_data[i];
        if (words[0] == "TEMPLATE")
        {
            markers.push_back("NAME");
            do {
                i++;
                words = input_data[i];
                for (int j=0; j<words.size()-1; ++j)
                    markers.push_back(words[j]);

                if (words[words.size() - 1] == "/")
                    break;

                markers.push_back(words[words.size() - 1]);
            } while (input_data[i+1][0].find('\'') == 0);
        }
        else if (words[0] == "WELSPECS")
        {
            do {
                i++;
                string well_name = input_data[i][1]; // well name
                // set perforations of the well above
                do {
                    vector<string> each_perf;
                    each_perf.push_back(well_name);
                    i++;
                    each_perf.insert(each_perf.end(), input_data[i].begin(),
                                     input_data[i].end()); // parameters
                    wells_data.push_back(each_perf);
                } while (input_data[i+1][0] != "NAME" && input_data[i+1].size() == markers.size()-1);
            } while (input_data[i+1][0] == "NAME");
        }
        else if (words[0] == "MSWOPT")
        {
            i++;
            cout << input_data[i][0] << endl;
        }
        else
        {
            OCP_MESSAGE("Found not supported keywords in WELL section: " << buff);
            OCP_ABORT("Wrong keywords!");
        }
    }
    // Construct wells using markers and wells_data
    std::map<std::string, std::vector<std::string>> well_data_map;
    for (int i=0; i<markers.size(); ++i)
    {
        for (int j=0; j<wells_data.size(); ++j)
            well_data_map[markers[i]].push_back(wells_data[j][i]);
    }
//    if (gridType == GridType::structured || gridType == GridType::orthogonal) /// fff
    {
        for (int i=0; i<wells_data.size(); ++i) // loop over all wells
        {
            // set well header
            string name = well_data_map["NAME"][i];
            int I = stoi(well_data_map["'I'"][i]);
            int J = stoi(well_data_map["'J'"][i]);
            // find if the well exists
            bool found_well = false;
            int idx = paramWell.well.size();
            for (int j=0; j<paramWell.well.size(); ++j)
            {
                auto& wl = paramWell.well[j];
                if (wl.name == name)
                {
                    found_well = true;
                    idx = j;
                    break;
                }
            }
            if (!found_well)
                paramWell.well.push_back(WellParam(name, I, J));

            // set perforations of well
            double diam = 1.0; // default
            if (well_data_map.find("'DIAM'") != well_data_map.end())
                diam = stod(well_data_map["'DIAM'"][i]);
            else if (well_data_map.find("'RW'") != well_data_map.end())
                diam = stod(well_data_map["'RW'"][i]) * 2.0;
            int k1 = stoi(well_data_map["'K1'"][i]);
            int k2 = stoi(well_data_map["'K2'"][i]);
            for (int k=k1; k<=k2; ++k)
                paramWell.well[idx].SetWellParams(I, J, k, diam);
        }
    }


    /// PROPS section
    for (int i=PROPS_start+1; i<PROPS_end; ++i)
    {
        words = input_data[i];
        if (words[0] == "STCOND")
        {
            i++;
            words = input_data[i];

            paramRs.Tsurf = stod(words[0]);
            paramRs.Tsurf = paramWell.Tsurf;

            paramWell.Psurf = stod(words[1]);
            paramRs.Psurf = paramWell.Psurf;
        }
        else if (words[0] == "NCOMPS")
        {
            paramRs.SetCOMPS(stoi(words[1]));
        }
        else if (words[0] == "EOS")
        {

        }
        else if (words[0] == "PRCORR")
        {

        }
        else if (words[0] == "CNAMES")
        {
            i++;
            words = input_data[i];
            paramRs.SetCNAMES(words);
        }
        else if (words[0] == "EOSCOEF")
        {
            continue;
        }
        else if (words[0] == "RTEMP")
        {
            i++;
            words = input_data[i];
            paramRs.rsTemp = stod(words[0]);
            cout << "rtemp: " << words[0] << endl;
        }
        else if (words[0] == "TCRIT")
        {
            i++;
            words = input_data[i];
            paramRs.InputCOMPONENTS("TCRIT", words);
        }
        else if (words[0] == "PCRIT")
        {
            i++;
            words = input_data[i];
            paramRs.InputCOMPONENTS("PCRIT", words);
        }
        else if (words[0] == "ZCRIT")
        {
            i++;
            words = input_data[i];
            paramRs.InputCOMPONENTS("ZCRIT", words);
        }
        else if (words[0] == "MW")
        {
            i++;
            words = input_data[i];
            paramRs.InputCOMPONENTS("MW", words);
        }
        else if (words[0] == "ACF")
        {
            i++;
            words = input_data[i];
            paramRs.InputCOMPONENTS("ACF", words);
        }
        else if (words[0] == "BIC")
        {
            vector<double> bic;
            for (int j=0; j<paramRs.numCom-1; ++j)
            {
                i++;
                words = input_data[i];
                for (auto& itm: words)
                    bic.push_back(stod(itm));
            }
            paramRs.SetBIC(bic);
        }
        else if (words[0] == "STONE2")
        {
            continue;
        }
        else if (words[0] == "SWOF")
        {
            TableSet* obj = paramRs.FindPtrTable("SWOF");
            const int num_cols = obj->colNum;
            vector<vector<OCP_DBL>> swof(num_cols);
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                    for (int l=0; l<num_cols; ++l)
                        swof[l].push_back(stod(words[l]));
                else
                {
                    i--;
                    break;
                }
            }
            obj->data.push_back(swof);
        }
        else if (words[0] == "SGOF")
        {
            TableSet* obj = paramRs.FindPtrTable("SGOF");
            const int num_cols = obj->colNum;
            vector<vector<OCP_DBL>> sgof(num_cols);
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                    for (int l=0; l<num_cols; ++l)
                        sgof[l].push_back(stod(words[l]));
                else
                {
                    i--;
                    break;
                }
            }
            obj->data.push_back(sgof);
        }
        else if (words[0] == "PVDG")
        {
            TableSet* obj = paramRs.FindPtrTable("PVDG");
            const int num_cols = obj->colNum;
            vector<vector<OCP_DBL>> pvdg(num_cols);
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    for (int l=0; l<num_cols; ++l)
                        pvdg[l].push_back(stod(words[l]));
                }
                else
                {
                    i--;
                    break;
                }
            }
            obj->data.push_back(pvdg);
        }
        else if (words[0] == "PVCO")
        {
            TableSet* obj = paramRs.FindPtrTable("PVCO");
            const int num_cols = obj->colNum;
            vector<vector<OCP_DBL>> pvco(num_cols);
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    for (int l=0; l<num_cols; ++l)
                        pvco[l].push_back(stod(words[l]));
                }
                else
                {
                    if (words[0] != "/")
                        i--;
                    break;
                }
            }
            obj->data.push_back(pvco);
        }
        else if (words[0] == "PVTW")
        {
            TableSet* obj = paramRs.FindPtrTable("PVTW");
            const int num_cols = obj->colNum;
            vector<vector<OCP_DBL>> pvtw(num_cols); /// TODO not string, double
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    for (int l=0; l<num_cols; ++l)
                        pvtw[l].push_back(stod(words[l]));
                }
                else
                {
                    i--;
                    break;
                }
            }
            obj->data.push_back(pvtw);
        }
        else if (words[0] == "ROCK")
        {
            i++;
            words = input_data[i];
            double p_ref = stod(words[0]);
            double cp1 = stod(words[1]);
            paramRs.SetROCK(p_ref, cp1);
        }
        else if (words[0] == "ROCKTAB")
        {
            int num_cols = 5; // 五列数据: 流体压力或有效应力, 孔隙度乘数, X,Y,Z方向渗透率乘数
            vector<vector<double>> rocktab(num_cols);
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    for (int j=0; j<num_cols; ++j)
                        rocktab[j].push_back(stod(words[j]));
                }
                else
                {
                    i--;
                    break;
                }
            }

//            for (int j=0; j<rocktab[0].size()-1; ++j)
            for (int j=0; j<1; ++j) /// fff
            {
                double pref = rocktab[0][j];
                double cp1 = (rocktab[1][j+1] - rocktab[1][j]) / (rocktab[0][j+1] - rocktab[0][j]);
                paramRs.SetROCK(pref, cp1);
            }
        }
        else if (words[0] == "WATERTAB")
        {
            vector<vector<double>> params(3); // 三列数据: 压力 p, 水的体积系数 Bw, 水的粘度
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    for (int j=0; j<3; ++j)
                        params[j].push_back(stod(words[j]));
                }
                else
                {
                    i--;
                    break;
                }
            }

            vector<vector<double>> pvtw(5); // fff
            pvtw[0].push_back(14.7);
            pvtw[1].push_back(1.0);
            pvtw[2].push_back(-1.0 * (params[1][1] - params[1][0]) / (params[0][1] - params[0][0]));
            pvtw[3].push_back(params[2][0]);
            pvtw[4].push_back(0.0);

            paramRs.SetPVTW(pvtw);
        }
        else if (words[0] == "DENSITY")
        {
            i++;
            words = input_data[i];
            paramRs.SetDENSITY(words);
        }
        else if (words[0] == "PVTO")
        {
            OCP_WARNING("Not support keyword: PVTO!");
//            TableSet* obj = paramRs.FindPtrTable("PVDG");
//            const int num_cols = obj->colNum;
            const int num_cols = 4;
            vector<vector<OCP_DBL>> pvto(num_cols);
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    for (int l=0; l<num_cols; ++l)
                    {
                        double val;
                        try {
                            val = stod(words[l]);
                        }
                        catch (std::exception)
                        {
                            val = 0.0;
                        }
                        pvto[l].push_back(val);
                    }
                }
                else
                {
                    i--;
                    break;
                }
            }
//            obj->data.push_back(pvdg);
        }
        else if (words[0] == "PVDO")
        {
            TableSet* obj = paramRs.FindPtrTable("PVDO");
            const int num_cols = obj->colNum;
            vector<vector<OCP_DBL>> pvdo(num_cols);
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    for (int l=0; l<num_cols; ++l)
                        pvdo[l].push_back(stod(words[l]));
                }
                else
                {
                    i--;
                    break;
                }
            }
            obj->data.push_back(pvdo);
        }
        else if (words[0] == "STONEII")
            OCP_WARNING("Not support keyword: STONEII!");
        else if (words[0] == "/")
            continue;
        else
        {
            OCP_ABORT("Not support keywords!");
        }
    }
    // tabdims fff
    paramRs.NTSFUN = 1;
    paramRs.NTPVT = 1;
    paramRs.NTROOC = 1;
    paramRs.comsParam.NTPVT = 1;

    /// SOLUTION section
    for (int i=SOLUTION_start+1; i<SOLUTION_end; ++i)
    {
        words = input_data[i];

        if (words[0] == "EQUILPAR") /// fff EQUIL
        {
            vector<string> equilpar; /// TODO not string, double
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    for (int j=0; j<words.size(); ++j)
                    {
                        if (j == 4) // 跳过: 初始化时的深度步长 dh
                            continue;
                        equilpar.push_back(words[j]);
                    }
                }
                else
                {
                    i--;
                    break;
                }
            }

            vector<string> tmp_equilpar;
            if (equilpar.size() >= 6)
                tmp_equilpar.assign(equilpar.begin(), equilpar.begin()+6); /// 暂时支持6个参数
            else
            {
                int rest = 6 - equilpar.size();
                tmp_equilpar.assign(equilpar.begin(), equilpar.end());
                for (int j=0; j<rest; ++j)
                    tmp_equilpar.push_back("NA");
            }
            paramRs.SetEQUIL(tmp_equilpar);
            paramRs.initType = "EQUIL";
        }
        else if (words[0] == "PBVD")
        {
            TableSet* obj = paramRs.FindPtrTable("PBVD");
            const int num_cols = obj->colNum;

            vector<vector<OCP_DBL>> pbvd(num_cols);
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    for (int j=0; j<num_cols; ++j)
                    {
                        pbvd[j].push_back(std::stod(words[j]));
                    }
                    cout << endl;
                }
                else
                {
                    if (words[0] != "/")
                        i--;
                    break;
                }
            }

            obj->data.push_back(pbvd);
        }
        else if (words[0] == "ZMFVD")
        {
            TableSet* obj = paramRs.FindPtrTable("ZMFVD");
            const int num_cols = obj->colNum;

            vector<vector<OCP_DBL>> zmfvd(num_cols);
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    for (int j=0; j<num_cols; ++j)
                    {
                        zmfvd[j].push_back(std::stod(words[j]));
                    }
                    cout << endl;
                }
                else
                {
                    i--;
                    break;
                }
            }

            obj->data.push_back(zmfvd);
        }
        else if (words[0] == "METHOD")
        {
            i++;
            words = input_data[i];
            paramControl.method.clear();
            paramControl.lsFile.clear();
            paramControl.method.push_back(words[0]);
            paramControl.lsFile.push_back(words[1]);
        }
        else if (words[0] == "SUMMARY") // read in ReadRestFile()
        {
            do {
                i++;
                words = input_data[i];
            } while (input_data[i+1][0] != "/");
        }
        else if (words[0] == "/")
            break;
        else
        {
            OCP_ABORT("Not support keywords!");
        }
    }

    /// TUNE section fff
    for (int i=TUNE_start+1; i<TUNE_end; ++i)
    {
        vector<vector<OCP_DBL>> tmp(paramControl.tuning);
        do {
            words = input_data[i];
            i++;
        } while (i < TUNE_end);
//        paramControl.tuning_T.push_back(TuningPair(0, tmp)); // fff
    }

    /// SCHEDULE section
    double now = 0.0;
    double is_first_date = true;
    string first_date;
    vector<double> all_tsteps;
    for (int i=SCHEDULE_start+1; i<SCHEDULE_end; ++i)
    {
        words = input_data[i];

        if (words[0] == "USEENDTIME" || words[0] == "USESTARTTIME")
        {
            if (words[0] == "USEENDTIME")
                recurrent_type = 1;
            else
                recurrent_type = 0; // same with OCP
        }
        else if (words[0] == "RPTSCHED")
            continue;
        else if (words[0] == "RECURRENT")
            continue;
        else if (words[0] == "BASIC")
            continue;
        else if (words[0] == "TIME")
        {
            vector<double> tsteps;
            if (is_first_date)
            {
                first_date = words[1];
                is_first_date = false;
                tsteps.push_back(0);
                all_tsteps.push_back(0.0);

                if (WhichDateFormat(words[1]) == 1)
                {
                    for (int j=1; j<words.size(); ++j)
                    {
                        now += stod(words[j]);
                        tsteps.push_back(now);
                        all_tsteps.push_back(now);
                    }
                }
                else
                {
                    OCP_ABORT("Not support this time format!");
                }
            }
            else
            {
                if (words.size() == 1)
                    OCP_WARNING("Not support this time format!");
                else
                {
                    for (int j=1; j<words.size(); ++j)
                    {
                        if (WhichDateFormat(words[j]) == 1)
                        {
                            now += stod(words[j]);
                            tsteps.push_back(now);
                            all_tsteps.push_back(now);
                        }
                        else
                        {
                            int days = NumDaysBetweenDates(first_date, words[j]);
                            now = days;
                            tsteps.push_back(days);
                            all_tsteps.push_back(now);
                        }
                    }
                }

            }

            while (input_data[i+1][0] == "WELL")
            {
                words = input_data[i+1];
                string well_name = words[1];
                string state = "OPEN"; /// fff
                string mode = words[2];
                string max_rate = words[3];

                if (TypeMap[mode] == WellOptParam::PROD)
                {
                    double min_bhp = stod(words[4]);
                    WellOptParam opt(state, mode, max_rate, min_bhp);
                    paramWell.well_oper_list.push_back(ParamWell::WellOperation(now,
                                                                                well_name,
                                                                                "PROD",
                                                                                opt));
                }
                else if (TypeMap[mode] == WellOptParam::INJ)
                {
                    string fluidType;
                    if (words[2] == "WIR")
                        fluidType = "WATER";
                    else if (words[2] == "GIR")
                    {
                        // fff
                        fluidType = "solvent" + std::to_string(paramWell.solSet.size());
                    }
                    else
                        OCP_ABORT("Not support fluidType!");

                    if (paramRs.blackOil)
                        fluidType = "GAS";

                    double max_bhp = stod(words[4]);

                    if (mode == "WIR" || mode == "GIR")
                        mode = "RATE";

                    WellOptParam opt(fluidType, state, mode, max_rate, max_bhp); /// fff stream
                    paramWell.well_oper_list.push_back(ParamWell::WellOperation(now,
                                                                                well_name,
                                                                                "INJ",
                                                                                opt));
                }
                else
                    OCP_ABORT("Wrong well operation type!");
                i++;
            }
        }
        else if (words[0] == "WELLSCHED")
        {
            double now_tmp = 0.0;
            double start_tmp=0.0;

            string well_name = words[1];
            while (1)
            {
                i++;
                words = input_data[i];
                if (words[0] == "TIME")
                {
                    int num_times=1;
                    start_tmp += stod(words[1]);
                    all_tsteps.push_back(start_tmp);
                    if (words.size() > 2)
                    {
                        for (int l=2; l<words.size(); ++l)
                        {
                            if (words[l] == words[1])
                            {
                                num_times++;
                                start_tmp += stod(words[l]);
                                all_tsteps.push_back(start_tmp);
                            }
                            else
                                break;
                        }
                    }

                    if (num_times == 1)
                    {
                        /// 生成一个WellOptParam对象
                        now_tmp += stod(words[1]); /// time

                        string state = "OPEN";
                        string mode = words[2];
                        string max_rate = words[3];

                        if (TypeMap[mode] == WellOptParam::PROD)
                        {
                            double min_bhp;
                            if (words[4] == "BHP")
                                min_bhp = stod(words[5]);
                            else
                                min_bhp = stod(words[4]);
                            WellOptParam opt(state, mode, max_rate, min_bhp);
                            paramWell.well_oper_list.push_back(ParamWell::WellOperation(now_tmp,
                                                                                        well_name,
                                                                                        "PROD",
                                                                                        opt));
                        }
                        else if (TypeMap[mode] == WellOptParam::INJ)
                        {
                            string fluidType;
                            if (words[2] == "WIR")
                                fluidType = "WATER";
                            else if (words[2] == "GIR")
                            {
                                // fff
                                if (paramWell.solSet.size() <= 1)
                                    fluidType = "solvent" + std::to_string(0);
                                else
                                    fluidType = "solvent" + std::to_string(paramWell.solSet.size());
                            }
                            else
                                OCP_ABORT("Not support fluidType!");

                            if (paramRs.blackOil)
                                fluidType = "GAS";

                            if (mode == "GIR" || mode == "WIR")
                                mode = "RATE";

                            double max_bhp;
                            if (words[4] == "BHP")
                                max_bhp = stod(words[5]);
                            else
                                max_bhp = stod(words[4]);
                            WellOptParam opt(fluidType, state, mode, max_rate, max_bhp);
                            paramWell.well_oper_list.push_back(ParamWell::WellOperation(now_tmp,
                                                                                        well_name,
                                                                                        "INJ",
                                                                                        opt));
                        }
                        else
                            OCP_ABORT("Wrong well operation type!");
                    }
                    else
                    {
                        for (int k=1; k<num_times; ++k)
                            now_tmp += stod(words[1]);

                        if (std::find(words.begin(), words.end(), "SHUT") != words.end())
                            continue;

                        /// 生成一个WellOptParam对象
                        string state = "OPEN";
                        string mode = words[num_times+1];
                        string max_rate = words[num_times+2];
                        double min_bhp;
                        if (words[num_times+3] == "BHP")
                            min_bhp = stoi(words[num_times+4]);
                        else
                            cout << "### ERROR: Prod Rate is missing in WCONINJE!" << endl;

                        if (TypeMap[mode] == WellOptParam::PROD)
                        {
                            WellOptParam opt(state, mode, max_rate, min_bhp);
                            paramWell.well_oper_list.push_back(ParamWell::WellOperation(now_tmp,
                                                                                        well_name,
                                                                                        "PROD",
                                                                                        opt));
                        }
                        else if (TypeMap[mode] == WellOptParam::INJ)
                        {
                            string fluidType;
                            if (words[2] == "WIR")
                                fluidType = "WATER";
                            else if (words[2] == "GIR")
                            {
                                // fff
                                fluidType = "solvent" + std::to_string(paramWell.solSet.size());
                            }
                            else
                                OCP_ABORT("Not support fluidType!");

                            if (paramRs.blackOil)
                                fluidType = "GAS";

                            if (mode == "GIR" || mode == "WIR")
                                mode = "RATE";

                            WellOptParam opt(fluidType, state, mode, max_rate, min_bhp);
                            paramWell.well_oper_list.push_back(ParamWell::WellOperation(now_tmp,
                                                                                        well_name,
                                                                                        "INJ",
                                                                                        opt));
                        }
                        else
                            OCP_ABORT("Wrong well operation type!");
                    }
                }
                else if (words[0] == "LIMIT") /// fff
                {
                    string opt = words[1];
                }
                else if (words[0] == "PERF") /// fff
                {
                    string opt = words[1];
                }
                else if (words[0] == "STREAM")
                {
                    string opt = words[1];
                    vector<double> ratios;
                    for (int l=0; l<paramRs.numCom; ++l)
                    {
                        if (l < words.size()-1)
                            ratios.push_back(stod(words[l+1]));
                        else
                            ratios.push_back(0.0);
                    }
                    paramWell.SetSTREAM(ratios);
                }
                else
                {
                    i--;
                    break;
                }
            }
        }
        else if (words[0] == "RESTART")
        {

        }
        else if (words[0] == "RPTSUM")
        {
            do {
                i++;
                string var = input_data[i][0];
            } while (input_data[i+1][0] != "/");
        }
        else if (words[0] == "/")
            break;
        else if (words[0] == "BINOUT" || words[0] == "POIL") // fff
            continue;
        else
        {
            OCP_MESSAGE("Found not supported keywords in SCHEDULE section: " << words[0]);
            OCP_ABORT("Wrong keywords!");
        }
    }

    /// Establish Time and Well Operations
    /// Sort and unique
    std::sort(all_tsteps.begin(), all_tsteps.end());
    all_tsteps.erase(std::unique(all_tsteps.begin(), all_tsteps.end()), all_tsteps.end());
    /// construct criticalTime
    if (std::abs(all_tsteps[0] - 0.0) > 1.0E-10)
        paramWell.criticalTime.push_back(all_tsteps[0]);
    //
    for (int i=1; i<all_tsteps.size(); ++i)
    {
        paramWell.criticalTime.push_back(all_tsteps[i]);
    }
    //
    paramControl.criticalTime = paramWell.criticalTime;
    /// construct
    for (int i=0; i<paramWell.well_oper_list.size(); ++i)
    {
        ParamWell::WellOperation opt = paramWell.well_oper_list[i];

        double t = opt.tstep;
        vector<double>::iterator iter = std::find(paramWell.criticalTime.begin(),
                                                  paramWell.criticalTime.end(), t);
        if (iter == paramWell.criticalTime.end())
            OCP_ABORT("Wrong time step!");
        int idx_criticalTime = std::distance(paramWell.criticalTime.begin(), iter);

        int idx_well = -1;
        string name = opt.name;
        for (int j=0; j<paramWell.well.size(); ++j)
        {
            if (paramWell.well[j].name == name)
                idx_well = j;
        }
        if (idx_well == -1)
            OCP_ABORT("Wrong well name!");

        paramWell.well[idx_well].optParam.push_back(WellOptPair(idx_criticalTime, opt.opt));
    }
}

void ParamRead::ReadRestFile(const string& filename)
{
    ifstream ifs(filename, ios::in);
    if (!ifs) {
        OCP_MESSAGE("Trying to open file: " << (filename));
        OCP_ABORT("Failed to open the input file!");
    }

    while (!ifs.eof()) {
        vector<string> vbuf;
        if (!ReadLine(ifs, vbuf)) break;
        string keyword = vbuf[0];

        switch (Map_Str2Int(&keyword[0], keyword.size()))
        {
            case Map_Str2Int("METHOD", 6):
                paramControl.InputMETHOD(ifs);
                break;

            case Map_Str2Int("TUNING", 6):
                paramControl.InputTUNING(ifs);
                break;

            case Map_Str2Int("SUMMARY", 7):
                paramOutput.InputSUMMARY(ifs);
                break;

            default: // skip non-keywords
                break;
        }
    }

    ifs.close();
}


/// Read INCLUDE files; these files should have identical format.
void ParamRead::ReadINCLUDE(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    DealDefault(vbuf);

    if (CURRENT_RANK == MASTER_PROCESS)
        cout << "begin to read " + workDir + vbuf[0] << endl;
    ReadFile(workDir + vbuf[0]);
    if (CURRENT_RANK == MASTER_PROCESS)
        cout << "finish reading " + workDir + vbuf[0] << endl;
}



/// Check parameters in paramRs and paramWell.
void ParamRead::CheckParam()
{
    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        OCP_INFO("Check reading parameters from input data -- begin");
    }

    paramRs.CheckParam();
    paramWell.CheckParam();

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        OCP_INFO("Check reading parameters from input data -- end");
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/08/2022      Test robustness for wrong keywords   */
/*----------------------------------------------------------------------------*/