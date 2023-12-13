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

#include "ParamRead.hpp"

/// Initialize paramRs, paramWell, and paramControl.
void ParamRead::Init()
{
    paramRs.Init();
    paramWell.Init();
    paramControl.Init(workDir, fileName);
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
void ParamRead::ReadInputFile(const string& filename)
{
    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Input Reservoir Param -- begin");
    }

    inputFile = filename;
    GetDirAndName();
    Init();
    ReadFile(inputFile);
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

            case Map_Str2Int("MAXSTIME", 8):
                paramControl.InpuMaxSimTime(ifs);
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
        cout << endl
            << "=========================================" << endl
            << "Check reading parameters from input data!" << endl
            << "=========================================" << endl;
    }

    paramRs.CheckParam();
    paramWell.CheckParam();
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