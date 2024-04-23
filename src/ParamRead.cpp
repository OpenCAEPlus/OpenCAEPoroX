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
        ReadFileHiSim(inputFile);
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

    /// TODO 第二处关键字读取
    while (!ifs.eof()) {
        vector<string> vbuf;
        if (!ReadLine(ifs, vbuf)) break;
        string keyword = vbuf[0];

        std::cout << "第二阶段 关键字: " << keyword << std::endl;
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
    string buff;
    std::vector<std::vector<std::string>> input_data;
    while (GetLineSkipComments(input, buff))
    {
        words = strip_split(buff);
        if (!words.empty())
            input_data.push_back(words);
    }
    input.close();

    /// Compute section starts
    int GRID_start = -1;
    int WELL_start = -1;
    int PROPS_start = -1;
    int SOLUTION_start = -1;
    int SCHEDULE_start = -1;
    int TUNE_start = -1;
    assert(input_data[0][0] == "MODELTYPE");
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

    /// MODEL section
    for (int i=MODEL_start; i<MODEL_end; ++i)
    {
        words = input_data[i];
        if (words[0] == "MODELTYPE")
        {
            string model_type = words[1]; // GASWATER, OILWATER, BLACKOIL, COMP
//            model = OCPModel::isothermal; /// fff HiSim only have isothermal?
            continue;
        }
        else if (words[0] == "SOLNT")
        {
            int omp = std::stoi(words[1]);
            cout << "并行线程数: " << omp << endl;
        }
        else if (words[0] == "FIELD")
        {
            cout << "单位制: " << words[0] << endl;
        }
        else
        {
            OCP_MESSAGE("Found not supported keywords: " << buff);
            OCP_ABORT("Wrong keywords!");
        }
    }

    /// GRID section
    for (int i=GRID_start+1; i<GRID_end; ++i)
    {
        words = input_data[i];
        if (words[0] == "DIMENS")
        {
//            if (words.size() > 1)
//            {
//                nx = stoi(words[1]);
//                ny = stoi(words[2]);
//                nz = stoi(words[3]);
//            }
//            else
//            {
//                i += 1;
//                words = input_data[i];
//                nx = stoi(words[0]);
//                ny = stoi(words[1]);
//                nz = stoi(words[2]);
//            }
//
//            numGridM = nx * ny * nz;
//            if (DUALPORO) numGridF = numGridM; /// fff how to set DUALPORO
//
//            numGrid = numGridM + numGridF;
            continue;
        }
        else if (words[0] == "DXV")
        {
            if (words.size() > 1)
            {
                cout << "DXV: " << words[1] << endl;
            }
            else
            {
                i++;
                words = input_data[i];
                cout << "DXV: " << words[0] << endl;
            }
        }
        else if (words[0] == "DYV")
        {
            if (words.size() > 1)
            {
                cout << "DYV: " << words[1] << endl;
            }
            else
            {
                i++;
                words = input_data[i];
                cout << "DYV: " << words[0] << endl;
            }
        }
        else if (words[0] == "DZV")
        {
            if (words.size() > 1)
                cout << "DZV: " << words[1] << ", " << words[2] << ", " << words[3] << endl;
            else
            {
                i += 1;
                words = input_data[i];
                cout << "DZV, " << words[0] << ", " << words[1] << ", " << words[2] << endl;
            }
        }
        else if (words[0] == "TOPS")
        {
            if (words.size() > 1)
            {
                cout << "TOPS: " << words[1] << endl;
            }
            else
            {
                i++;
                words = input_data[i];
                cout << "TOPS: " << words[0] << endl;
            }
        }
        else if (words[0] == "PORO")
        {
            if (words.size() > 1)
            {
                cout << "PORO: " << words[1] << endl;
            }
            else
            {
                i++;
                words = input_data[i];
                cout << "PORO: " << words[0] << endl;
            }
        }
        else if (words[0] == "PERMX")
        {
            if (words.size() > 1)
            {
                cout << "PERMX: " << words[1] << " " << words[2] << " " << words[3] << endl;
            }
            else
            {
                i++;
                words = input_data[i];
                cout << "PERMX: " << words[0] << " " << words[1] << " " << words[2] << endl;
            }
        }
        else if (words[0] == "PERMY")
        {
            if (words.size() > 1)
            {
                cout << "PERMY: " << words[1] << " " << words[2] << " " << words[3] << endl;
            }
            else
            {
                i++;
                words = input_data[i];
                cout << "PERMY: " << words[0] << " " << words[1] << " " << words[2] << endl;
            }
        }
        else if (words[0] == "PERMZ")
        {
            if (words.size() > 1)
            {
                cout << "PERMZ: " << words[1] << " " << words[2] << " " << words[3] << endl;
            }
            else
            {
                i++;
                words = input_data[i];
                cout << "PERMZ: " << words[0] << " " << words[1] << " " << words[2] << endl;
            }
        }
        else
        {
            OCP_ABORT("Not support keywords!");
        }
    }

    /// WELL section
    for (int i=WELL_start+1; i<WELL_end; ++i)
    {
        words = input_data[i];
        if (words[0] == "TEMPLATE")
        {
            do {
                i++;
                words = input_data[i];
                if (words[words.size() - 1] == "/")
                    break;
            } while (input_data[i+1][0].find('\'') == 0);
        }
        else if (words[0] == "WELSPECS")
        {
            do {
                i++;
                words = input_data[i];
                i++;
                words = input_data[i];
            } while (input_data[i+1][0] == "NAME");
        }
        else
        {
            OCP_ABORT("Not support keywords!");
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
            cout << "STCOND: " << words[0] << ", " << words[1] << endl;
        }
        else if (words[0] == "NCOMPS")
        {
            cout << "NCOMPS: " << words[1] << endl;
        }
        else if (words[0] == "EOS")
        {
            cout << "EOS: " << words[1] << endl;
        }
        else if (words[0] == "PRCORR")
        {
            cout << "PRCORR" << endl;
        }
        else if (words[0] == "CNAMES")
        {
            i++;
            words = input_data[i];
            cout << "cnames: " << words[0] << ", " << words[5] << endl;
        }
        else if (words[0] == "EOSCOEF")
        {
            continue;
        }
        else if (words[0] == "RTEMP")
        {
            i++;
            words = input_data[i];
            cout << "rtemp: " << words[0] << endl;
        }
        else if (words[0] == "TCRIT")
        {
            i++;
            words = input_data[i];
            cout << "TCRIT: " << words[0] << endl;
        }
        else if (words[0] == "PCRIT")
        {
            i++;
            words = input_data[i];
            cout << "PCRIT: " << words[0] << endl;
        }
        else if (words[0] == "ZCRIT")
        {
            i++;
            words = input_data[i];
            cout << "ZCRIT: " << words[0] << endl;
        }
        else if (words[0] == "MW")
        {
            i++;
            words = input_data[i];
            cout << "MW: " << words[0] << endl;
        }
        else if (words[0] == "ACF")
        {
            i++;
            words = input_data[i];
            cout << "ACF: " << words[0] << endl;
        }
        else if (words[0] == "BIC")
        {
            vector<vector<string>> bic;  /// TODO not string, is double
            int NCOMPS = 6;
            for (int j=0; j<NCOMPS-1; ++j)
            {
                i++;
                words = input_data[i];
                cout << "bic: " << words[0] << endl;
                bic.push_back(words);
            }
        }
        else if (words[0] == "STONE2")
        {
            continue;
        }
        else if (words[0] == "SWOF")
        {
            vector<vector<string>> swof; /// TODO not string, double
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    cout << "swof: " << words[0] << endl;
                    swof.push_back(words);
                }
                else
                {
                    i--;
                    break;
                }
            }
        }
        else if (words[0] == "SGOF")
        {
            vector<vector<string>> sgof; /// TODO not string, double
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    cout << "sgof: " << words[0] << endl;
                    sgof.push_back(words);
                }
                else
                {
                    i--;
                    break;
                }
            }
        }
        else if (words[0] == "PVDG")
        {
            vector<vector<string>> pvdg; /// TODO not string, double
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    cout << "pvdg: " << words[0] << endl;
                    pvdg.push_back(words);
                }
                else
                {
                    i--;
                    break;
                }
            }
        }
        else if (words[0] == "PVCO")
        {
            vector<vector<string>> pvco; /// TODO not string, double
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    pvco.push_back(words);
                }
                else
                {
                    if (words[0] != "/")
                        i--;
                    break;
                }
            }
        }
        else if (words[0] == "PVTW")
        {
            vector<vector<string>> pvtw; /// TODO not string, double
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    pvtw.push_back(words);
                }
                else
                {
                    i--;
                    break;
                }
            }
        }
        else if (words[0] == "ROCK")
        {
            vector<vector<string>> rock; /// TODO not string, double
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    rock.push_back(words);
                }
                else
                {
                    i--;
                    break;
                }
            }
        }
        else if (words[0] == "ROCKTAB")
        {
            vector<vector<string>> rocktab; /// TODO not string, double
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    cout << "rocktab: " << words[0] << endl;
                    rocktab.push_back(words);
                }
                else
                {
                    i--;
                    break;
                }
            }
        }
        else if (words[0] == "WATERTAB")
        {
            vector<vector<string>> watertab; /// TODO not string, double
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    cout << "watertab: " << words[0] << endl;
                    watertab.push_back(words);
                }
                else
                {
                    i--;
                    break;
                }
            }
        }
        else if (words[0] == "DENSITY")
        {
            i++;
            words = input_data[i];
            cout << "density: " << words[0] << ", " << words[1] << ", " << words[2] << endl;
        }
        else
        {
            OCP_ABORT("Not support keywords!");
        }
    }

    /// SOLUTION section
    for (int i=SOLUTION_start+1; i<SOLUTION_end; ++i)
    {
        words = input_data[i];

        if (words[0] == "EQUILPAR")
        {
            vector<vector<string>> equilpar; /// TODO not string, double
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    equilpar.push_back(words);
                }
                else
                {
                    i--;
                    break;
                }
            }
            cout << "EQUILPAR: " << equilpar[0][0] << endl;
        }
        else if (words[0] == "PBVD")
        {
            vector<vector<string>> pbvd; /// TODO not string, double
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    pbvd.push_back(words);
                }
                else
                {
                    if (words[0] != "/")
                        i--;
                    break;
                }
            }
            cout << "EQUILPAR: " << pbvd[1][1] << endl;
        }
        else if (words[0] == "ZMFVD")
        {
            vector<vector<string>> zmfvd; /// TODO not string, double
            while (1)
            {
                i++;
                words = input_data[i];
                if (!isRegularString(words[0]))
                {
                    cout << "zmfvd: " << words[0] << endl;
                    zmfvd.push_back(words);
                }
                else
                {
                    i--;
                    break;
                }
            }
        }
        else
        {
            OCP_ABORT("Not support keywords!");
        }
    }

    /// TUNE section
    for (int i=TUNE_start+1; i<TUNE_end; ++i)
    {
        do {
            words = input_data[i];
            cout << words[0] << endl;
            i++;
        } while (i < TUNE_end);
    }

    /// SCHEDULE section
    for (int i=SCHEDULE_start+1; i<SCHEDULE_end; ++i)
    {
        words = input_data[i];

        if (words[0] == "USEENDTIME")
            continue;
        else if (words[0] == "USESTARTTIME")
            continue;
        else if (words[0] == "RPTSCHED")
            continue;
        else if (words[0] == "RECURRENT")
            continue;
        else if (words[0] == "BASIC")
            continue;
        else if (words[0] == "TIME")
        {
            cout << "time step: " << words[1] << endl;
            while (input_data[i+1][0] == "WELL")
            {
                words = input_data[i+1];
                cout << "well name: " << words[1] << endl;
                i++;
            }
        }
        else if (words[0] == "WELLSCHED")
        {
            cout << "wellsched: " << words[1] << endl;
            while (1)
            {
                i++;
                words = input_data[i];
                if (words[0] == "TIME")
                {
                    string opt = words[1];
                    cout << "TIME: " << opt << endl;
                }
                else if (words[0] == "LIMIT")
                {
                    string opt = words[1];
                    cout << "LIMIT: " << opt << endl;
                }
                else if (words[0] == "PERF")
                {
                    string opt = words[1];
                    cout << "PERF: " << opt << endl;
                }
                else if (words[0] == "STREAM")
                {
                    string opt = words[1];
                    cout << "STREAM: " << opt << endl;
                }
                else
                {
                    i--;
                    break;
                }
            }
        }
        else
        {
            OCP_MESSAGE("Found not supported keywords: " << words[0]);
            OCP_ABORT("Wrong keywords!");
        }
    }

    cout << "Good 2" << endl;
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