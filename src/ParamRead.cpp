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
    out << paramRs.unitType << indent;

    //
    out << paramRs.blackOil << indent;

    // InputCOMPS
    out << paramRs.comps << indent
        << paramRs.numCom << indent
        << paramRs.comsParam.numCom << indent;
    for (auto val: paramRs.comsParam.LBCcoef)
        out << val << indent;
    for (auto val: paramRs.comsParam.SSMparamSTA)
        out << val << indent;
    for (auto val: paramRs.comsParam.NRparamSTA)
        out << val << indent;
    for (auto val: paramRs.comsParam.SSMparamSP)
        out << val << indent;
    for (auto val: paramRs.comsParam.NRparamSP)
        out << val << indent;
    for (auto val: paramRs.comsParam.RRparam)
        out << val << indent;

    // THERMAL
    out << paramRs.thermal << indent
        << paramWell.thermal << indent
        << indent;
    //
    out << paramRs.oil << indent
        << paramRs.gas <<indent
        << paramRs.water << indent
        << paramRs.disGas << indent
        << paramRs.GRAVDR << indent
        << paramRs.initType << indent
        << paramRs.initType << indent
        << '\n';

    // RTEMP
    out << paramRs.rsTemp << '\n';

    // InputTABLE
    


    for (auto val: dx)
        out << val << indent;

}


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


    /***************************************************************************
     *
     *                    Read and Set parameters.
     *
     ***************************************************************************/

    /// MODEL section
    for (int i=MODEL_start; i<MODEL_end; ++i)
    {
        words = input_data[i];
        if (words[0] == "MODELTYPE" || words[0] == "SOLNT")
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
                std::vector<std::string> each_well;
                i++;
                each_well.push_back(input_data[i][1]); // name
                i++;
                each_well.insert(each_well.end(), input_data[i].begin(), input_data[i].end()); // parameters
                wells_data.push_back(each_well);
            } while (input_data[i+1][0] == "NAME");
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
            string name = well_data_map["NAME"][i];
            int I = stoi(well_data_map["'I'"][i]);
            int J = stoi(well_data_map["'J'"][i]);
            paramWell.well.push_back(WellParam(name, I, J));

            double diam = stod(well_data_map["'DIAM'"][i]);
            int k1 = stoi(well_data_map["'K1'"][i]);
            int k2 = stoi(well_data_map["'K2'"][i]);
            for (int k=k1; k<=k2; ++k)
                paramWell.well[paramWell.well.size()-1].SetWellParams(I, J, k, diam);
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

            paramWell.Psurf = stod(words[0]);
            paramRs.Psurf = paramWell.Psurf;

            paramRs.Tsurf = stod(words[1]);
            paramRs.Tsurf = paramWell.Tsurf;
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
            TableSet* obj = paramRs.FindPtrTable("PVCO");
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
        else
        {
            OCP_ABORT("Not support keywords!");
        }
    }

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
                    for (int j=0; j<6; ++j) /// 暂时支持6个参数
                        equilpar.push_back(words[j]);
                }
                else
                {
                    i--;
                    break;
                }
            }

            paramRs.SetEQUIL(equilpar);
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
            double now = 0.0;

            string well_name = words[1];
            cout << "wellsched: " << well_name << endl;
            while (1)
            {
                i++;
                words = input_data[i];
                if (words[0] == "TIME")
                {
                    vector<vector<string>> tmp = ExpandWellOptions(words);

                    if (tmp.size() == 1)
                    {
                        /// 生成一个WellOptParam对象
                        now += stod(words[1]); /// time

                        string state = "OPEN";
                        string mode = words[2];
                        string max_rate = words[3];
                        double min_bhp;
                        if (words[4] == "BHP")
                            min_bhp = stoi(words[5]);
                        else
                            cout << "### ERROR: Prod Rate is missing in WCONINJE!" << endl;
                        WellOptParam opt(state, mode, max_rate, min_bhp);

                        paramWell.well_oper_list.push_back(ParamWell::WellOperation(now,
                                                                                    well_name,
                                                                                    "PROD", /// fff怎么知道是注入还是生产?
                                                                                    opt));
                    }
                    else
                    {
                        for (auto words: tmp)
                        {
                            if (std::find(words.begin(), words.end(), "SHUT") != words.end())
                                continue;

                            /// 生成一个WellOptParam对象
                            now += stod(words[1]);

                            string state = "OPEN";
                            string mode = words[2];
                            string max_rate = words[3];
                            double min_bhp;
                            if (words[4] == "BHP")
                                min_bhp = stoi(words[5]);
                            else
                                cout << "### ERROR: Prod Rate is missing in WCONINJE!" << endl;
                            WellOptParam opt(state, mode, max_rate, min_bhp);

                            paramWell.well_oper_list.push_back(ParamWell::WellOperation(now,
                                                                                        well_name,
                                                                                        "PROD", /// fff怎么知道是注入还是生产?
                                                                                        opt));
                        }
                    }
                }
                else if (words[0] == "LIMIT") /// fff
                {
                    string opt = words[1];
                    cout << "LIMIT: " << opt << endl;
                }
                else if (words[0] == "PERF") /// fff
                {
                    string opt = words[1];
                    cout << "PERF: " << opt << endl;
                }
                else if (words[0] == "STREAM")
                {
                    string opt = words[1];
                    cout << "STREAM: " << opt << endl;
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
        else
        {
            OCP_MESSAGE("Found not supported keywords in SCHEDULE section: " << words[0]);
            OCP_ABORT("Wrong keywords!");
        }
    }

    /// Establish Time and Well Operations
    vector<double> all_tsteps;
    for (int i=0; i<paramWell.well_oper_list.size(); ++i)
    {
        ParamWell::WellOperation opt = paramWell.well_oper_list[i];
        all_tsteps.push_back(opt.tstep);
    }
    /// Sort and unique
    std::sort(all_tsteps.begin(), all_tsteps.end());
    all_tsteps.erase(std::unique(all_tsteps.begin(), all_tsteps.end()), all_tsteps.end());
    /// construct criticalTime
    if (std::abs(all_tsteps[0] - 0.0) > 1.0E-10)
        paramWell.criticalTime.push_back(all_tsteps[0]);
    for (int i=1; i<all_tsteps.size(); ++i)
    {
        paramWell.criticalTime.push_back(all_tsteps[i]);
    }
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