/*! \file    ParamWell.cpp
 *  \brief   ParamWell class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */
#include <algorithm>

#include "ParamWell.hpp"

WellOptParam::WellOptParam(string intype, vector<string>& vbuf)
{
    type = intype;
    if (type == "INJ") {
        fluidType = vbuf[1];
        state     = vbuf[2];
        mode   = vbuf[3];
        if (vbuf[4] == "DEFAULT") {
            if (mode == "BHP")
                maxRate = 1E10;
            else {
                cout << "### ERROR: Inj Rate is missing in WCONINJE!" << endl;
            }
        } else {
            maxRate = stod(vbuf[4]);
        }
        maxBHP = stod(vbuf[5]);
    } else if (type == "PROD") {
        state   = vbuf[1];
        mode = vbuf[2];
        if (vbuf[3] == "DEFAULT") {
            if (mode == "BHP")
                maxRate = 1E10;
            else {
                cout << "### ERROR: Prod Rate is missing in WCONINJE!" << endl;
            }
        } else {
            maxRate = stod(vbuf[3]);
        }
        minBHP = stod(vbuf[4]);
    } else {
        OCP_ABORT("Wrong Well Type!");
    }
}


WellParam::WellParam(vector<string>& info)
{
    gridType = GridType::structured;
    name = info[0];
    if (info[1] != "DEFAULT") group = info[1];
    I = stoi(info[2]);
    J = stoi(info[3]);
    if (info[4] != "DEFAULT") depth = stod(info[4]);
}

WellParam::WellParam(const std::string& name_)
{
    gridType = GridType::structured;
    name = name_;
}


WellParam::WellParam(vector<string>& info, const string& unstructured)
{
    gridType = GridType::unstructured;
    name = info[0];
    if (info[1] != "DEFAULT") group = info[1];
    X = stod(info[2]);
    Y = stod(info[3]);
    Z = stod(info[4]);
}


void WellParam::InputCOMPDAT(vector<string>& vbuf)
{
    if (gridType == GridType::structured) {
        InputCOMPDATS(vbuf);
    }
    else if (gridType == GridType::unstructured) {
        InputCOMPDATUS(vbuf);
    }
    else {
        OCP_ABORT("INAVAILABLE GRID TYPE!");
    }
}


void WellParam::InputCOMPDATS(vector<string>& vbuf)
{
    const USI k1 = stoi(vbuf[3]);
    const USI k2 = stoi(vbuf[4]);

    for (USI k = k1; k <= k2; k++) {
        if (vbuf[1] == "DEFAULT" || vbuf[2] == "DEFAULT") {
            I_perf.push_back(I);
            J_perf.push_back(J);
        }
        else {
            I_perf.push_back(stoi(vbuf[1]));
            J_perf.push_back(stoi(vbuf[2]));
        }
        K_perf.push_back(k);

        if (vbuf[5] != "DEFAULT")
            WI.push_back(stod(vbuf[5]));
        else
            WI.push_back(-1.0);

        if (vbuf[6] != "DEFAULT")
            diameter.push_back(stod(vbuf[6]));
        else
            diameter.push_back(1.0);

        if (vbuf[7] != "DEFAULT")
            kh.push_back(stod(vbuf[7]));
        else
            kh.push_back(-1.0);

        if (vbuf[8] != "DEFAULT")
            skinFactor.push_back(stod(vbuf[8]));
        else
            skinFactor.push_back(0.0);

        if (vbuf[9] != "DEFAULT")
            direction.push_back(vbuf[9]);
        else
            direction.push_back("z");
    }
}


void WellParam::InputCOMPDATUS(vector<string>& vbuf)
{
    if (vbuf[1] == "DEFAULT")  X_perf.push_back(X);
    else                       X_perf.push_back(stod(vbuf[1]));
    if (vbuf[2] == "DEFAULT")  Y_perf.push_back(Y);
    else                       Y_perf.push_back(stod(vbuf[2]));
    if (vbuf[3] == "DEFAULT")  Z_perf.push_back(Z);
    else                       Z_perf.push_back(stod(vbuf[3]));

    if (vbuf[4] != "DEFAULT")  WI.push_back(stod(vbuf[4]));
    else                       WI.push_back(-1.0);

    if (vbuf[5] != "DEFAULT")  diameter.push_back(stod(vbuf[5]));
    else                       diameter.push_back(1.0);

    if (vbuf[6] != "DEFAULT")  kh.push_back(stod(vbuf[6]));
    else                       kh.push_back(-1.0);

    if (vbuf[7] != "DEFAULT")  skinFactor.push_back(stod(vbuf[7]));
    else                       skinFactor.push_back(0.0);

    direction.push_back("usg");
}


Solvent::Solvent(const vector<string>& vbuf)
{
    name    = vbuf[0];
    USI len = vbuf.size();
    for (USI i = 0; i < len - 1; i++) {
        if (vbuf[i + 1] == "/") break;
        comRatio.push_back(stod(vbuf[i + 1]));
    }
}

int ParamWell::FindWellByName(const string& name)
{
    for (int i=0; i<well.size(); ++i)
    {
        if (well[i].name == name)
            return i;
    }

    OCP_ABORT("Not found well!");
}


void ParamWell::Init() 
{ 
    // Field default
    Psurf = FIELD_PRESSURE_STD;
    Tsurf = FIELD_TEMPERATURE_STD;
    InitTime(); 
};

void ParamWell::InputTEMPLATE(ifstream &ifs)
{
    vector<string> vbuf;
    int index = 0;
    do {
        ReadLine(ifs, vbuf);

        int size = vbuf.size();
        for (int i=0; i<size-1; ++i)
        {
            templates.insert(std::pair<std::string, int>(vbuf[i], index));
            index++;
        }

        if (vbuf[size-1] != "/")
        {
            templates.insert(std::pair<std::string, int> (vbuf[size-1], index));
        }
        else
            break;

    } while (true);
}

void ParamWell::InputWELSPECS(ifstream& ifs)
{
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        DealDefault(vbuf);
        const USI len = vbuf.size();

        if (vbuf[len - 1] == "COORDINATE" || vbuf[len - 2] == "COORDINATE") {
            well.push_back(WellParam(vbuf, "unstructrued"));
        }
        else {
            well.push_back(WellParam(vbuf));
        }
    }
}

void ParamWell::InputWELSPECS(ifstream& ifs, int type)
{
    assert(type == 1);
    vector<string> vbuf;

    ReadLine(ifs, vbuf);
    assert(vbuf[0] == "NAME");
    int well_index = well.size();
    well.push_back(WellParam(vbuf[1]));

    bool first_perf = true;
    int i_, j_, k_;
    std::map<std::string, int>::iterator iter;
    int pos = ifs.tellg();
    while (ReadLine(ifs, vbuf))
    {
        if (isRegularString(vbuf[0]))
        {
            ifs.seekg(pos);
            break;
        }

        pos = ifs.tellg();

        DealDefault(vbuf);
        const USI len = vbuf.size();

        /// Obtain perforations data
        int i = stoi(vbuf[templates["I"]]);
        int j = stoi(vbuf[templates["J"]]);
        if (first_perf)
        {
            i_ = i;
            j_ = j;
            first_perf = false;
        }
        else
        { /// 一口井只有一个井头fff: (i,j)定义一个井头
//            assert(i == i_);
//            assert(j == j_);
        }
        int k1 = stoi(vbuf[templates["K1"]]);
        int k2 = stoi(vbuf[templates["K2"]]);
        double diam = stod(vbuf[templates["RW"]]) * 2;

        for (int k=k1; k<=k2; ++k)
        {
            well[well_index].I_perf.push_back(i);
            well[well_index].J_perf.push_back(j);
            well[well_index].K_perf.push_back(k);
            well[well_index].diameter.push_back(diam);
            well[well_index].skinFactor.push_back(0.0);
            well[well_index].kh.push_back(-1.0);
            well[well_index].WI.push_back(-1.0); // fff
            well[well_index].direction.push_back("z");
        }
    }
}


void ParamWell::InputCOMPDAT(ifstream& ifs)
{
    USI            num = well.size();
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        DealDefault(vbuf);
        string            src   = vbuf[0];
        string::size_type pos   = src.find("*");
        OCP_BOOL          match = (pos != string::npos);
        if (match) {
            src.erase(pos);
        }
        OCP_BOOL tmp = OCP_FALSE;

        for (USI w = 0; w < num; w++) {
            if (match)
                tmp = (well[w].name.substr(0, pos) == src);
            else
                tmp = (well[w].name == src);

            if (tmp) {
                well[w].InputCOMPDAT(vbuf);
            }
        }
    }
}

void ParamWell::InputWCONINJE(ifstream& ifs)
{
    assert(criticalTime.size() > 0);

    const USI      d   = criticalTime.size() - 1;
    const USI      num = well.size();
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        DealDefault(vbuf);
        string            src        = vbuf[0];
        string::size_type pos        = src.find("*");
        const OCP_BOOL    fuzzyMatch = (pos != string::npos);
        if (fuzzyMatch) {
            src.erase(pos);
        }

        if (fuzzyMatch) {
            for (USI w = 0; w < num; w++) {
                if (well[w].name.find(src) != string::npos) {
                    well[w].optParam.push_back(WellOptPair(d, "INJ", vbuf));
                }
            }
        } else {
            for (USI w = 0; w < num; w++) {
                if (well[w].name == src) {
                    well[w].optParam.push_back(WellOptPair(d, "INJ", vbuf));
                }
            }
        }
    }
    // cout << "WCONINJE" << endl;
}

void ParamWell::InputWCONPROD(ifstream& ifs)
{
    assert(criticalTime.size() > 0);

    const USI      d   = criticalTime.size() - 1;
    const USI      num = well.size();
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        DealDefault(vbuf);
        string            src        = vbuf[0];
        string::size_type pos        = src.find("*");
        const OCP_BOOL    fuzzyMatch = (pos != string::npos);
        if (fuzzyMatch) {
            src.erase(pos);
        }

        if (fuzzyMatch) {
            for (USI w = 0; w < num; w++)
                if (well[w].name.find(src) != string::npos)
                    well[w].optParam.push_back(WellOptPair(d, "PROD", vbuf));
        } else {
            for (USI w = 0; w < num; w++)
                if (well[w].name == src)
                    well[w].optParam.push_back(WellOptPair(d, "PROD", vbuf));
        }
    }
    // cout << "WCONPROD" << endl;
}

void ParamWell::InputTSTEP(ifstream& ifs)
{
    assert(criticalTime.size() > 0);

    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        DealDefault(vbuf);
        OCP_INT len = vbuf.size();
        for (OCP_INT i = 0; i < len - 1; i++) {
            OCP_DBL t = criticalTime.back() + stod(vbuf[i]);
            criticalTime.push_back(t);
        }
        if (vbuf.back() != "/") {
            OCP_DBL t = criticalTime.back() + stod(vbuf.back());
            criticalTime.push_back(t);
        }
    }
}

void ParamWell::InputTIME(const std::vector<std::string>& vbuf)
{
    assert(criticalTime.size() > 0);
    assert(vbuf[0] == "TIME");

    if (WhichDateFormat(vbuf[1]) == 0)
    {
        if (is_first_time)
        {
            first_time = vbuf[1];
            is_first_time = OCP_FALSE;

            criticalTime.push_back(0.0);
        }
        else
        {
            int days = NumDaysBetweenDates(first_time, vbuf[1]);
            criticalTime.push_back(days);
        }
    }
    else
        OCP_ABORT("Not support double time format now!");
}

void ParamWell::InputWELLOperations(vector<string> str, ifstream &ifs)
{
    string well_name = str[1];
    int well_index = FindWellByName(well_name);
    string well_state = "OPEN";
    string mode = "";

    if (str.size() > 2)
        OCP_ABORT("Not support now!");

//    if (well_name == "YA002-2")
//        cout << "Goog" << endl;

    const USI time_index = criticalTime.size() - 1;

    vector<string> vbuf;
    int pos = ifs.tellg();
    while (ReadLine(ifs, vbuf) && (vbuf[0] != "WELL" || vbuf[0] != "TIME"))
    {
        if (vbuf[0] == "WELL" || vbuf[0] == "TIME")
        {
            ifs.seekg(pos);
            break;
        }
        pos = ifs.tellg(); // save position

        if (std::find(vbuf.begin(), vbuf.end(), "SHUT") != vbuf.end())
        {
            well_state = "SHUT";

        }
        else if (std::find(vbuf.begin(), vbuf.end(), "GRAT") != vbuf.end())
        {
            int index = std::distance(vbuf.begin(), std::find(vbuf.begin(), vbuf.end(), "GRAT"));
            // 生产井

        }
        else if (vbuf[0] == "PERF")
        {
            continue;
        }
        else if (vbuf[0] == "/")
        {
            break;
        }
        else
        {
            string msg = "Not support WELL opertion: " + vbuf[0];
            OCP_ABORT(msg);
        }

    }
}

void ParamWell::InputWELTARG(ifstream& ifs)
{
    assert(criticalTime.size() > 0);

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT)
        cout << "\n---------------------" << endl
            << "WELTARG"
            << "\n---------------------" << endl;

    const USI      d   = criticalTime.size() - 1;
    const USI      num = well.size();
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
            for (const auto& s : vbuf) {
                cout << s << "   ";
            }
            cout << endl;
        }

        string            src        = vbuf[0];
        string::size_type pos        = src.find("*");
        const OCP_BOOL    fuzzyMatch = (pos != string::npos);
        if (fuzzyMatch) {
            src.erase(pos);
        }
        OCP_BOOL succMatch = OCP_FALSE;

        for (USI w = 0; w < num; w++) {
            if (fuzzyMatch)
                succMatch = (well[w].name.find(src) != string::npos);
            else
                succMatch = (well[w].name == src);

            if (succMatch) {
                if (well[w].optParam.size() == 0) {
                    OCP_ABORT("No Well Control Defined in Advance!");
                }
                WellOptPair tar = well[w].optParam.back();
                tar.d           = d;
                tar.opt.mode = vbuf[1];
                OCP_DBL val     = stod(vbuf[2]);
                if (vbuf[1] == "BHP") {
                    if (tar.opt.type == "INJ")
                        tar.opt.maxBHP = val;
                    else
                        tar.opt.minBHP = val;
                } else {
                    tar.opt.maxRate = val;
                }
                well[w].optParam.push_back(tar);
            }
        }
    }
}

void ParamWell::InputWTEMP(ifstream& ifs)
{
    assert(criticalTime.size() > 0);

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT)
        cout << "\n---------------------" << endl
            << "WTEMP"
            << "\n---------------------" << endl;

    const USI      d       = criticalTime.size() - 1;
    const USI      numWell = well.size();
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
            for (const auto& s : vbuf) {
                cout << s << "   ";
            }
            cout << endl;
        }

        string            src        = vbuf[0];
        string::size_type pos        = src.find("*");
        const OCP_BOOL    fuzzyMatch = (pos != string::npos);
        if (fuzzyMatch) {
            src.erase(pos);
        }
        OCP_BOOL succMatch = OCP_FALSE;

        for (USI w = 0; w < numWell; w++) {
            if (fuzzyMatch)
                succMatch = (well[w].name.find(src) != string::npos);
            else
                succMatch = (well[w].name == src);

            if (succMatch) {
                if (well[w].optParam.size() == 0) {
                    OCP_ABORT("No Well Control Defined in Advance!");
                }
                WellOptPair tar = well[w].optParam.back();
                tar.d           = d;
                tar.opt.injTemp = stod(vbuf[1]);
                well[w].optParam.push_back(tar);
            }
        }
    }
}

void ParamWell::InputUNWEIGHT(ifstream& ifs)
{
    assert(criticalTime.size() > 0);

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        cout << "\n---------------------" << endl
            << "UNWEIGHT"
            << "\n---------------------" << endl;
    }

    const USI      num = well.size();
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
            for (const auto& s : vbuf) {
                cout << s << "   ";
            }
            cout << endl;
        }

        DealDefault(vbuf);
        string            src        = vbuf[0];
        string::size_type pos        = src.find("*");
        const OCP_BOOL    fuzzyMatch = (pos != string::npos);
        if (fuzzyMatch) {
            src.erase(pos);
        }

        if (fuzzyMatch) {
            for (USI w = 0; w < num; w++)
                if (well[w].name.find(src) != string::npos)
                    well[w].ifUseUnweight = OCP_TRUE;
        } else {
            for (USI w = 0; w < num; w++)
                if (well[w].name == src) well[w].ifUseUnweight = OCP_TRUE;
        }
    }
}

void ParamWell::InputWELLSTRE(ifstream& ifs)
{
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;
        solSet.push_back(Solvent(vbuf));
    }
    // cout << "WELLSTRE" << endl;
}

void ParamWell::InputPSURF(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);

    Psurf = stod(vbuf[0]);

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        cout << "\n---------------------" << endl
            << "PSURF"
            << "\n---------------------" << endl;
        cout << Psurf << endl;
    }
}

void ParamWell::InputTSURF(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);

    Tsurf = stod(vbuf[0]);


    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        cout << "\n---------------------" << endl
            << "TSURF"
            << "\n---------------------" << endl;
        cout << Tsurf << endl;
    }
}


void ParamWell::InputWELINITP(ifstream& ifs)
{
    const USI      num = well.size();
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        DealDefault(vbuf);
        string            src = vbuf[0];
        string::size_type pos = src.find("*");
        const OCP_BOOL    fuzzyMatch = (pos != string::npos);
        if (fuzzyMatch) {
            src.erase(pos);
        }

        if (fuzzyMatch) {
            for (USI w = 0; w < num; w++)
                if (well[w].name.find(src) != string::npos)
                    well[w].initP = stod(vbuf[1]);
        }
        else {
            for (USI w = 0; w < num; w++)
                if (well[w].name == src)
                    well[w].initP = stod(vbuf[1]);
        }
    }
}


// check
void ParamWell::CheckParam() const
{
    CheckPerf();
}


void ParamWell::CheckPerf() const
{
    USI wellnum = well.size();
    USI perfnum;
    for (USI w = 0; w < wellnum; w++) {
        perfnum = well[w].GetPerfNum();
        if (well[w].diameter.size() != perfnum) {
            OCP_ABORT("Wrong perforation diameter size!");
        }
        if (well[w].WI.size() != perfnum) {
            OCP_ABORT("Wrong perforation WI size!");
        }
        if (well[w].kh.size() != perfnum) {
            OCP_ABORT("Wrong perforation kh size!");
        }
        if (well[w].skinFactor.size() != perfnum) {
            OCP_ABORT("Wrong perforation skinFactor size!");
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
/*  Chensong Zhang      Oct/27/2021      Unify error messages                 */
/*----------------------------------------------------------------------------*/