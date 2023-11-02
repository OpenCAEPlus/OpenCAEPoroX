/*! \file    ParamReservoir.cpp
 *  \brief   ParamReservoir class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "ParamReservoir.hpp"


/// Find pointer to the specified table.
TableSet* ParamReservoir::FindPtrTable(const string& varName)
{
    TableSet* myPtr = nullptr;

    switch (Map_Str2Int(&varName[0], varName.size())) {
        case Map_Str2Int("SWFN", 4):
            myPtr = &SWFN_T;
            break;

        case Map_Str2Int("SWOF", 4):
            myPtr = &SWOF_T;
            break;

        case Map_Str2Int("SGFN", 4):
            myPtr = &SGFN_T;
            break;

        case Map_Str2Int("SGOF", 4):
            myPtr = &SGOF_T;
            break;

        case Map_Str2Int("SOF3", 4):
            myPtr = &SOF3_T;
            break;

        case Map_Str2Int("PBVD", 4):
            myPtr = &PBVD_T;
            break;

        case Map_Str2Int("PVCO", 4):
            myPtr = &PVCO_T;
            break;

        case Map_Str2Int("PVDO", 4):
            myPtr = &PVDO_T;
            break;

        case Map_Str2Int("PVCDO", 5):
            myPtr = &PVCDO_T;
            break;

        case Map_Str2Int("PVDG", 4):
            myPtr = &PVDG_T;
            break;

        case Map_Str2Int("PVTW", 4):
            myPtr = &PVTW_T;
            break;

        case Map_Str2Int("ZMFVD", 5):
            if (numCom == 0) OCP_ABORT("Number of Components has not been specified!");

            ZMFVD_T.colNum = numCom + 1;
            myPtr = &ZMFVD_T;
            break;

        case Map_Str2Int("TEMPVD", 6):
            myPtr = &TEMPVD_T;
            break;

        default:
            OCP_ABORT(varName + " is Inavailable!");
    }

    return myPtr;
}


Table2Set* ParamReservoir::FindPtrTable2(const string& varName)
{

    Table2Set* myPtr = nullptr;

    switch (Map_Str2Int(&varName[0], varName.size())) 
    {
        case Map_Str2Int("VISCTAB", 7):
            if (numCom == 0) OCP_ABORT("Number of Components has not been specified!");

            myPtr = &comsParam.viscTab;
            break;

        case Map_Str2Int("PVTH2O", 6):
            myPtr = &PVTH2O;
            break;

        case Map_Str2Int("PVTCO2", 6):
            myPtr = &PVTCO2;
            break;

        default:
            OCP_ABORT(varName + " is Inavailable!");
    }

    return myPtr;
}


/// Initialize tables and other reservoir parameters.
void ParamReservoir::Init()
{
    // Field default
    unitType = "FIELD";
    Psurf = FIELD_PRESSURE_STD;
    Tsurf = FIELD_TEMPERATURE_STD;


    InitTable();

    gravity.data.resize(3);
    gravity.data[0] = 45.5;   // oil
    gravity.data[1] = 1.0;    // pure water
    gravity.data[2] = 0.7773; // air

    density.data.resize(3);
    density.data[0] = 37.457;    // The density of oil at surface conditions: lb/ft3
    density.data[1] = 62.366416; // The density of water at surface conditions: lb/ft3
    density.data[2] = 0.062428;  // The density of gas at surface conditions: lb/ft3

    rsTemp = 60.0;
}

/// Initialize tables.
void ParamReservoir::InitTable()
{
    SWFN_T.name     = "SWFN";
    SWFN_T.colNum   = 3;
    SWOF_T.name     = "SWOF";
    SWOF_T.colNum   = 4;
    SGFN_T.name     = "SGFN";
    SGFN_T.colNum   = 3;
    SGOF_T.name     = "SGOF";
    SGOF_T.colNum   = 4;
    SOF3_T.name     = "SOF3";
    SOF3_T.colNum   = 3;
    PBVD_T.name     = "PBVD";
    PBVD_T.colNum   = 2;
    PVCO_T.name     = "PVCO";
    PVCO_T.colNum   = 6;
    PVDO_T.name     = "PVDO";
    PVDO_T.colNum   = 3;
    PVCDO_T.name    = "PVCDO";
    PVCDO_T.colNum  = 5;
    PVDG_T.name     = "PVDG";
    PVDG_T.colNum   = 3;
    PVTW_T.name     = "PVTW";
    PVTW_T.colNum   = 5;
    ZMFVD_T.name    = "ZMFVD";  // colnum equals numCom(hydrocarbon) + 1
    TEMPVD_T.name   = "TEMPVD"; // colnum equals 2
    TEMPVD_T.colNum = 2;

    comsParam.viscTab.name = "VISCTAB"; // colnum equals numCom(hydrocarbon) + 1
    PVTH2O.name = "PVTH2O";
    PVTCO2.name = "PVTCO2";
}

/// TODO: Add Doxygen
void ParamReservoir::InputCOMPS(ifstream& ifs)
{
    comps = OCP_TRUE;
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    numCom           = stoi(vbuf[0]);
    comsParam.numCom = numCom;
    comsParam.Init();

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        cout << endl << "COMPS" << endl;
        cout << numCom << endl;
    }
}


/// TODO: Add Doxygen
void ParamReservoir::InputRTEMP(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (vbuf[0] == "/") return;
    rsTemp = stod(vbuf[0]);

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT)
        cout << "RTEMP\n" << rsTemp << endl << endl;
}


/// TODO: Add Doxygen
void ParamReservoir::InputTABLE(ifstream& ifs, const string& tabName)
{
    TableSet* obj = FindPtrTable(tabName);

    const USI col = obj->colNum;
    vector<vector<OCP_DBL>> tmpTab(col);

    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        for (USI i = 0; i < col; i++) {
            tmpTab[i].push_back(stod(vbuf[i]));
        }

        if (vbuf.back() == "/") {
            obj->data.push_back(tmpTab);
            for (USI j = 0; j < col; j++) {
                tmpTab[j].clear();
            }
        }
    }
    if (!tmpTab[0].empty()) obj->data.push_back(tmpTab);

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT)
        obj->DisplayTable();
}


void ParamReservoir::InputTABLE2(ifstream& ifs, const string& tabName)
{
    Table2Set*      obj = FindPtrTable2(tabName);
    Table2          tmpTab(1);
    vector<string>  vbuf;

    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;
        if (vbuf[0][0] == '*') {
            // there is reference data
            tmpTab.refName = vbuf[0];
            tmpTab.refData.push_back(stod(vbuf[0]));
            tmpTab.data.resize(tmpTab.refData.size());
            continue;
        }
        auto& data = tmpTab.data.back();

        if (vbuf.back() == "/")        vbuf.pop_back();
        if (data.size() < vbuf.size()) data.resize(vbuf.size());
        for (USI i = 0; i < vbuf.size(); i++) {
            data[i].push_back(stod(vbuf[i]));
        }
    }
    tmpTab.SetColNum();
    obj->data.push_back(tmpTab);
}


/// Read data from the ROCK keyword.
void ParamReservoir::InputROCK(ifstream& ifs)
{
    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT)
        cout << "\n---------------------" << endl
            << "ROCK"
            << "\n---------------------" << endl;

    vector<string> vbuf;
    while (true) {
        ReadLine(ifs, vbuf);
        if (vbuf[0] == "/") break;

        RockParam rock;
        rock.type = vbuf[0];
        rock.Pref = stod(vbuf[1]);
        rock.cp1  = stod(vbuf[2]);

        if (rock.type == "LINEAR02") {
            if (vbuf.size() > 3 && vbuf[3] != "/") {
                rock.cp2 = stod(vbuf[3]);
            } else {
                rock.cp2 = rock.cp1;
            }
        }
        rockSet.push_back(rock);

        if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT)
            cout << "   " << rock.type << "   " << rock.Pref << "   " << rock.cp1 << "   "
                << rock.cp2 << endl;
    }
}

/// Read data from the ROCK keyword.
void ParamReservoir::InputROCKT(ifstream& ifs)
{
    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT)
        cout << "\n---------------------" << endl
            << "ROCKT"
            << "\n---------------------" << endl;

    RockParam      rock;
    vector<string> vbuf;
    while (true) {
        ReadLine(ifs, vbuf);
        if (vbuf[0] == "/") break;

        USI index = 0;
        USI len   = vbuf.size();
        while (index < len) {
            if (vbuf[index] == "*PORFORM") {
                rock.type = vbuf[index + 1];
            } else if (vbuf[index] == "*PRPOR") {
                rock.Pref = stod(vbuf[index + 1]);
            } else if (vbuf[index] == "*TRPOR") {
                rock.Tref = stod(vbuf[index + 1]);
            } else if (vbuf[index] == "*CPOR") {
                rock.cp1 = stod(vbuf[index + 1]);
            } else if (vbuf[index] == "*CTPOR") {
                rock.ct = stod(vbuf[index + 1]);
            } else if (vbuf[index] == "*CPTPOR") {
                rock.cpt = stod(vbuf[index + 1]);
            } else if (vbuf[index] == "*VOLCONST") {
                if (vbuf[index + 1] == "BULK") rock.ConstRock = OCP_FALSE;
            } else if (vbuf[index] == "*CP1") {
                rock.HCP1 = stod(vbuf[index + 1]);
            } else if (vbuf[index] == "*CP2") {
                rock.HCP2 = stod(vbuf[index + 1]);
            }
            index += 2;
        }
    }
    rockSet.push_back(rock);

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        cout << "*PORFORM   " << rock.type << endl;
        cout << "*PRPOR     " << rock.Pref << endl;
        cout << "*TRPOR     " << rock.Tref << endl;
        cout << "*CPOR      " << rock.cp1 << endl;
        cout << "*CTPOR     " << rock.ct << endl;
        cout << "*CPTPOR    " << rock.cpt << endl;
        cout << "*VOLCONST  " << (rock.ConstRock ? "ROCK" : "BULK") << endl;
        cout << "*CP1       " << rock.HCP1 << endl;
        cout << "*CP2       " << rock.HCP2 << endl;
    }
}

void ParamReservoir::InputHLOSS(ifstream& ifs)
{
    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT)
        cout << "\n---------------------" << endl
            << "HLOSSPROR"
            << "\n---------------------" << endl;

    hLoss.ifHLoss = OCP_TRUE;

    vector<string> vbuf;
    while (true) {
        ReadLine(ifs, vbuf);
        if (vbuf[0] == "/") break;

        USI index = 0;
        USI len   = vbuf.size();
        while (index < len) {
            if (vbuf[index] == "*OVERBUR") {
                hLoss.obUse = OCP_TRUE;
                hLoss.obC   = stod(vbuf[index + 1]);
                hLoss.obK   = stod(vbuf[index + 2]);
            } else if (vbuf[index] == "*UNDERBUR") {
                hLoss.ubUse = OCP_TRUE;
                hLoss.ubC   = stod(vbuf[index + 1]);
                hLoss.ubK   = stod(vbuf[index + 2]);
            }
            index += 3;
        }
    }
    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        cout << "*OVERBUR   " << hLoss.obC << "   " << hLoss.obK << endl;
        cout << "*UNDERBUR  " << hLoss.ubC << "   " << hLoss.ubK << endl;
    }
}


void ParamReservoir::InputBrooksCorey(ifstream& ifs)
{
    BrooksCoreyParam bc;
    vector<string>   vbuf;
    while (true) {
        ReadLine(ifs, vbuf);
        if (vbuf[0] == "/") break;

        USI     index = 0;
        USI     len   = vbuf.size();
        OCP_DBL val   = -1;
        while (index < len) {
            if (vbuf[index + 1] != "NAN")  val = stod(vbuf[index + 1]);
            else                           val = -1;
                
            if (vbuf[index] == "*SWIMM") {
                bc.sw_imm = val;
            }
            else if (vbuf[index] == "*SNIMM") {
                bc.sn_imm = val;
            }
            else if (vbuf[index] == "*PENTRY") {
                bc.Pentry = val;
            }
            else if (vbuf[index] == "*PCMAX") {
                bc.Pcmax = val;
            }
            else if (vbuf[index] == "*CWREPERM") {
                bc.Cw_kr = val;
            }
            else if (vbuf[index] == "*CNREPERM") {
                bc.Cn_kr = val;
            }
            else if (vbuf[index] == "*CPC") {
                bc.C_pc = val;
            }
            index += 2;
        }
    }
    BCparam.push_back(bc);
}

/// Read data from the MISCSTR keyword.
void ParamReservoir::InputMISCSTR(ifstream& ifs)
{
	miscstr.ifMiscible = OCP_TRUE;
	vector<string> vbuf;
	ReadLine(ifs, vbuf);
	if (vbuf[0] == "/") return;
	if (vbuf.back() == "/") vbuf.pop_back();
    DealDefault(vbuf);

	const USI len = vbuf.size();
    if (len > 0) {
        if (vbuf[0] != "DEFAULT")
            miscstr.surTenRef = stod(vbuf[0]);
    }
    if (len > 1) {
        if (vbuf[1] != "DEFAULT")
            miscstr.surTenEpt = stod(vbuf[1]);
    }
    if (len > 2) {
        if (vbuf[2] != "DEFAULT")
            miscstr.surTenPc = stod(vbuf[2]);
    }
    if (len > 3) {
        if (vbuf[3] != "DEFAULT")
            miscstr.surTenExp = stod(vbuf[3]);
    }

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        cout << "\n---------------------" << endl
            << "MISCSTR"
            << "\n---------------------" << endl;
        cout << miscstr.surTenRef << "   "
            << miscstr.surTenEpt << "   "
            << miscstr.surTenPc << "   "
            << miscstr.surTenExp << endl;
    }
}

/// Read data from the GRAVITY keyword.
void ParamReservoir::InputGRAVITY(ifstream& ifs)
{
    gravity.activity = OCP_TRUE;
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (vbuf[0] == "/") return;
    DealDefault(vbuf);
    if (vbuf.back() == "/") vbuf.pop_back();
    OCP_ASSERT(vbuf.size() == 3, "Wrong Keyword GRAVITY!");
    for (USI i = 0; i < 3; i++) {
        if (vbuf[i] != "DEFAULT") {
            gravity.data[i] = stod(vbuf[i]);
        }
    }

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        cout << "\n---------------------" << endl
            << "GRAVITY"
            << "\n---------------------" << endl;
        cout << "   " << gravity.data[0] << "  " << gravity.data[1] << "  "
            << gravity.data[2] << endl;
    }
}

/// Read data from the DENSITY keyword.
void ParamReservoir::InputDENSITY(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (vbuf[0] == "/") return;

    DealDefault(vbuf);
    if (vbuf.back() == "/") vbuf.pop_back();
    OCP_ASSERT(vbuf.size() == 3, "Wrong Keyword DENSITY!");
    for (USI i = 0; i < 3; i++) {
        if (vbuf[i] != "DEFAULT") {
            density.activity = OCP_TRUE;
            density.data[i]  = stod(vbuf[i]);
        }
    }

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        cout << "\n---------------------" << endl
            << "DENSITY"
            << "\n---------------------" << endl;
        cout << density.data[0] << "  " << density.data[1] << "  " << density.data[2]
            << endl;
    }
}

/// Read data from the THCONO, THCONG, THCONW, THCONR
void ParamReservoir::InputTHCON(ifstream& ifs, const string& keyword)
{
    ifThcon = OCP_TRUE;

    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (keyword == "THCONO") {
        thcono = stod(vbuf[0]);
    } else if (keyword == "THCONG") {
        thcong = stod(vbuf[0]);
    } else if (keyword == "THCONW") {
        thconw = stod(vbuf[0]);
    } else if (keyword == "THCONR") {
        thconr = stod(vbuf[0]);
    }

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        cout << "THCONO\n" << thcono << endl << endl;
        cout << "THCONG\n" << thcong << endl << endl;
        cout << "THCONW\n" << thconw << endl << endl;
        cout << "THCONR\n" << thconr << endl << endl;
    }
}

/// Read data from the EQUIL keyword.
void ParamReservoir::InputEQUIL(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (vbuf[0] == "/") return;

    EQUILParam tmpEQUIL;
    tmpEQUIL.data.resize(6, 0);
    DealDefault(vbuf);
    for (USI i = 0; i < 6; i++) {
        if (vbuf[i] != "DEFAULT") tmpEQUIL.data[i] = stod(vbuf[i]);
    }

    EQUIL.push_back(tmpEQUIL);

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        cout << "\n---------------------" << endl
            << "EQUIL"
            << "\n---------------------" << endl;
        cout << "   ";
        for (USI i = 0; i < 6; i++) cout << tmpEQUIL.data[i] << "  ";
        cout << endl;
    }
}

/// Read data from the TABDIMS keyword.
void ParamReservoir::InputTABDIMS(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);

    if (vbuf.size() < 3) {
        OCP_ABORT("Input the number of Saturation tables, PVT tables, and Rock tables "
                  "in turn!");
    }

    NTSFUN = stoi(vbuf[0]);
    comsParam.NTPVT = NTPVT  = stoi(vbuf[1]);
    NTROOC = stoi(vbuf[2]);

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        cout << "\n---------------------" << endl
            << "TABDIMS"
            << "\n---------------------" << endl;
        cout << "   " << NTSFUN << "   " << NTPVT << "   " << NTROOC << endl;
    }
}


void ParamReservoir::InputBoundary(ifstream& ifs)
{   
    vector<string> vbuf;
    while (OCP_TRUE) {
        ReadLine(ifs, vbuf);
        if (vbuf[0] == "BOUNDARYEND") break;
        while (OCP_TRUE) {
            if (vbuf[0] == "/") break;

            if (vbuf[0] == "*NAME") {   
                string name;
                for (USI i = 1; i < vbuf.size(); i++) {
                    if (vbuf[i] != "/") {
                        name += vbuf[i];
                    }
                }
                BDparam.push_back(BoundaryParam(name));
            }
            else if (vbuf[0] == "CONSTP") {
                BDparam.back().constP = OCP_TRUE;
                BDparam.back().P      = stod(vbuf[1]);
            }

            ReadLine(ifs, vbuf);
        }     
    } 
}


/// Check consistency of input parameters.
void ParamReservoir::CheckParam()
{
    CheckRock();
    CheckCPL();
    CheckCPG();
}


/// Check rock keyword.
void ParamReservoir::CheckRock()
{
    if (rockSet.size() != NTROOC) {
        OCP_ABORT("Wrong ROCK or ROCKT!");
    }
}


/// Check cpl1, cpl2, cpl3, cpl4
void ParamReservoir::CheckCPL()
{
    if (comsParam.cpl1.activity) {
        const USI l = comsParam.cpl1.data.size();
        const vector<OCP_DBL> tmp(comsParam.numCom, 0);
        for (USI i = 0; i < l; i++) {
            if (comsParam.cpl2.data.size() < l) {
                comsParam.cpl2.data.push_back(tmp);
            }
            if (comsParam.cpl3.data.size() < l) {
                comsParam.cpl3.data.push_back(tmp);
            }
            if (comsParam.cpl4.data.size() < l) {
                comsParam.cpl4.data.push_back(tmp);
            }
        }
    } 
}

/// Check cpg1, cpg2, cpg3, cpg4
void ParamReservoir::CheckCPG()
{
    if (comsParam.cpg1.activity) {
        const USI l = comsParam.cpg1.data.size();
        const vector<OCP_DBL> tmp(l, 0);
        for (USI i = 0; i < l; i++) {
            if (comsParam.cpg2.data.size() < i) {
                comsParam.cpg2.data.push_back(tmp);
            }
            if (comsParam.cpg3.data.size() < i) {
                comsParam.cpg3.data.push_back(tmp);
            }
            if (comsParam.cpg4.data.size() < i) {
                comsParam.cpg4.data.push_back(tmp);
            }
            if (comsParam.hvapr.data.size() < i) {
                comsParam.hvapr.data.push_back(tmp);
            }
            if (comsParam.hvr.data.size() < i) {
                comsParam.hvr.data.push_back(tmp);
            }
            if (comsParam.ev.data.size() < i) {
                comsParam.ev.data.push_back(tmp);
            }
        }
    }
}


/// TODO: Add Doxygen
void TableSet::DisplayTable() const
{
    cout << "\n---------------------" << endl
         << name << "\n---------------------" << endl;

    for (USI n = 0; n < data.size(); n++) {
        const USI len = data[n][0].size();
        for (USI i = 0; i < len; i++) {
            for (USI j = 0; j < colNum; j++) {
                cout << setw(10) << data[n][j][i];
            }
            cout << "\n";
        }
    }
}


void ComponentParam::Init()
{
    // Init LBC coefficient
    LBCcoef.resize(5);
    LBCcoef[0] = 0.1023;
    LBCcoef[1] = 0.023364;
    LBCcoef[2] = 0.058533;
    LBCcoef[3] = -0.040758;
    LBCcoef[4] = 0.0093324;
}

Type_A_r<vector<OCP_DBL>>* ComponentParam::FindPtr01(const string& varName)
{
    Type_A_r<vector<OCP_DBL>>* myPtr = nullptr;

    switch (Map_Str2Int(&varName[0], varName.size())) {
        case Map_Str2Int("TCRIT", 5):
            myPtr = &Tc;
            break;

        case Map_Str2Int("PCRIT", 5):
            myPtr = &Pc;
            break;

        case Map_Str2Int("VCRIT", 5):
            myPtr = &Vc;
            break;

        case Map_Str2Int("ZCRIT", 5):
            myPtr = &Zc;
            break;

        case Map_Str2Int("MW", 2):
            myPtr = &MW;
            break;

        case Map_Str2Int("ACF", 3):
            myPtr = &Acf;
            break;

        case Map_Str2Int("OMEGAA", 6):
            myPtr = &OmegaA;
            break;

        case Map_Str2Int("OMEGAB", 6):
            myPtr = &OmegaB;
            break;

        case Map_Str2Int("SSHIFT", 6):
            myPtr = &Vshift;
            break;

        case Map_Str2Int("PARACHOR", 8):
            myPtr = &parachor;
            break;

        case Map_Str2Int("VCRITVIS", 8):
            myPtr = &Vcvis;
            break;

        case Map_Str2Int("ZCRITVIS", 8):
            myPtr = &Zcvis;
            break;

        case Map_Str2Int("MOLDEN", 6):
            myPtr = &molden;
            break;

        case Map_Str2Int("CP", 2):
            myPtr = &cp;
            break;

        case Map_Str2Int("CT1", 3):
            myPtr = &ct1;
            break;

        case Map_Str2Int("CT2", 3):
            myPtr = &ct2;
            break;

        case Map_Str2Int("CPT", 3):
            myPtr = &cpt;
            break;

        case Map_Str2Int("CPL1", 4):
            myPtr = &cpl1;
            break;

        case Map_Str2Int("CPL2", 4):
            myPtr = &cpl2;
            break;

        case Map_Str2Int("CPL3", 4):
            myPtr = &cpl3;
            break;

        case Map_Str2Int("CPL4", 4):
            myPtr = &cpl4;
            break;

        case Map_Str2Int("CPG1", 4):
            myPtr = &cpg1;
            break;

        case Map_Str2Int("CPG2", 4):
            myPtr = &cpg2;
            break;

        case Map_Str2Int("CPG3", 4):
            myPtr = &cpg3;
            break;

        case Map_Str2Int("CPG4", 4):
            myPtr = &cpg4;
            break;

        case Map_Str2Int("HVAPR", 5):
            myPtr = &hvapr;
            break;

        case Map_Str2Int("HVR", 3):
            myPtr = &hvr;
            break;

        case Map_Str2Int("EV", 2):
            myPtr = &ev;
            break;

        case Map_Str2Int("AVSIC", 5):
            myPtr = &avisc;
            break;

        case Map_Str2Int("BVSIC", 5):
            myPtr = &bvisc;
            break;

        case Map_Str2Int("AVG", 3):
            myPtr = &avg;
            break;

        case Map_Str2Int("BVG", 3):
            myPtr = &bvg;
            break;
    }

    return myPtr;
}

void ComponentParam::InputRefPR(ifstream& ifs, const string& keyword)
{
    OCP_ASSERT(NTPVT > 0, "NTPVT has not been set!");

    vector<OCP_DBL>* objPtr = nullptr;
    objPtr                  = FindPtr02(keyword);
    if (objPtr == nullptr) {
        OCP_ABORT("Unknown keyword!");
    }

    vector<string> vbuf;
    while (OCP_TRUE) {
        ReadLine(ifs, vbuf);

        if (vbuf[0] == "/") break;

        for (auto& v : vbuf) {
            if (v != "/") objPtr->push_back(stod(v));
            if (objPtr->size() >= NTPVT) break;
        }
        if (objPtr->size() >= NTPVT) break;
    }

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        cout << keyword << endl;
        for (USI i = 0; i < NTPVT; i++) {
            cout << objPtr->at(i) << "   ";
        }
        cout << "\n/" << endl << endl;
    }
}

vector<OCP_DBL>* ComponentParam::FindPtr02(const string& varName)
{
    vector<OCP_DBL>* myPtr = nullptr;

    switch (Map_Str2Int(&varName[0], varName.size())) {
        case Map_Str2Int("PRSR", 4):
            myPtr = &Pref;
            break;

        case Map_Str2Int("TEMR", 4):
            myPtr = &Tref;
            break;
    }

    return myPtr;
}

void ComponentParam::InputCOMPONENTS(ifstream& ifs, const string& keyword)
{
    OCP_ASSERT((numCom > 0) && (NTPVT > 0), "number of components has not been input!");

    Type_A_r<vector<OCP_DBL>>* objPtr = nullptr;
    objPtr                            = FindPtr01(keyword);
    if (objPtr == nullptr) {
        OCP_ABORT("Unknown keyword!");
    }
    objPtr->activity = OCP_TRUE;

    vector<string>  vbuf;
    vector<OCP_DBL> tmp;
    USI             nReg = 0;

    while (OCP_TRUE) {
        ReadLine(ifs, vbuf);
        if (vbuf[0] == "/") {
            nReg++;
            tmp.resize(numCom);
            objPtr->data.push_back(tmp);
            if (nReg >= NTPVT) break;
            tmp.clear();
            continue;
        }
        for (auto& v : vbuf) {
            if (v != "/") {
                tmp.push_back(stod(v));
            }
        }
        if (vbuf.back() == "/") {
            nReg++;
            tmp.resize(numCom);
            objPtr->data.push_back(tmp);
            tmp.clear();
            if (nReg >= NTPVT) break;
        }
    }

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        cout << keyword << endl;
        for (USI i = 0; i < NTPVT; i++) {
            for (auto& v : objPtr->data[i]) {
                cout << v << endl;
            }
            cout << "/" << endl;
        }
        cout << endl;
    }
}

void ComponentParam::InputCNAMES(ifstream& ifs)
{
    OCP_ASSERT(numCom > 0, "numCom has not been set!");

    vector<string> vbuf;
    while (OCP_TRUE) {
        ReadLine(ifs, vbuf);
        if (vbuf[0] == "/") {
            break;
        }
        for (auto& v : vbuf) {
            if (v != "/") Cname.push_back(v);
        }
        if (vbuf.back() == "/") break;
    }

    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        cout << "CNAMES" << endl;
        for (USI i = 0; i < numCom; i++) {
            cout << Cname[i] << "   ";
        }
        cout << endl << endl;
    }
}

void ComponentParam::InputLBCCOEF(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    DealDefault(vbuf);
    for (USI i = 0; i < 5; i++) {
        if (vbuf[i] != "DEFAULT") LBCcoef[i] = stod(vbuf[i]);
    }


    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        cout << "LBCCOEF" << endl;
        for (USI i = 0; i < 5; i++) {
            cout << LBCcoef[i] << "   ";
        }
        cout << endl << endl;
    }  
}

/// Input Binary Interaction Coefficients Matrix
void ComponentParam::InputBIC(ifstream& ifs)
{
    OCP_ASSERT((numCom > 0) && (NTPVT > 0), "numCom or NTPVT has not been set!");

    BIC.resize(NTPVT);

    vector<string> vbuf;
    USI            nReg = 0;
    while (OCP_TRUE) {
        ReadLine(ifs, vbuf);
        if (vbuf[0] == "/") {
            nReg++;
            if (nReg >= NTPVT) break;
            continue;
        }
        for (auto& v : vbuf) {
            if (v != "/") {
                BIC[nReg].push_back(stod(v));
                if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT)
                    cout << setw(10) << BIC[nReg].back();
            }
        }
        if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT)
            cout << endl;
        if (vbuf.back() == "/") {
            nReg++;
            if (nReg >= NTPVT) break;
        }
    }
}


/// TODO: Add Doxygen
void ComponentParam::InputSSMSTA(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    int len = vbuf.size();
    for (int i = 0; i < len; i++) {
        SSMparamSTA.push_back(vbuf[i]);
    }
    
    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        for (int i = 0; i < len; i++) {
            cout << SSMparamSTA[i] << "   ";
        }
        cout << endl << endl;
    }
}

/// TODO: Add Doxygen
void ComponentParam::InputNRSTA(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    for (int i = 0; i < 2; i++) {
        NRparamSTA.push_back(vbuf[i]);
    }
    
    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        for (int i = 0; i < 2; i++) {
            cout << NRparamSTA[i] << "   ";
        }
        cout << endl << endl;
    }
}

/// TODO: Add Doxygen
void ComponentParam::InputSSMSP(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    for (USI i = 0; i < 2; i++) {
        SSMparamSP.push_back(vbuf[i]);
    }
    
    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        for (USI i = 0; i < 2; i++) {
            cout << SSMparamSP[i] << "   ";
        }
        cout << endl << endl;
    }
}

/// TODO: Add Doxygen
void ComponentParam::InputNRSP(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    for (USI i = 0; i < 2; i++) {
        NRparamSP.push_back(vbuf[i]);
    }
    
    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        for (USI i = 0; i < 2; i++) {
            cout << NRparamSP[i] << "   ";
        }
        cout << endl << endl;
    }
}

/// TODO: Add Doxygen
void ComponentParam::InputRR(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    for (USI i = 0; i < 2; i++) {
        RRparam.push_back(vbuf[i]);
    }
    
    if (CURRENT_RANK == MASTER_PROCESS && PRINTINPUT) {
        for (USI i = 0; i < 2; i++) {
            cout << RRparam[i] << "   ";
        }
        cout << endl << endl;
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/09/2022      Update output and Doxygen            */
/*----------------------------------------------------------------------------*/