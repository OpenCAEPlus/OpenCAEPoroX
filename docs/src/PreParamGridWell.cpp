/*! \file    PreParamGridWell.cpp
 *  \brief   PreParamGridWell class definition
 *  \author  Shizhe Li
 *  \date    Feb/15/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "PreParamGridWell.hpp"


void PreParamGridWell::InputFile(const string& myFile, const string& myWorkdir)
{
    OCP_INFO("Input Grid File -- begin");

    workdir = myWorkdir;
    Input(myFile);
    CheckInput();
    PostProcessInput();

    OCP_INFO("Input Grid File -- end");
}


void PreParamGridWell::Input(const string& myFilename)
{
    ifstream ifs(workdir + myFilename, ios::in);
    if (!ifs) {
        OCP_MESSAGE("Trying to open file: " << (workdir + myFilename));
        OCP_ABORT("Failed to open the input file!");
    }

    while (!ifs.eof()) {
        vector<string> vbuf;
        if (!ReadLine(ifs, vbuf)) break;
        string keyword = vbuf[0];

        switch (Map_Str2Int(&keyword[0], keyword.size())) 
        {
            case Map_Str2Int("MODEL", 5):
                InputMODEL(ifs);
                break;

            case Map_Str2Int("DUALPORO", 8):
                InputDUALPORO();
                break;

            case Map_Str2Int("DPGRID", 6):
                InputDPGRID();
                break;

            case Map_Str2Int("DIMENS", 6):
                InputDIMENS(ifs);
                break;

            case Map_Str2Int("EQUALS", 6):
                InputEQUALS(ifs);
                break;

            case Map_Str2Int("COPY", 4):
                InputCOPY(ifs);
                break;

            case Map_Str2Int("MULTIPLY", 8):
                InputMULTIPLY(ifs);
                break;

            case Map_Str2Int("DX", 2):
            case Map_Str2Int("DY", 2):
            case Map_Str2Int("DZ", 2):
            case Map_Str2Int("TOPS", 4):
            case Map_Str2Int("COORD", 5):
            case Map_Str2Int("ZCORN", 5):
            case Map_Str2Int("NTG", 3):
            case Map_Str2Int("PORO", 4): 
            case Map_Str2Int("PERMX", 5):
            case Map_Str2Int("PERMY", 5):
            case Map_Str2Int("PERMZ", 5):
            case Map_Str2Int("ACTNUM", 6):
            case Map_Str2Int("SATNUM", 6):
            case Map_Str2Int("PVTNUM", 6):           
            case Map_Str2Int("ROCKNUM", 7):
            case Map_Str2Int("SWAT", 4):
            case Map_Str2Int("SWATINIT", 8):
            case Map_Str2Int("SIGMAV", 6):
            case Map_Str2Int("MULTZ", 5):
            case Map_Str2Int("DZMTRXV", 7):
                InputGrid(ifs, keyword);
                break;

            case Map_Str2Int("INCLUDE", 7):
                InputINCLUDE(ifs);
                break;
#ifdef WITH_GMSH
            case Map_Str2Int("GMSH", 4):
                InputGMSH(ifs);
                break;

            case Map_Str2Int("GMSHPRO", 7):
                InputGMSHPRO(ifs);
                break;
#endif

            case Map_Str2Int("WELSPECS", 8):
                InputWELSPECS(ifs);
                break;

            case Map_Str2Int("COMPDAT", 7):
                InputCOMPDAT(ifs);
                break;

            case Map_Str2Int("VTKSCHED", 8):
                ifUseVtk = OCP_TRUE;
                break;

            default: // skip non-keywords
                break;
        }
    }

    ifs.close();
}


void PreParamGridWell::CheckInput()
{
    cout << endl << "-------------------------------------" << endl;
    cout << "Check Grid param ... begin" << endl;

    if (model == OCPModel::none)           OCP_ABORT("WRONG MODEL!");


    if (gridType == GridType::corner) {
        if (nx == 0 || ny == 0 || nz == 0) OCP_ABORT("WRONG DIMENS!");
        if (poro.size() != numGrid)        OCP_ABORT("WRONG PORO!");
        if (zcorn.empty())                 OCP_ABORT("WRONG ZCORN!");
        if (coord.empty())                 OCP_ABORT("WRONG ZCORN!");
    }
    else if (gridType == GridType::orthogonal) {
        if (nx == 0 || ny == 0 || nz == 0) OCP_ABORT("WRONG DIMENS!");
        if (poro.size() != numGrid)        OCP_ABORT("WRONG PORO!");
        if (dx.size() != numGrid)          OCP_ABORT("WRONG DX!");
        if (dy.size() != numGrid)          OCP_ABORT("WRONG DY!");
        if (dz.size() != numGrid)          OCP_ABORT("WRONG DZ!");
        if (tops.size() != nx * ny)        OCP_ABORT("WRONG TOPS!");
    }
#ifdef WITH_GMSH
    else if (gridType == GridType::gmsh) {
        if (gmshGrid.edges.empty())       OCP_ABORT("WRONG GMSH!");
        if (gmshGrid.elements.empty())    OCP_ABORT("WRONG GMSH!");
    }
#endif
    else                                  OCP_ABORT("WRONG Grid Type!");

    cout << "Check Grid param ... done";
    cout << endl << "-------------------------------------" << endl;
}



void PreParamGridWell::PostProcessInput()
{
    if (ntg.size() != numGrid) {
        OCP_WARNING("NTG will be set to 1 !");
        vector<OCP_DBL>().swap(ntg);
    }
    else {
        for (OCP_ULL n = 0; n < numGrid; n++) {
            poro[n] *= ntg[n];
        }
    }

    if (actGC.ACTNUM.size() != numGrid) {
        OCP_WARNING("ACTNUM will be set to 1 !");
        vector<USI>().swap(actGC.ACTNUM);
        actGC.allAct = OCP_TRUE;
        
    }
    if (!sigma.empty())  sigma.resize(numGrid, 0);
    if (!dzMtrx.empty()) dzMtrx.resize(numGrid, 0);
}


void PreParamGridWell::InputMODEL(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (vbuf[0] == "THERMAL") {
        model = OCPModel::thermal;
    }
    else if (vbuf[0] == "ISOTHERMAL")
    {
        model = OCPModel::isothermal;
    }
    else {
        OCP_ABORT("WRONG MODEL in keyword MODEL!");
    }

    if (PRINTINPUT) {
        cout << "MODEL" << endl;
        cout << vbuf[0] << endl << endl;
    }
}


void PreParamGridWell::InputDIMENS(ifstream& ifs)
{

    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    nx = stoi(vbuf[0]);
    ny = stoi(vbuf[1]);
    nz = stoi(vbuf[2]);
    numGridM = nx * ny * nz;

    if (DUALPORO) numGridF = numGridM;

    numGrid = numGridM + numGridF;

    if (PRINTINPUT) {
        cout << "DIMENS" << endl;
        cout << setw(6) << nx << setw(6) << ny << setw(6) << nz << endl << endl;
    }
}


void PreParamGridWell::InputEQUALS(ifstream& ifs)
{

    if (PRINTINPUT) {
        cout << "\n---------------------" << endl
            << "EQUALS"
            << "\n---------------------" << endl;
    }

    vector<USI>    index(6, 0);
    vector<string> vbuf;

    USI nzTmp = nz;
    if (DUALPORO) nzTmp *= 2;

    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        index[0] = 0, index[1] = nx - 1;
        index[2] = 0, index[3] = ny - 1;
        index[4] = 0, index[5] = nzTmp - 1;

        const string  objName = vbuf[0];
        const OCP_DBL val     = stod(vbuf[1]);

        DealDefault(vbuf);

        for (USI n = 2; n < 8; n++) {
            if (vbuf[n] != "DEFAULT") index[n - 2] = stoi(vbuf[n]) - 1;
        }

        if (index[0] < 0 || index[2] < 0 || index[4] < 0 || index[1] > nx - 1 ||
            index[3] > ny - 1 || index[5] > nzTmp - 1) {
            OCP_ABORT("WRONG Range in " + objName + " in EQUALS!");
        }

        if (PRINTINPUT) {
            cout << setw(8) << vbuf[0] << setw(16) << vbuf[1];
            for (USI i = 0; i < 6; i++) {
                cout << setw(6) << index[i] + 1;
            }
            cout << endl;
        }

        
        {
            auto objPtr = FindPtr(objName, (OCP_DBL)0);
            if (objPtr != nullptr) {
                objPtr->resize(objPtr->capacity());
                if (objName == "TOPS") {
                    index[4] = index[5] = 0;
                }
                setVal(*objPtr, val, index);
                continue;
            }
        }

        {
            auto objPtr = FindPtr(objName, (USI)0);
            if (objPtr != nullptr) {
                objPtr->resize(objPtr->capacity());
                setVal(*objPtr, (USI)val, index);
                continue;
            }
        }


        OCP_ABORT("WRONG Item " + objName + " in EQUALS!");
    }

    if (PRINTINPUT) {
        cout << "/" << endl;
    }
}


void PreParamGridWell::InputCOPY(ifstream& ifs)
{
    if (PRINTINPUT) {
        cout << "\n---------------------" << endl
            << "COPY"
            << "\n---------------------" << endl;
    }

    vector<string> vbuf;
    vector<USI>    index(6, 0);

    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        index[0] = 0, index[1] = nx - 1;
        index[2] = 0, index[3] = ny - 1;
        index[4] = 0, index[5] = nz - 1;

        string srcName = vbuf[0];
        string objName = vbuf[1];
        DealDefault(vbuf);
        for (USI n = 2; n < 8; n++) {
            if (vbuf[n] != "DEFAULT") index[n - 2] = stoi(vbuf[n]) - 1;
        }

        if (PRINTINPUT) {
            cout << setw(8) << vbuf[0] << setw(8) << vbuf[1];
            for (USI i = 0; i < 6; i++) {
                cout << setw(6) << index[i] + 1;
            }
            cout << endl;
        }

        {
            auto srcPtr = FindPtr(srcName, (OCP_DBL)0);
            auto objPtr = FindPtr(objName, (OCP_DBL)0);
            if (srcPtr != nullptr && objPtr != nullptr) {
                objPtr->resize(srcPtr->size());
                CopyVal(*objPtr, *srcPtr, index);
                continue;
            }
        }

        {
            auto srcPtr = FindPtr(srcName, (USI)0);
            auto objPtr = FindPtr(objName, (USI)0);
            if (srcPtr != nullptr && objPtr != nullptr) {
                objPtr->resize(srcPtr->size());
                CopyVal(*objPtr, *srcPtr, index);
                continue;
            }
        }

        OCP_ABORT("WRONG Item " + srcName + "  " + objName + " in EQUALS!");
    }
}


void PreParamGridWell::InputMULTIPLY(ifstream& ifs)
{
    if (PRINTINPUT) {
        cout << "\n---------------------" << endl
            << "MULTIPLY"
            << "\n---------------------" << endl;
    }


    vector<string> vbuf;
    vector<USI>    index(6, 0);

    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        index[0] = 0, index[1] = nx - 1;
        index[2] = 0, index[3] = ny - 1;
        index[4] = 0, index[5] = nz - 1;

        string  objName = vbuf[0];
        OCP_DBL val = stod(vbuf[1]);

        DealDefault(vbuf);
        for (USI n = 2; n < 8; n++) {
            if (vbuf[n] != "DEFAULT") index[n - 2] = stoi(vbuf[n]) - 1;
        }

        auto objPtr = FindPtr(objName, (OCP_DBL)0);
        if (objPtr != nullptr) {
            if (objName == "TOPS") {
                index[4] = index[5] = 0;
            }
            MultiplyVal(*objPtr, val, index);
        }
        else {
            OCP_ABORT("Wrong object name: " + objName);
        }

        if (PRINTINPUT) {
            cout << setw(8) << vbuf[0] << setw(8) << vbuf[1];
            for (USI i = 0; i < 6; i++) {
                cout << setw(6) << index[i] + 1;
            }
            cout << endl;
        }
    }
}


void PreParamGridWell::InputGrid(ifstream& ifs, string& keyword)
{
    vector<string> vbuf;

    {
        auto objPtr = FindPtr(keyword, (OCP_DBL)0);
        if (objPtr != nullptr) {
            while (ReadLine(ifs, vbuf)) {
                if (vbuf[0] == "/") break;

                for (auto& str : vbuf) {
                    // if m*n occurs, then push back n  m times
                    auto pos = str.find('*');
                    if (pos == string::npos) {
                        objPtr->push_back(stod(str));
                    }
                    else {
                        const USI     len = str.size();
                        const OCP_ULL num = stoi(str.substr(0, pos));
                        const OCP_DBL val = stod(str.substr(pos + 1, len - (pos + 1)));
                        for (OCP_ULL i = 0; i < num; i++) objPtr->push_back(val);
                    }
                }
            }
            return;
        }        
    }

    {
        auto objPtr = FindPtr(keyword, (USI)0);
        if (objPtr != nullptr) {
            while (ReadLine(ifs, vbuf)) {
                if (vbuf[0] == "/") break;

                for (auto& str : vbuf) {
                    // if m*n occurs, then push back n  m times
                    auto pos = str.find('*');
                    if (pos == string::npos) {
                        objPtr->push_back(stod(str));
                    }
                    else {
                        USI     len = str.size();
                        OCP_ULL num = stoi(str.substr(0, pos));
                        USI val = stoi(str.substr(pos + 1, len - (pos + 1)));
                        for (OCP_ULL i = 0; i < num; i++) objPtr->push_back(val);
                    }
                }
            }
            return;
        }
    }



    OCP_ABORT("Unknown keyword!");
}


void PreParamGridWell::InputINCLUDE(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    DealDefault(vbuf);
    Input(vbuf[0]);
}


#ifdef WITH_GMSH
void PreParamGridWell::InputGMSH(ifstream& ifs)
{
    gridType = GridType::gmsh;
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    DealDefault(vbuf);
    gmshGrid.InputGrid(workdir + vbuf[0]);
}

void PreParamGridWell::InputGMSHPRO(ifstream& ifs)
{
    gmshGrid.InputProperty(ifs);

    // input params
    numGridM = gmshGrid.elements.size();
    numGrid  = numGridM;

    poro.resize(numGrid);
    kx.resize(numGrid);
    ky.resize(numGrid);
    kz.resize(numGrid);
    for (OCP_ULL n = 0; n < numGrid; n++) {
        poro[n] = gmshGrid.facies[gmshGrid.elements[n].phyIndex].poro;
        kx[n]   = gmshGrid.facies[gmshGrid.elements[n].phyIndex].kx;
        ky[n]   = gmshGrid.facies[gmshGrid.elements[n].phyIndex].ky;
        kz[n]   = gmshGrid.facies[gmshGrid.elements[n].phyIndex].kz;
    }
}
#endif

void PreParamGridWell::InputWELSPECS(ifstream& ifs)
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


void PreParamGridWell::InputCOMPDAT(ifstream& ifs)
{
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
        OCP_BOOL tmp = OCP_FALSE;

        for (USI w = 0; w < num; w++) {
            if (fuzzyMatch)
                tmp = (well[w].name.substr(0, pos) == src);
            else
                tmp = (well[w].name == src);

            if (tmp) {
                well[w].InputCOMPDAT(vbuf);
            }
        }
    }
}


vector<OCP_DBL>* PreParamGridWell::FindPtr(const string& varName, const OCP_DBL&)
{
    vector<OCP_DBL>* myPtr = nullptr;

    switch (Map_Str2Int(&varName[0], varName.size())) 
    {
    case Map_Str2Int("DX", 2):
        dx.reserve(numGrid);
        myPtr = &dx;
        break;

    case Map_Str2Int("DY", 2):
        dy.reserve(numGrid);
        myPtr = &dy;
        break;

    case Map_Str2Int("DZ", 2):
        dz.reserve(numGrid);
        myPtr = &dz;
        break;

    case Map_Str2Int("TOPS", 4):
        tops.reserve(nx * ny);
        myPtr = &tops;
        break;

    case Map_Str2Int("COORD", 5):
        gridType = GridType::corner;
        coord.reserve((nx + 1) * (ny + 1) * 6);
        myPtr = &coord;
        break;

    case Map_Str2Int("ZCORN", 5):
        gridType = GridType::corner;
        zcorn.reserve(numGridM * 8);
        myPtr = &zcorn;
        break;

    case Map_Str2Int("PORO", 4):
        poro.reserve(numGrid);
        myPtr = &poro;
        break;

    case Map_Str2Int("NTG", 3):
        ntg.reserve(numGrid);
        myPtr = &ntg;
        break;

    case Map_Str2Int("PERMX", 5):
        kx.reserve(numGrid);
        myPtr = &kx;
        break;

    case Map_Str2Int("PERMY", 5):
        ky.reserve(numGrid);
        myPtr = &ky;
        break;

    case Map_Str2Int("PERMZ", 5):
        kz.reserve(numGrid);
        myPtr = &kz;
        break;

    case Map_Str2Int("SWAT", 4):
        initR.swat.reserve(numGrid);
        myPtr = &initR.swat;
        break;

    case Map_Str2Int("SWATINIT", 8):
        initR.swatInit.reserve(numGrid);
        myPtr = &initR.swatInit;
        initR.scalePcow = OCP_TRUE;
        break;

    case Map_Str2Int("SIGMAV", 6):
        sigma.reserve(numGrid);
        myPtr = &sigma;
        break;

    case Map_Str2Int("MULTZ", 5):
        multZ.reserve(numGrid);
        myPtr = &multZ;
        break;

    case Map_Str2Int("DZMTRXV", 7):
        dzMtrx.reserve(numGrid);
        myPtr = &dzMtrx;
        break;
    }

    return myPtr;
}


vector<USI>* PreParamGridWell::FindPtr(const string& varName, const USI&)
{
    vector<USI>* myPtr = nullptr;

    switch (Map_Str2Int(&varName[0], varName.size())) 
    {
    case Map_Str2Int("ACTNUM", 6):
        actGC.ACTNUM.reserve(numGrid);
        myPtr = &actGC.ACTNUM;
        break;

    case Map_Str2Int("SATNUM", 6):
        SATNUM.reserve(numGrid);
        myPtr = &SATNUM;
        break;

    case Map_Str2Int("PVTNUM", 6):
        PVTNUM.reserve(numGrid);
        myPtr = &PVTNUM;
        break;

    case Map_Str2Int("ROCKNUM", 7):
        ROCKNUM.reserve(numGrid);
        myPtr = &ROCKNUM;
        break;
    }
    return myPtr;
}


template <typename T>
void PreParamGridWell::setVal(vector<T>& obj, const T& val, const vector<USI>& index)
{
    const OCP_ULL nxny = nx * ny;
    OCP_ULL id = 0;

    for (USI k = index[4]; k <= index[5]; k++) {
        for (USI j = index[2]; j <= index[3]; j++) {
            for (USI i = index[0]; i <= index[1]; i++) {
                id = k * nxny + j * nx + i;
                obj[id] = val;
            }
        }
    }
}


template <typename T>
void PreParamGridWell::CopyVal(vector<T>& obj,
    const vector<T>& src,
    const vector<USI>& index)
{

    const OCP_ULL nxny = nx * ny;
    OCP_ULL id = 0;

    for (USI k = index[4]; k <= index[5]; k++) {
        for (USI j = index[2]; j <= index[3]; j++) {
            for (USI i = index[0]; i <= index[1]; i++) {
                id = k * nxny + j * nx + i;
                obj[id] = src[id];
            }
        }
    }
}


void PreParamGridWell::MultiplyVal(vector<OCP_DBL>& obj,
    const OCP_DBL& val,
    const vector<USI>& index)
{
    const OCP_ULL nxny = nx * ny;
    OCP_ULL id = 0;

    for (USI k = index[4]; k <= index[5]; k++) {
        for (USI j = index[2]; j <= index[3]; j++) {
            for (USI i = index[0]; i <= index[1]; i++) {
                id = k * nxny + j * nx + i;
                obj[id] *= val;
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////
// check grid activity
/////////////////////////////////////////////////////////////////////

OCP_ULL ActiveGridCheck::CheckActivity(const OCPModel& Model, const OCP_DBL& ev, const OCP_DBL& ep,
    const vector<OCP_DBL>& v, const vector<OCP_DBL>& poro)
{
    model         = Model;

    eV            = ev;
    eP            = ep;
    activeGridNum = 0;
    numGrid       = v.size();
    map_Act2All.reserve(numGrid);
    map_All2Act.resize(numGrid, -1);

    if (model == OCPModel::isothermal)    CheckActivityIsoT(v, poro);
    else if (model == OCPModel::thermal)  CheckActivityT(v, poro);
    else                                  OCP_ABORT("INAVAILABLE MODEL!");

    cout << "  Number of inactive cells is " << (numGrid - activeGridNum) << " ("
        << (numGrid - activeGridNum) * 100.0 / numGrid << "%)" << endl;

    vector<USI>().swap(ACTNUM);

    return activeGridNum;
}


void ActiveGridCheck::CheckActivityIsoT(const vector<OCP_DBL>& v, const vector<OCP_DBL>& poro)
{
    OCP_ULL activeCount = 0;
	for (OCP_ULL n = 0; n < numGrid; n++) {
		if (poro[n] < eP || v[n] < eV)  continue;
		if (!allAct) {
			if (ACTNUM[n] == 0) continue;
		}

		map_Act2All.push_back(n);
		map_All2Act[n] = activeCount;
		activeCount++;
	}

    activeGridNum = activeCount;
}

void ActiveGridCheck::CheckActivityT(const vector<OCP_DBL>& v, const vector<OCP_DBL>& poro)
{
    OCP_ULL activeCount = 0;
    for (OCP_ULL n = 0; n < numGrid; n++) {
        if (v[n] < eV)  continue;
        if (!allAct) {
            if (ACTNUM[n] == 0) continue;
        }

        map_Act2All.push_back(n);
        map_All2Act[n] = activeCount;
        activeCount++;
    }

    activeGridNum = activeCount;
}


OCP_BOOL ActiveGridCheck::IfFluid(const OCP_ULL& n, const OCP_DBL& poro)
{
    if (map_All2Act[n] >= 0 && poro > eP)  return OCP_TRUE;
    else                                   return OCP_FALSE;
}


void ActiveGridCheck::FreeSomeMemory() 
{
    vector<OCP_SLL>().swap(map_All2Act);
}


/////////////////////////////////////////////////////////////////////
// Initial reservoir data
/////////////////////////////////////////////////////////////////////

void InitialReservoir::CheckData(const OCP_USI& numGrid)
{
    if (swat.size() == 1) {
        swat.resize(numGrid, swat[0]);
    }
    else if (swat.size() != 0 && swat.size() != numGrid) {
        OCP_ABORT("SWAT is not given correctly!");
    }
}


/////////////////////////////////////////////////////////////////////
// Generate active grids' connections
/////////////////////////////////////////////////////////////////////

void PreParamGridWell::Setup()
{
    OCP_INFO("Setup Grid and Well -- begin");

    SetupGrid();
    SetupConnWellGrid();
    initR.CheckData(numGrid);
    actGC.FreeSomeMemory();

    OCP_INFO("Setup Grid and Well -- end");
}


void PreParamGridWell::SetupGrid()
{
    OCP_INFO("Setup Grid -- begin");

    switch (gridType) 
    {
    case GridType::orthogonal:
        SetupOrthogonalGrid();    
        OutputBaiscInfo();
        SetLocationStructral();
        break;
    case GridType::corner:
        SetupCornerGrid();
        OutputBaiscInfo();
        SetLocationStructral();
        break;
#ifdef WITH_GMSH
    case GridType::gmsh:
        SetupGmshGrid();
        break;
#endif
    default:
        OCP_ABORT("WRONG Grid Type!");
    }

    SetupTransMult();

    OCP_INFO("Setup Grid -- end");
}


void PreParamGridWell::SetupOrthogonalGrid()
{
    OCP_INFO("Setup Orthogonal Grid -- begin");

    // x -> y -> z
    CalDepthVOrthogonalGrid();
    CalActiveGrid(1E-6, 1E-6);
    SetupActiveConnOrthogonalGrid();

    OutputPointsOrthogonalGrid();

    OCP_INFO("Setup Orthogonal Grid -- end");
}

void PreParamGridWell::CalDepthVOrthogonalGrid()
{
    OCP_INFO("Calculate Depth and Volume -- begin");

    depth.resize(numGrid, 0);
    const OCP_ULL nxny = nx * ny;
    OCP_ULL id;
    // 0th layer
    for (USI j = 0; j < ny; j++) {
        for (USI i = 0; i < nx; i++) {
            id = j * nx + i;
            depth[id] = tops[id] + dz[id] / 2;
        }
    }
    // 1th - (nz-1)th layer
    for (USI k = 1; k < nz; k++) {
        OCP_ULL knxny = k * nxny;
        for (USI j = 0; j < ny; j++) {
            for (USI i = 0; i < nx; i++) {
                id = knxny + j * nx + i;
                depth[id] = depth[id - nxny] + dz[id - nxny] / 2 + dz[id] / 2;
            }
        }
    }

    if (DUALPORO) {
        copy(&depth[0], &depth[numGridM], &depth[numGridM]);
    }

    v.resize(numGrid);
    for (OCP_ULL i = 0; i < numGrid; i++) v[i] = dx[i] * dy[i] * dz[i];

    vector<OCP_DBL>().swap(tops);

    OCP_INFO("Calculate Depth and Volume -- end");
}

void PreParamGridWell::SetupActiveConnOrthogonalGrid()
{
    OCP_INFO("Setup Active Grid Connections -- begin");

    if (DUALPORO) {
        SetupActiveConnOrthogonalGridDP();
    }
    else {
        SetupActiveConnOrthogonalGridSM();
    }

    OCP_INFO("Setup Active Grid Connections -- end");
}


void PreParamGridWell::SetupActiveConnOrthogonalGridSM()
{
    gNeighbor.resize(activeGridNum);
    // PreAllocate
    for (OCP_ULL n = 0; n < activeGridNum; n++) {
        gNeighbor[n].reserve(6);
    }

    // Begin Id and End Id in Grid, bIdg < eIdg
    OCP_SLL       bIdg, eIdg, bIdb, eIdb;
    OCP_DBL       areaB, areaE;
    const OCP_ULL nxny = nx * ny;

    for (USI k = 0; k < nz; k++) {
        for (USI j = 0; j < ny; j++) {
            for (USI i = 0; i < nx; i++) {

                bIdg = k * nxny + j * nx + i;

                bIdb = actGC.map_All2Act[bIdg];
                if (bIdb < 0)  continue;

                // right  --  x-direction
                if (i < nx - 1) {
                    eIdg = bIdg + 1;
                    eIdb = actGC.map_All2Act[eIdg];
                    if (eIdb < 0)  continue;

					areaB = 2 * dy[bIdg] * dz[bIdg] / dx[bIdg];
					areaE = 2 * dy[eIdg] * dz[eIdg] / dx[eIdg];
					gNeighbor[bIdb].push_back(ConnPair(eIdb, WEIGHT_GG, ConnDirect::xp, areaB, areaE));
					gNeighbor[eIdb].push_back(ConnPair(bIdb, WEIGHT_GG, ConnDirect::xm, areaE, areaB));
                }
                // front  --  y-direction
                if (j < ny - 1) {
                    eIdg = bIdg + nx;
                    eIdb = actGC.map_All2Act[eIdg];
                    if (eIdb < 0)  continue;

					areaB = 2 * dz[bIdg] * dx[bIdg] / dy[bIdg];
					areaE = 2 * dz[eIdg] * dx[eIdg] / dy[eIdg];
					gNeighbor[bIdb].push_back(ConnPair(eIdb, WEIGHT_GG, ConnDirect::yp, areaB, areaE));
					gNeighbor[eIdb].push_back(ConnPair(bIdb, WEIGHT_GG, ConnDirect::ym, areaE, areaB));
                }
                // down --   z-direction
                if (k < nz - 1) {
                    eIdg = bIdg + nxny;
                    eIdb = actGC.map_All2Act[eIdg];
                    if (eIdb < 0)  continue;

					areaB = 2 * dx[bIdg] * dy[bIdg] / dz[bIdg];
					areaE = 2 * dx[eIdg] * dy[eIdg] / dz[eIdg];
					gNeighbor[bIdb].push_back(ConnPair(eIdb, WEIGHT_GG, ConnDirect::zp, areaB, areaE));
					gNeighbor[eIdb].push_back(ConnPair(bIdb, WEIGHT_GG, ConnDirect::zm, areaE, areaB));
                }
            }
        }
    }
}


void PreParamGridWell::SetupActiveConnOrthogonalGridDP()
{

    // for fractures connection

    gNeighbor.resize(activeGridNum);
    // PreAllocate
    for (OCP_ULL n = 0; n < activeGridNum; n++) {
        gNeighbor[n].reserve(6);
    }

    // Begin Id and End Id in Grid, bIdg < eIdg
    OCP_SLL       bIdg, eIdg, bIdb, eIdb;
    OCP_DBL       areaB, areaE;
    const OCP_ULL nxny = nx * ny;

    for (USI k = nz; k < 2 * nz; k++) {
        for (USI j = 0; j < ny; j++) {
            for (USI i = 0; i < nx; i++) {

                bIdg = k * nxny + j * nx + i;
                bIdb = actGC.map_All2Act[bIdg];
                if (bIdb < 0)  continue;

                // right  --  x-direction
                if (i < nx - 1) {
                    eIdg = bIdg + 1;
                    eIdb = actGC.map_All2Act[eIdg];
                    if (eIdb < 0)  continue;

					areaB = 2 * dy[bIdg] * dz[bIdg] / dx[bIdg];
					areaE = 2 * dy[eIdg] * dz[eIdg] / dx[eIdg];
					gNeighbor[bIdb].push_back(ConnPair(eIdb, WEIGHT_GG, ConnDirect::xp, areaB, areaE));
					gNeighbor[eIdb].push_back(ConnPair(bIdb, WEIGHT_GG, ConnDirect::xm, areaE, areaB));
                }
                // front  --  y-direction
                if (j < ny - 1) {
                    eIdg = bIdg + nx;
                    eIdb = actGC.map_All2Act[eIdg];
                    if (eIdb < 0)  continue;

					areaB = 2 * dz[bIdg] * dx[bIdg] / dy[bIdg];
					areaE = 2 * dz[eIdg] * dx[eIdg] / dy[eIdg];
					gNeighbor[bIdb].push_back(ConnPair(eIdb, WEIGHT_GG, ConnDirect::yp, areaB, areaE));
					gNeighbor[eIdb].push_back(ConnPair(bIdb, WEIGHT_GG, ConnDirect::ym, areaE, areaB));
                }
                // down --   z-direction
                if (k < 2*nz - 1) {
                    eIdg = bIdg + nxny;
                    eIdb = actGC.map_All2Act[eIdg];
                    if (eIdb < 0)  continue;

					areaB = 2 * dx[bIdg] * dy[bIdg] / dz[bIdg];
					areaE = 2 * dx[eIdg] * dy[eIdg] / dz[eIdg];
					gNeighbor[bIdb].push_back(ConnPair(eIdb, WEIGHT_GG, ConnDirect::zp, areaB, areaE));
					gNeighbor[eIdb].push_back(ConnPair(bIdb, WEIGHT_GG, ConnDirect::zm, areaE, areaB));
                }
            }
        }
    }

    // for fracture-matrix connection
    for (bIdg = 0; bIdg < numGridM; bIdg++) {
        bIdb = actGC.map_All2Act[bIdg];
        if (bIdb < 0)  continue;
    
        eIdg = bIdg + numGridM;
        eIdb = actGC.map_All2Act[eIdg];
        if (eIdb < 0)  continue;
    
		gNeighbor[bIdb].push_back(ConnPair(eIdb, WEIGHT_GG, ConnDirect::mf, 0.0, 0.0));
		gNeighbor[eIdb].push_back(ConnPair(bIdb, WEIGHT_GG, ConnDirect::fm, 0.0, 0.0));
    }
}


void PreParamGridWell::OutputPointsOrthogonalGrid()
{
    if (!ifUseVtk)  return;

    OCP_INFO("Ouput Orthogonal Grid points for vtk -- begin");

    vector<OCP_DBL> points_xyz;  ///< x,y,z coordinates
    vector<OCP_ULL> cell_points; ///< numpoints, points index
    vector<USI>     cell_type;   ///< type of cell
    points_xyz.reserve(activeGridNum * 8 * 3);
    cell_points.reserve(activeGridNum * 9);
    cell_type.resize(activeGridNum, VTK_HEXAHEDRON);

    OCP_ULL         pIndex = 0;
    OCP_DBL         tmpX, tmpY;
    OCP_ULL         id;
    for (USI k = 0; k < nz; k++) {
        tmpY = 0;
        for (USI j = 0; j < ny; j++) {
            tmpX = 0;
            for (USI i = 0; i < nx; i++) {
                id = k * nx * ny + j * nx + i;
                if (actGC.map_All2Act[id] >= 0) {
                    points_xyz.push_back(tmpX);
                    points_xyz.push_back(tmpY);
                    points_xyz.push_back(depth[id] + dz[id] / 2);

                    points_xyz.push_back(tmpX + dx[id]);
                    points_xyz.push_back(tmpY);
                    points_xyz.push_back(depth[id] + dz[id] / 2);

                    points_xyz.push_back(tmpX + dx[id]);
                    points_xyz.push_back(tmpY + dy[id]);
                    points_xyz.push_back(depth[id] + dz[id] / 2);

                    points_xyz.push_back(tmpX);
                    points_xyz.push_back(tmpY + dy[id]);
                    points_xyz.push_back(depth[id] + dz[id] / 2);

                    points_xyz.push_back(tmpX);
                    points_xyz.push_back(tmpY);
                    points_xyz.push_back(depth[id] - dz[id] / 2);

                    points_xyz.push_back(tmpX + dx[id]);
                    points_xyz.push_back(tmpY);
                    points_xyz.push_back(depth[id] - dz[id] / 2);

                    points_xyz.push_back(tmpX + dx[id]);
                    points_xyz.push_back(tmpY + dy[id]);
                    points_xyz.push_back(depth[id] - dz[id] / 2);

                    points_xyz.push_back(tmpX);
                    points_xyz.push_back(tmpY + dy[id]);
                    points_xyz.push_back(depth[id] - dz[id] / 2);


                    cell_points.push_back(8);
                    for (USI p = 0; p < 8; p++) {
                        cell_points.push_back(pIndex++);
                    }
                }
                tmpX += dx[id];
            }
            tmpY += dy[id];
        }
    }

    OCP_ASSERT(points_xyz.size() == activeGridNum * 8 * 3, "WRONG OutputPointsOrthogonalGrid!");

    Output4Vtk::OutputGridInfo(workdir, activeGridNum, points_xyz, cell_points, cell_type);

    OCP_INFO("Ouput Orthogonal Grid points for vtk -- end");
}


void PreParamGridWell::SetupCornerGrid()
{
    OCP_COORD coordTmp;
    coordTmp.Allocate(nx, ny, nz);
    coordTmp.InputData(coord, zcorn);
    coordTmp.SetupCornerPoints();
    SetupBasicCornerGrid(coordTmp);
    CalActiveGrid(1E-6, 1E-6);
    SetupActiveConnCornerGrid(coordTmp);

    OutputPointsCornerGrid(coordTmp);

    vector<OCP_DBL>().swap(coord);
    vector<OCP_DBL>().swap(zcorn);
}

void PreParamGridWell::SetupBasicCornerGrid(const OCP_COORD& CoTmp)
{
    if (DUALPORO) {
                                  
        dx.resize(numGrid);
        copy(CoTmp.dx.begin(), CoTmp.dx.end(), &dx[0]);
        copy(CoTmp.dx.begin(), CoTmp.dx.end(), &dx[numGridM]);

        dy.resize(numGrid);
        copy(CoTmp.dy.begin(), CoTmp.dy.end(), &dy[0]);
        copy(CoTmp.dy.begin(), CoTmp.dy.end(), &dy[numGridM]);

        dz.resize(numGrid);
        copy(CoTmp.dz.begin(), CoTmp.dz.end(), &dz[0]);
        copy(CoTmp.dz.begin(), CoTmp.dz.end(), &dz[numGridM]);

        v.resize(numGrid);
        copy(CoTmp.v.begin(), CoTmp.v.end(), &v[0]);
        copy(CoTmp.v.begin(), CoTmp.v.end(), &v[numGridM]);

        depth.resize(numGrid);
        copy(CoTmp.depth.begin(), CoTmp.depth.end(), &depth[0]);
        copy(CoTmp.depth.begin(), CoTmp.depth.end(), &depth[numGridM]);
    }
    else {
        dx    = CoTmp.dx;
        dy    = CoTmp.dy;
        dz    = CoTmp.dz;
        v     = CoTmp.v;
        depth = CoTmp.depth;
    }
}

void PreParamGridWell::SetupActiveConnCornerGrid(const OCP_COORD& CoTmp)
{
    if (DUALPORO) {
        SetupActiveConnCornerGridDP(CoTmp);
    }
    else {
        SetupActiveConnCornerGridSM(CoTmp);
    }
}


void PreParamGridWell::SetupActiveConnCornerGridSM(const OCP_COORD& CoTmp)
{
    gNeighbor.resize(activeGridNum);
    // PreAllocate
    for (OCP_ULL n = 0; n < activeGridNum; n++) {
        gNeighbor[n].reserve(10);
    }

    OCP_SLL bIdg, eIdg, bIdb, eIdb;
    OCP_DBL areaB, areaE;
    for (OCP_ULL n = 0; n < CoTmp.numConn; n++) {
        const GeneralConnect& ConnTmp = CoTmp.connect[n];

        bIdg = ConnTmp.begin;
        eIdg = ConnTmp.end;

        bIdb = actGC.map_All2Act[bIdg];
        eIdb = actGC.map_All2Act[eIdg];

        if (bIdb >= 0 && eIdb >= 0) {          
            areaB = ConnTmp.Ad_dd_begin;
            areaE = ConnTmp.Ad_dd_end;
            gNeighbor[bIdb].push_back(ConnPair(eIdb, WEIGHT_GG, ConnTmp.directionType, areaB, areaE));
        }
    }
}


void PreParamGridWell::SetupActiveConnCornerGridDP(const OCP_COORD& CoTmp)
{
    // for fracture connections
    gNeighbor.resize(activeGridNum);
    // PreAllocate
    for (OCP_ULL n = 0; n < activeGridNum; n++) {
        gNeighbor[n].reserve(10);
    }

    OCP_SLL bIdg, eIdg, bIdb, eIdb;
    OCP_DBL areaB, areaE;
    for (OCP_ULL n = 0; n < CoTmp.numConn; n++) {
        const GeneralConnect& ConnTmp = CoTmp.connect[n];

        bIdg = ConnTmp.begin + numGridM;
        eIdg = ConnTmp.end + numGridM;

        bIdb = actGC.map_All2Act[bIdg];
        eIdb = actGC.map_All2Act[eIdg];

        if (bIdb >= 0 && eIdb >= 0) {
            areaB = ConnTmp.Ad_dd_begin;
            areaE = ConnTmp.Ad_dd_end;
            gNeighbor[bIdb].push_back(ConnPair(eIdb, WEIGHT_GG, ConnTmp.directionType, areaB, areaE));
        }
    }

    // for fracture-matrix connection
    for (bIdg = 0; bIdg < numGridM; bIdg++) {
        bIdb = actGC.map_All2Act[bIdg];
        if (bIdb < 0)  continue;

        eIdg = bIdg + numGridM;
        eIdb = actGC.map_All2Act[eIdg];
        if (eIdb >= 0) {

            gNeighbor[bIdb].push_back(ConnPair(eIdb, WEIGHT_GG, ConnDirect::mf, 0.0, 0.0));
            gNeighbor[eIdb].push_back(ConnPair(bIdb, WEIGHT_GG, ConnDirect::fm, 0.0, 0.0));
        }
    }
}


void PreParamGridWell::OutputPointsCornerGrid(const OCP_COORD& mycord)
{
    if (!ifUseVtk)  return;

    vector<OCP_DBL> points_xyz;  ///< x,y,z coordinates
    vector<OCP_ULL> cell_points; ///< numpoints, points index
    vector<USI>     cell_type;   ///< type of cell
    points_xyz.reserve(activeGridNum * 8 * 3);
    cell_points.reserve(activeGridNum * 9);
    cell_type.resize(activeGridNum, VTK_HEXAHEDRON);

    OCP_ULL pIndex = 0;
    OCP_ULL id;
    for (USI k = 0; k < nz; k++) {
        for (USI j = 0; j < ny; j++) {
            for (USI i = 0; i < nx; i++) {
                id = k * nx * ny + j * nx + i;
                if (actGC.map_All2Act[id] >= 0) {
                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p4.x);
                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p4.y);
                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p4.z);

                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p5.x);
                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p5.y);
                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p5.z);

                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p6.x);
                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p6.y);
                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p6.z);

                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p7.x);
                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p7.y);
                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p7.z);

                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p0.x);
                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p0.y);
                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p0.z);

                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p1.x);
                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p1.y);
                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p1.z);

                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p2.x);
                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p2.y);
                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p2.z);

                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p3.x);
                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p3.y);
                    points_xyz.push_back(mycord.cornerPoints[i][j][k].p3.z);

                    cell_points.push_back(8);
                    for (USI p = 0; p < 8; p++) {
                        cell_points.push_back(pIndex++);
                    }
                }
            }
        }
    }
    OCP_ASSERT(points_xyz.size() == activeGridNum * 8 * 3, "WRONG OutputPointsOrthogonalGrid!");

    Output4Vtk::OutputGridInfo(workdir, activeGridNum, points_xyz, cell_points, cell_type);
}


#ifdef WITH_GMSH
void PreParamGridWell::SetupGmshGrid()
{
    SetupBasicGmshGrid();
    CalActiveGrid(1E-12, 1E-12);
    SetupActiveConnGmshGrid();
    OutputPointsGmshGrid();
}


void PreParamGridWell::SetupBasicGmshGrid()
{
    numGridM = gmshGrid.elements.size();
    numGrid  = numGridM;
    v.resize(numGrid);
    depth.resize(numGrid);
    boundIndex.resize(numGrid);
    boundArea.resize(numGrid);
    SATNUM.resize(numGrid);
    PVTNUM.resize(numGrid);
    ROCKNUM.resize(numGrid);

    if (gmshGrid.dimen == 2) {
        for (OCP_ULL n = 0; n < numGridM; n++) {
            v[n]          = gmshGrid.elements[n].area * gmshGrid.thickness;
            depth[n]      = gmshGrid.elements[n].center.y; /// Use y-coordinate
            SATNUM[n]     = gmshGrid.elements[n].phyIndex;
            PVTNUM[n]     = gmshGrid.elements[n].phyIndex;
            ROCKNUM[n]    = gmshGrid.elements[n].phyIndex;
            boundIndex[n] = gmshGrid.elements[n].boundIndex;      
            boundArea[n]  = gmshGrid.elements[n].boundArea;
        }
    }
}


void PreParamGridWell::SetupActiveConnGmshGrid()
{

    OCP_DBL thickNess = 1.0;
    if (gmshGrid.dimen == 2) {
        thickNess = gmshGrid.thickness;
    }

    gNeighbor.resize(activeGridNum);
    // PreAllocate
    for (OCP_ULL n = 0; n < activeGridNum; n++) {
        gNeighbor[n].reserve(10);
    }

    OCP_SLL bIdg, eIdg, bIdb, eIdb;
    OCP_DBL areaB, areaE;
    for (const auto& e : gmshGrid.edges) {
       
        if (e.faceIndex.size() <= 2)  continue;  // boundary

        bIdg = e.faceIndex[0];
        eIdg = e.faceIndex[2];

        bIdb = actGC.map_All2Act[bIdg];
        eIdb = actGC.map_All2Act[eIdg];

        if (bIdb >= 0 && eIdb >= 0) {
            
            areaB = e.area[0] * thickNess;
            areaE = e.area[1] * thickNess;
            gNeighbor[bIdb].push_back(ConnPair(eIdb, WEIGHT_GG, ConnDirect::usg, areaB, areaE));
            gNeighbor[eIdb].push_back(ConnPair(bIdb, WEIGHT_GG, ConnDirect::usg, areaE, areaB));
        }
    }
}


/// Output grid points for a gmsh grid
void PreParamGridWell::OutputPointsGmshGrid()
{
    if (!ifUseVtk)  return;

    const vector<OCP_DBL>& points_xyz = gmshGrid.points;  ///< x,y,z coordinates
    vector<OCP_ULL>        cell_points;                   ///< numpoints, points index
    vector<USI>            cell_type;                     ///< type of cell
    cell_points.reserve(activeGridNum * 5);
    cell_type.reserve(activeGridNum);

    for (OCP_ULL n = 0; n < numGrid; n++) {
        if (actGC.map_All2Act[n] >= 0) {

            const auto& ep = gmshGrid.elements[n].p;

            if (ep.size() == 3) {
                cell_type.push_back(VTK_TRIANGLE);
                cell_points.push_back(3);
            }
            else if (ep.size() == 4) {
                cell_type.push_back(VTK_QUAD);
                cell_points.push_back(4);
            }
            for (const auto& p : ep) {
                cell_points.push_back(p);
            }
        }
    }


    Output4Vtk::OutputGridInfo(workdir, activeGridNum, points_xyz, cell_points, cell_type);
}
#endif


void PreParamGridWell::SetLocationStructral()
{
    boundIndex.resize(numGrid);
    const OCP_ULL uplim   = nx * ny;
    const OCP_ULL downlim = nx * ny * (nz - 1);
    for (OCP_ULL n = 0; n < uplim; n++) {
        boundIndex[n] = 1;
    }
    for (OCP_ULL n = downlim; n < nx * ny * nz; n++) {
        boundIndex[n] = 2;
    }
}


void PreParamGridWell::CalActiveGrid(const OCP_DBL& ev, const OCP_DBL& ep)
{
    OCP_INFO("Select Active Grid -- begin");

    activeGridNum = actGC.CheckActivity(model, ev, ep, v, poro);

    OCP_INFO("Select Active Grid -- end");
}


void PreParamGridWell::SetupTransMult()
{
    if (!multZ.empty()) {
        for (OCP_ULL n = 0; n < numGrid; n++) {
            for (auto& c : gNeighbor[n]) {
                if (c.ID() < n) continue;
                if (c.Direct() == ConnDirect::zp) {
                    c.SetTransMult(multZ[n]);
                    for (auto& c1 : gNeighbor[c.ID()]) {
                        if (c1.ID() == n) {
                            c1.SetTransMult(multZ[n]);
                            break;
                        }
                    }
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////
// Generate connections between active grids and wells
/////////////////////////////////////////////////////////////////////


void PreParamGridWell::SetupConnWellGrid()
{
    OCP_INFO("Setup Connection between Grid and Well -- begin");

    // Attention that all wells should be active -- own at least one connections to active grid
    numWell = well.size();
    connWellGrid.resize(numWell);
    for (USI w = 0; w < numWell; w++) {
        const USI numPerf = well[w].GetPerfNum();
        for (USI p = 0; p < numPerf; p++) {
            const OCP_ULL pId = GetPerfLocation(well[w], p);
            if (actGC.IfFluid(pId, poro[pId])) {
                connWellGrid[w].push_back(actGC.map_All2Act[pId]);
                // for well-connection, areaB and areaE contains its active perforation index and trans if necessary
                gNeighbor[actGC.map_All2Act[pId]].push_back(ConnPair(w + activeGridNum, WEIGHT_GW, ConnDirect::n, p, 0));
            }
        }
        if (connWellGrid[w].empty()) {
            OCP_ABORT("All perforations of Well " + well[w].name + " are in Inactive grid!");
        }
    }


    const OCP_ULL numTotal = activeGridNum + numWell;
    // Add well-grid connections to gNeighbor
    gNeighbor.resize(numTotal);
    for (USI w = 0; w < numWell; w++) {
        for (const auto& b : connWellGrid[w]) {
            gNeighbor[w + activeGridNum].push_back(ConnPair(b, WEIGHT_GW, ConnDirect::n, 0, 0));
        }
    }
    gNeighbor.shrink_to_fit();
    
    OCP_INFO("Setup Connection between Grid and Well -- end");
}


OCP_ULL PreParamGridWell::GetPerfLocation(const WellParam& well, const USI& p)
{
    if (gridType >= GridType::structured && gridType < GridType::unstructured) {
        return (well.K_perf[p] - 1) * (nx * ny) + (well.J_perf[p] - 1) * nx + (well.I_perf[p] - 1);
    }
#ifdef WITH_GMSH
    else if (gridType >= GridType::unstructured) {
        // find the element whose center is closest to the perforation first
        OCP_DBL mindis = 1E8;
        OCP_ULL bId    = 0;
        const Point3D&& pl = Point3D(well.X_perf[p], well.Y_perf[p], well.Z_perf[p]);
        for (OCP_ULL n = 0; n < gmshGrid.elements.size(); n++) {
            const auto& e = gmshGrid.elements[n];
            const OCP_DBL dis = (e.center - pl) * (e.center - pl);
            if (dis < mindis) {
                mindis = dis;
                bId    = n;
            }
        }

        // check if perforation is in the found element
        const OCP_BOOL flag = gmshGrid.elements[bId].IfPointInElement(pl, gmshGrid.points);
        if (!flag) OCP_ABORT("NEED MORE CHECK!");
        return bId;
    }
#endif
    else {
        OCP_ABORT("INAVAILABLE WELL TYPE!");
    }
}


/////////////////////////////////////////////////////////////////////
// Output basic grid information and grid connections
/////////////////////////////////////////////////////////////////////



void PreParamGridWell::OutputBaiscInfo() const
{
    OCP_DBL depthMax = 0;
    OCP_DBL depthMin = 1E8;
    OCP_DBL dxMax = 0;
    OCP_DBL dxMin = 1E8;
    OCP_DBL dyMax = 0;
    OCP_DBL dyMin = 1E8;
    OCP_DBL dzMax = 0;
    OCP_DBL dzMin = 1E8;

    for (OCP_ULL n = 0; n < numGrid; n++) {
        if (depthMax < depth[n]) {
            depthMax = depth[n];
        }
        if (depthMin > depth[n]) {
            depthMin = depth[n];
        }
        if (dxMax < dx[n]) {
            dxMax = dx[n];
        }
        if (dxMin > dx[n]) {
            dxMin = dx[n];
        }
        if (dyMax < dy[n]) {
            dyMax = dy[n];
        }
        if (dyMin > dy[n]) {
            dyMin = dy[n];
        }
        if (dzMax < dz[n]) {
            dzMax = dz[n];
        }
        if (dzMin > dz[n]) {
            dzMin = dz[n];
        }
    }

    cout << "\n---------------------" << endl
        << "GRID"
        << "\n---------------------" << endl;
    cout << "  depthMax = " << depthMax << endl
        << "  depthMin = " << depthMin << endl
        << "  dxMax    = " << dxMax << endl
        << "  dxMin    = " << dxMin << endl
        << "  dyMax    = " << dyMax << endl
        << "  dyMin    = " << dyMin << endl
        << "  dzMax    = " << dzMax << endl
        << "  dzMin    = " << dzMin << endl;
}



void PreParamGridWell::FreeMemory()
{
    vector<OCP_DBL>().swap(dx);
    vector<OCP_DBL>().swap(dy);
    vector<OCP_DBL>().swap(dz);
    vector<OCP_DBL>().swap(ntg);
    vector<OCP_DBL>().swap(poro);
    vector<OCP_DBL>().swap(kx);
    vector<OCP_DBL>().swap(ky);
    vector<OCP_DBL>().swap(kz);
    vector<OCP_DBL>().swap(sigma);
    vector<OCP_DBL>().swap(dzMtrx);
    vector<OCP_DBL>().swap(initR.swat);
    vector<OCP_DBL>().swap(initR.swatInit);
    vector<OCP_DBL>().swap(multZ);
    vector<USI>().swap(SATNUM);
    vector<USI>().swap(PVTNUM);
    vector<USI>().swap(ROCKNUM);
    vector<USI>().swap(boundIndex);
    vector<OCP_DBL>().swap(boundArea);

    vector<WellParam>().swap(well);

    vector<OCP_DBL>().swap(v);
    vector<OCP_DBL>().swap(depth);
    vector<vector<ConnPair>>().swap(gNeighbor);
    vector<OCP_ULL>().swap(actGC.map_Act2All);

    vector<vector<OCP_ULL>>().swap(connWellGrid);
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Feb/15/2023      Create file                          */
/*----------------------------------------------------------------------------*/