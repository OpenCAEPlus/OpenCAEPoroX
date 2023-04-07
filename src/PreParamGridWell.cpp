/*! \file    PreParamGridWell.cpp
 *  \brief   PreParamGridWell class definition
 *  \author  Shizhe Li
 *  \date    Feb/15/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "PreParamGridWell.hpp"


void PreParamGridWell::InputFile(const string& myFile, const string& myWorkdir)
{
    workdir = myWorkdir;
    Input(myFile);
    CheckInput();
}


void PreParamGridWell::Input(const string& myFilename)
{
    cout << workdir + myFilename << endl;

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
                InputBasic(ifs, keyword);
                break;

            case Map_Str2Int("ACTNUM", 6):
                InputRegion(ifs, keyword);
                break;

            case Map_Str2Int("INCLUDE", 7):
                InputINCLUDE(ifs);
                break;

            case Map_Str2Int("WELSPECS", 8):
                InputWELSPECS(ifs);
                break;

            case Map_Str2Int("COMPDAT", 7):
                InputCOMPDAT(ifs);
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
    if (model == 0)                       OCP_ABORT("WRONG MODEL!");
    if (nx == 0 || ny == 0 || nz == 0)    OCP_ABORT("WRONG DIMENS!");
    if (poro.size() != numGrid)           OCP_ABORT("WRONG PORO!");
    if (ntg.size() != numGrid) {
        OCP_WARNING("NTG will be set to 1 !");
    }
    else {
        for (OCP_USI n = 0; n < numGrid; n++) {
            poro[n] *= ntg[n];
        }
    }
    if (ACTNUM.size() != numGrid) {
        OCP_WARNING("ACTNUM will be set to 1 !");
        ACTNUM.clear();
        ACTNUM.resize(numGrid, 1);
    }
    if (gridType == CORNER_GRID) {
        if (zcorn.empty())                OCP_ABORT("WRONG ZCORN!");
        if (coord.empty())                OCP_ABORT("WRONG ZCORN!");
    }
    else if (gridType == ORTHOGONAL_GRID) {
        if (dx.size() != numGrid)         OCP_ABORT("WRONG DX!");
        if (dy.size() != numGrid)         OCP_ABORT("WRONG DY!");
        if (dz.size() != numGrid)         OCP_ABORT("WRONG DZ!");
        if (tops.size() != nx * ny)       OCP_ABORT("WRONG TOPS!");
    }
    else                                  OCP_ABORT("WRONG Grid Type!");

    cout << "Check Grid param ... done";
    cout << endl << "-------------------------------------" << endl;
}


void PreParamGridWell::InputMODEL(ifstream& ifs)
{
    cout << "MODEL" << endl;

    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (vbuf[0] == "THERMAL") {
        model = THERMAL;
    }
    else if (vbuf[0] == "ISOTHERMAL")
    {
        model = ISOTHERMALMODEL;
    }
    else {
        OCP_ABORT("WRONG MODEL in keyword MODEL!");
    }

    cout << vbuf[0] << endl << endl;
}


void PreParamGridWell::InputDIMENS(ifstream& ifs)
{
    cout << "DIMENS" << endl;

    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    nx = stoi(vbuf[0]);
    ny = stoi(vbuf[1]);
    nz = stoi(vbuf[2]);
    numGrid = nx * ny * nz;

    cout << setw(6) << nx << setw(6) << ny << setw(6) << nz << endl << endl;
}


void PreParamGridWell::InputEQUALS(ifstream& ifs)
{
    cout << "\n---------------------" << endl
        << "EQUALS"
        << "\n---------------------" << endl;

    vector<USI>    index(6, 0);
    vector<string> vbuf;

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

        if (index[0] < 0 || index[2] < 0 || index[4] < 0 || index[1] > nx - 1 ||
            index[3] > ny - 1 || index[5] > nz - 1) {
            OCP_ABORT("WRONG Range in " + objName + " in EQUALS!");
        }

        vector<OCP_DBL>* objPtr = FindPtr(objName);

        if (objPtr != nullptr) {
            objPtr->resize(objPtr->capacity());
            if (objName == "TOPS") {
                // objPtr->resize(nx * ny);               
                index[4] = index[5] = 0;
            }
            else {
                // objPtr->resize(numGrid);
            }
            setVal(*objPtr, val, index);
        }
        else {
            continue;
            // OCP_ABORT("Wrong object name: " + objName);
        }

        cout << setw(8) << vbuf[0] << setw(16) << vbuf[1];
        for (USI i = 0; i < 6; i++) {
            cout << setw(6) << index[i] + 1;
        }
        cout << endl;
    }

    cout << "/" << endl;
}


void PreParamGridWell::InputCOPY(ifstream& ifs)
{
    cout << "\n---------------------" << endl
        << "COPY"
        << "\n---------------------" << endl;

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

        auto srcPtr = FindPtr(srcName);
        auto objPtr = FindPtr(objName);
        if (srcPtr != nullptr && objPtr != nullptr) {
            objPtr->resize(srcPtr->size());
            CopyVal(*objPtr, *srcPtr, index);
        }
        else {
            continue;
            // OCP_ABORT("Wrong object names: " + srcName + ", " + objName);
        }

        cout << setw(8) << vbuf[0] << setw(8) << vbuf[1];
        for (USI i = 0; i < 6; i++) {
            cout << setw(6) << index[i] + 1;
        }
        cout << endl;
    }
}


void PreParamGridWell::InputMULTIPLY(ifstream& ifs)
{
    cout << "\n---------------------" << endl
        << "MULTIPLY"
        << "\n---------------------" << endl;


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

        auto objPtr = FindPtr(objName);
        if (objPtr != nullptr) {
            if (objName == "TOPS") {
                index[4] = index[5] = 0;
            }
            MultiplyVal(*objPtr, val, index);
        }
        else {
            continue;
            // OCP_ABORT("Wrong object name: " + objName);
        }

        cout << setw(8) << vbuf[0] << setw(8) << vbuf[1];
        for (USI i = 0; i < 6; i++) {
            cout << setw(6) << index[i] + 1;
        }
        cout << endl;
    }
}


void PreParamGridWell::InputBasic(ifstream& ifs, string& keyword)
{
    auto objPtr = FindPtr(keyword);
    if (objPtr == nullptr) {
        OCP_ABORT("Unknown keyword!");
    }

    vector<string> vbuf;
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
                OCP_USI num = stoi(str.substr(0, pos));
                OCP_DBL val = stod(str.substr(pos + 1, len - (pos + 1)));
                for (USI i = 0; i < num; i++) objPtr->push_back(val);
            }
        }
    }
}


void PreParamGridWell::InputRegion(ifstream& ifs, const string& keyword)
{
    auto ptr = &ACTNUM;

    ptr->reserve(numGrid);
    vector<string>  vbuf;
    vector<OCP_USI> obj;
    vector<USI>     region;

    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        for (auto& str : vbuf) {
            // if m*n occurs, then push back n m times
            auto pos = str.find('*');
            if (pos == string::npos) {
                ptr->push_back(stod(str));
            }
            else {
                USI     len = str.size();
                OCP_USI num = stoi(str.substr(0, pos));
                OCP_DBL val = stod(str.substr(pos + 1, len - (pos + 1)));
                for (USI i = 0; i < num; i++) ptr->push_back(val);
            }
        }
    }
}


void PreParamGridWell::InputINCLUDE(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    DealDefault(vbuf);
    Input(vbuf[0]);
}


void PreParamGridWell::InputWELSPECS(ifstream& ifs)
{
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        DealDefault(vbuf);
        well.push_back(PreParamWell(vbuf));
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

                const USI k1 = stoi(vbuf[3]);
                const USI k2 = stoi(vbuf[4]);

                for (USI k = k1; k <= k2; k++) {
                    if (vbuf[1] == "DEFAULT" || vbuf[2] == "DEFAULT") {
                        well[w].I_perf.push_back(well[w].I);
                        well[w].J_perf.push_back(well[w].J);
                    }
                    else {
                        well[w].I_perf.push_back(stoi(vbuf[1]));
                        well[w].J_perf.push_back(stoi(vbuf[2]));
                    }
                    well[w].K_perf.push_back(k);

                    if (vbuf[5] != "DEFAULT")
                        well[w].WI.push_back(stod(vbuf[5]));
                    else
                        well[w].WI.push_back(-1.0);

                    if (vbuf[6] != "DEFAULT")
                        well[w].diameter.push_back(stod(vbuf[6]));
                    else
                        well[w].diameter.push_back(1.0);

                    if (vbuf[7] != "DEFAULT")
                        well[w].kh.push_back(stod(vbuf[7]));
                    else
                        well[w].kh.push_back(-1.0);

                    if (vbuf[8] != "DEFAULT")
                        well[w].skinFactor.push_back(stod(vbuf[8]));
                    else
                        well[w].skinFactor.push_back(0.0);

                    if (vbuf[9] != "DEFAULT")
                        well[w].direction.push_back(vbuf[9]);
                    else
                        well[w].direction.push_back("z");
                }
            }
        }
    }
}


/// Find pointer to the specified variable.
vector<OCP_DBL>* PreParamGridWell::FindPtr(const string& varName)
{
    vector<OCP_DBL>* myPtr = nullptr;

    switch (Map_Str2Int(&varName[0], varName.size())) {
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

    case Map_Str2Int("COORD", 5):
        gridType = CORNER_GRID;
        coord.reserve((nx + 1) * (ny + 1) * 6);
        myPtr = &coord;
        break;

    case Map_Str2Int("ZCORN", 5):
        gridType = CORNER_GRID;
        zcorn.reserve(numGrid * 8);
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

    case Map_Str2Int("TOPS", 4):
        tops.reserve(nx * ny);
        myPtr = &tops;
        break;
    }

    return myPtr;
}


template <typename T>
void PreParamGridWell::setVal(vector<T>& obj, const T& val, const vector<USI>& index)
{
    const OCP_USI nxny = nx * ny;
    OCP_USI id = 0;

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

    const OCP_USI nxny = nx * ny;
    OCP_USI id = 0;

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

    const OCP_USI nxny = nx * ny;
    OCP_USI id = 0;

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
// Generate active grids' connections
/////////////////////////////////////////////////////////////////////

void PreParamGridWell::Setup()
{
    SetupGrid();
    SetupConnWellGrid();
}


void PreParamGridWell::SetupGrid()
{
    switch (gridType) 
    {
    case ORTHOGONAL_GRID:
        SetupOrthogonalGrid();
        break;
    case CORNER_GRID:
        SetupCornerGrid();
        break;
    default:
        OCP_ABORT("WRONG Grid Type!");
    }

    OutputBaiscInfo();

}


void PreParamGridWell::SetupInitGrid()
{
    switch (gridType) {
    case ORTHOGONAL_GRID:
        SetupOrthogonalGrid();
        break;
    case CORNER_GRID:
        SetupCornerGrid();
        break;
    default:
        OCP_ABORT("WRONG Grid Type!");
    }

    OutputBaiscInfo();
}

void PreParamGridWell::SetupOrthogonalGrid()
{
    // x -> y -> z
    CalDepthVOrthogonalGrid();
    CalActiveGrid(1E-6, 1E-6);
    SetupActiveConnOrthogonalGrid();
}

void PreParamGridWell::CalDepthVOrthogonalGrid()
{
    depth.resize(numGrid, 0);
    const OCP_USI nxny = nx * ny;
    // 0th layer
    for (USI j = 0; j < ny; j++) {
        for (USI i = 0; i < nx; i++) {
            OCP_USI id = j * nx + i;
            depth[id] = tops[id] + dz[id] / 2;
        }
    }
    // 1th - (nz-1)th layer
    for (USI k = 1; k < nz; k++) {
        OCP_USI knxny = k * nxny;
        for (USI j = 0; j < ny; j++) {
            for (USI i = 0; i < nx; i++) {
                OCP_USI id = knxny + j * nx + i;
                depth[id] = depth[id - nxny] + dz[id - nxny] / 2 + dz[id] / 2;
            }
        }
    }

    v.resize(numGrid);
    for (OCP_USI i = 0; i < numGrid; i++) v[i] = dx[i] * dy[i] * dz[i];
}

void PreParamGridWell::SetupActiveConnOrthogonalGrid()
{
    gNeighbor.resize(activeGridNum);
    // PreAllocate
    for (OCP_USI n = 0; n < activeGridNum; n++) {
        gNeighbor[n].reserve(6);
    }

    // Begin Id and End Id in Grid, bIdg < eIdg
    OCP_USI       bIdg, eIdg, bIdb, eIdb;
    OCP_DBL       areaB, areaE;
    USI           direction;
    const OCP_USI nxny = nx * ny;

    for (USI k = 0; k < nz; k++) {
        for (USI j = 0; j < ny; j++) {
            for (USI i = 0; i < nx; i++) {

                bIdg = k * nxny + j * nx + i;

                if (!map_All2Act[bIdg].IsAct()) {
                    continue;
                }
                bIdb = map_All2Act[bIdg].GetId();

                // right  --  x-direction
                if (i < nx - 1) {
                    eIdg = bIdg + 1;
                    if (map_All2Act[eIdg].IsAct()) {
                        eIdb = map_All2Act[eIdg].GetId();

                        direction = 1;
                        areaB = 2 * dy[bIdg] * dz[bIdg] / dx[bIdg];
                        areaE = 2 * dy[eIdg] * dz[eIdg] / dx[eIdg];
                        gNeighbor[bIdb].push_back(GPair(eIdb, WEIGHT_GG, direction, areaB, areaE));
                        gNeighbor[eIdb].push_back(GPair(bIdb, WEIGHT_GG, direction, areaE, areaB));
                    }
                }
                // front  --  y-direction
                if (j < ny - 1) {
                    eIdg = bIdg + nx;
                    if (map_All2Act[eIdg].IsAct()) {
                        eIdb = map_All2Act[eIdg].GetId();

                        direction = 2;
                        areaB = 2 * dz[bIdg] * dx[bIdg] / dy[bIdg];
                        areaE = 2 * dz[eIdg] * dx[eIdg] / dy[eIdg];
                        gNeighbor[bIdb].push_back(GPair(eIdb, WEIGHT_GG, direction, areaB, areaE));
                        gNeighbor[eIdb].push_back(GPair(bIdb, WEIGHT_GG, direction, areaE, areaB));
                    }
                }
                // down --   z-direction
                if (k < nz - 1) {
                    eIdg = bIdg + nxny;
                    if (map_All2Act[eIdg].IsAct()) {
                        eIdb = map_All2Act[eIdg].GetId();

                        direction = 3;
                        areaB = 2 * dx[bIdg] * dy[bIdg] / dz[bIdg];
                        areaE = 2 * dx[eIdg] * dy[eIdg] / dz[eIdg];
                        gNeighbor[bIdb].push_back(GPair(eIdb, WEIGHT_GG, direction, areaB, areaE));
                        gNeighbor[eIdb].push_back(GPair(bIdb, WEIGHT_GG, direction, areaE, areaB));
                    }
                }
            }
        }
    }
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
}

void PreParamGridWell::SetupBasicCornerGrid(const OCP_COORD& CoTmp)
{
    dx = CoTmp.dx;
    dy = CoTmp.dy;
    dz = CoTmp.dz;
    v = CoTmp.v;
    depth = CoTmp.depth;
}

void PreParamGridWell::SetupActiveConnCornerGrid(const OCP_COORD& CoTmp)
{
    gNeighbor.resize(activeGridNum);
    // PreAllocate
    for (OCP_USI n = 0; n < activeGridNum; n++) {
        gNeighbor[n].reserve(10);
    }

    OCP_USI bIdg, eIdg, bIdb, eIdb;
    OCP_DBL areaB, areaE;
    USI     direction;
    for (OCP_USI n = 0; n < CoTmp.numConn; n++) {
        const GeneralConnect& ConnTmp = CoTmp.connect[n];

        bIdg = ConnTmp.begin;
        eIdg = ConnTmp.end;

        if (map_All2Act[bIdg].IsAct() && map_All2Act[eIdg].IsAct()) {
            bIdb = map_All2Act[bIdg].GetId();
            eIdb = map_All2Act[eIdg].GetId();

            direction = ConnTmp.directionType;
            areaB = ConnTmp.Ad_dd_begin;
            areaE = ConnTmp.Ad_dd_end;
            gNeighbor[bIdb].push_back(GPair(eIdb, WEIGHT_GG, direction, areaB, areaE));
            gNeighbor[eIdb].push_back(GPair(bIdb, WEIGHT_GG, direction, areaE, areaB));
        }
    }
}


void PreParamGridWell::CalActiveGrid(const OCP_DBL& e1, const OCP_DBL& e2)
{
    switch (model)
    {
    case ISOTHERMALMODEL:
        CalActiveGridIsoT(e1, e2);
        break;
    case THERMALMODEL:
        CalActiveGridT(e1, e2);
        break;
    default:
        OCP_ABORT("WRONG Grid Model!");
    }
}


/// If porosity or volume of the grid cell is too small, then the cell is inactive.
//  Note: Inactive cells do NOT participate simumlation; other rules can be given.
void PreParamGridWell::CalActiveGridIsoT(const OCP_DBL& e1, const OCP_DBL& e2)
{
    map_Act2All.reserve(numGrid);
    map_All2Act.resize(numGrid);
    OCP_USI count = 0;
    for (OCP_USI n = 0; n < numGrid; n++) {
        if (ACTNUM[n] == 0 || poro[n] < e1 || v[n] < e2) {
            map_All2Act[n] = GB_Pair(OCP_FALSE, 0);
            ACTNUM[n] = 0;
            continue;
        }
        map_Act2All.push_back(n);
        map_All2Act[n] = GB_Pair(OCP_TRUE, count);
        count++;
    }
    activeGridNum = count;
	cout << "  Number of inactive cells is " << (numGrid - activeGridNum) << " ("
		<< (numGrid - activeGridNum) * 100.0 / numGrid << "%)" << endl;

    // fluid grid = active grid
    fluidGridNum = activeGridNum;
    map_All2Flu = map_All2Act;
}

void PreParamGridWell::CalActiveGridT(const OCP_DBL& e1, const OCP_DBL& e2)
{
    map_Act2All.reserve(numGrid);
    map_All2Act.resize(numGrid);
    map_All2Flu.resize(numGrid);
    OCP_USI activeCount = 0;
    OCP_USI fluidCount = 0;
    for (OCP_USI n = 0; n < numGrid; n++) {
        if (ACTNUM[n] == 0 || v[n] < e1) {
            map_All2Act[n] = GB_Pair(OCP_FALSE, 0);
            map_All2Flu[n] = GB_Pair(OCP_FALSE, 0);
            ACTNUM[n] = 0;
        }
        else {
            if (poro[n] < e2) {
                map_All2Flu[n] = GB_Pair(OCP_FALSE, 0);
            }
            else {
                map_All2Flu[n] = GB_Pair(OCP_TRUE, fluidCount);
                fluidCount++;
            }
            map_Act2All.push_back(n);
            map_All2Act[n] = GB_Pair(OCP_TRUE, activeCount);
            activeCount++;
        }
    }
    activeGridNum = activeCount;
    fluidGridNum = fluidCount;

    cout << "  Number of inactive cells is " << (numGrid - activeGridNum) << " ("
        << (numGrid - activeGridNum) * 100.0 / numGrid << "%)" << endl;
}

/////////////////////////////////////////////////////////////////////
// Generate connections between active grids and wells
/////////////////////////////////////////////////////////////////////


void PreParamGridWell::SetupConnWellGrid()
{

    // Attention that all wells should be active -- own at least one connections to active grid
    numWell = well.size();
    connWellGrid.resize(numWell);
    for (USI w = 0; w < numWell; w++) {
        const USI numPerf = well[w].I_perf.size();
        for (USI p = 0; p < numPerf; p++) {
            const OCP_USI pId = (well[w].K_perf[p] - 1) * (nx * ny) + (well[w].J_perf[p] - 1) * nx + (well[w].I_perf[p] - 1);
            if (map_All2Flu[pId].IsAct()) {
                connWellGrid[w].push_back(map_All2Act[pId].GetId());
                gNeighbor[map_All2Act[pId].GetId()].push_back(GPair(w + activeGridNum, WEIGHT_GW, 0, 0, 0));
            }
        }
        if (connWellGrid[w].empty()) {
            OCP_ABORT("All perforations of Well " + to_string(w) + " are in Inactive grid!");
        }
    }


    const OCP_USI numTotal = activeGridNum + numWell;
    // Add well-grid connections to gNeighbor
    gNeighbor.resize(numTotal);
    for (USI w = 0; w < numWell; w++) {
        for (const auto& p : connWellGrid[w]) {
            gNeighbor[w + activeGridNum].push_back(GPair(p, WEIGHT_GW, 0, 0, 0));
        }
    }

    numNeighbor.resize(numTotal);
    for (OCP_USI n = 0; n < numTotal; n++) {
        numNeighbor[n] = gNeighbor[n].size();
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

    for (OCP_USI n = 0; n < numGrid; n++) {
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
    vector<OCP_DBL>().swap(tops);
    vector<OCP_DBL>().swap(coord);
    vector<OCP_DBL>().swap(zcorn);
    vector<OCP_DBL>().swap(ntg);
    vector<OCP_DBL>().swap(poro);
    vector<USI>().swap(ACTNUM);
    vector<PreParamWell>().swap(well);

    vector<OCP_DBL>().swap(v);
    vector<OCP_DBL>().swap(depth);
    vector<vector<GPair>>().swap(gNeighbor);
    vector<USI>().swap(numNeighbor);
    vector<OCP_USI>().swap(map_Act2All);
    vector<GB_Pair>().swap(map_All2Act);
    vector<GB_Pair>().swap(map_All2Flu);

    vector<vector<OCP_USI>>().swap(connWellGrid);
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Feb/15/2023      Create file                          */
/*----------------------------------------------------------------------------*/