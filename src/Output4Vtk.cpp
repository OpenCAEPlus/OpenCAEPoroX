/*! \file    Output4Vtk.cpp
 *  \brief   Output reservoir information in vtk format
 *  \author  Shizhe Li
 *  \date    Oct/19/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "Output4Vtk.hpp"

const string Output4Vtk::tmpFile = "grid.tmpinfo";



OCP_ULL Output4Vtk::Init(const string& dir, const string& myFile, const string& shortInfo) const
{
    if (ifASCII) {
        return InitASCII(dir, myFile, shortInfo);
    }
    else {
        return InitBINARY(dir, myFile, shortInfo);
    }
}


OCP_ULL Output4Vtk::InitASCII(const string& dir,
                              const string& myFile,
                              const string& shortInfo) const
{

    // Input Grid info from disk, which was stored before
    OCP_ULL nG, nP;
    vector<OCP_SIN> points_xyz;  ///< x,y,z coordinates
    vector<OCP_ULL> cell_points; ///< nP, points index
    vector<USI>     cell_type;   ///< type of cell

    InputGridInfo(dir, nG, nP, points_xyz, cell_points, cell_type);

    // output Grid info

    ofstream outVtk(myFile);
    outVtk << VTK_HEADER << "\n";
    outVtk << shortInfo  << "\n";
    outVtk << VTK_ASCII  << "\n";
    outVtk << VTK_DATASET << " " << VTK_UNSTRUCTURED_GRID << "\n";
    // Output points
    outVtk << VTK_POINTS << " " << nP << " " << VTK_FLOAT << "\n";
    OCP_ULL iterP = 0;
    for (OCP_ULL n = 0; n < nP; n++) {
        outVtk << setw(12) << points_xyz[iterP]
            << setw(12) << points_xyz[iterP + 1]
            << setw(12) << points_xyz[iterP + 2] << "\n";
        iterP += 3;
    }
    
    // Output cells
    outVtk << "\n" << VTK_CELLS << " " << nG << " " << cell_points.size() << "\n";
    OCP_ULL iterC = 0;
    for (OCP_ULL n = 0; n < nG; n++) {
        const USI nP = cell_points[iterC++];
        outVtk << nP;
        for (OCP_ULL i = 0; i < nP; i++)
            outVtk << setw(8) << cell_points[iterC++];
        outVtk << "\n";
    }
    // OutPut cell types
    outVtk << "\n" << VTK_CELL_TYPES << " " << nG << "\n";
    for (OCP_ULL n = 0; n < nG; n++)
        outVtk << cell_type[n] << "\n";
    outVtk.close();


    vector<OCP_SIN>().swap(points_xyz);
    vector<OCP_ULL>().swap(cell_points);
    vector<USI>().swap(cell_type);

    return nG;
}


OCP_ULL Output4Vtk::InitBINARY(const string& dir, const string& myFile, const string& shortInfo) const
{
    // Input Grid info from disk, which was stored before
    OCP_ULL nG, nP;
    vector<OCP_SIN> points_xyz;  ///< x,y,z coordinates
    vector<OCP_ULL> cell_points; ///< nP, points index
    vector<USI>     cell_type;   ///< type of cell

    InputGridInfo(dir, nG, nP, points_xyz, cell_points, cell_type);

    // output Grid info

    ofstream outVtk(myFile, ios::binary);
    outVtk << VTK_HEADER << "\n";
    outVtk << shortInfo << "\n";
    outVtk << VTK_BINARY << "\n";
    outVtk << VTK_DATASET << " " << VTK_UNSTRUCTURED_GRID << "\n";
    // Output points
    outVtk << VTK_POINTS << " " << nP << " " << VTK_FLOAT << "\n";
    SwapEnd(points_xyz.data(), points_xyz.size());
    outVtk.write((const char*)&points_xyz[0], points_xyz.size() * sizeof(points_xyz[0]));
    outVtk << "\n";

    // Output cells
    outVtk << VTK_CELLS << " " << nG << " " << cell_points.size() << "\n";
    SwapEnd(cell_points.data(), cell_points.size());
    outVtk.write((const char*)&cell_points[0], cell_points.size() * sizeof(cell_points[0]));
    outVtk << "\n";

    // OutPut cell types
    outVtk << VTK_CELL_TYPES << " " << nG << "\n";
    SwapEnd(cell_type.data(), cell_type.size());
    outVtk.write((const char*)&cell_type[0], cell_type.size() * sizeof(cell_type[0]));
    outVtk << "\n";

    outVtk.close();

    vector<OCP_SIN>().swap(points_xyz);
    vector<OCP_ULL>().swap(cell_points);
    vector<USI>().swap(cell_type);

    return nG;
}


void Output4Vtk::OutputGridInfo(const string& dir, const OCP_ULL& nG, const vector<OCP_SIN>& points_xyz,
    const vector<OCP_ULL>& cell_points, const vector<USI>& cell_type)
{
    const string myFile = dir + tmpFile;
    ofstream outF(myFile, ios::out | ios::binary);
    if (!outF.is_open()) {
        OCP_ABORT("Can not open " + myFile);
    }

    const OCP_ULL nP = points_xyz.size() / 3;

    outF.write((const char*)&nG, sizeof(nG));
    outF.write((const char*)&nP, sizeof(nP));
    outF.write((const char*)&points_xyz[0], points_xyz.size() * sizeof(points_xyz[0]));
    const OCP_ULL len_cell_points = cell_points.size();
    outF.write((const char*)&len_cell_points, sizeof(len_cell_points));
    outF.write((const char*)&cell_points[0], cell_points.size() * sizeof(cell_points[0]));
    outF.write((const char*)&cell_type[0], cell_type.size() * sizeof(cell_type[0]));
    outF.close();
}


void Output4Vtk::InputGridInfo(const string& dir, OCP_ULL& nG, OCP_ULL& nP, vector<OCP_SIN>& points_xyz, vector<OCP_ULL>& cell_points, vector<USI>& cell_type) const
{
    const string gridFile = dir + tmpFile;
    ifstream inP(gridFile, ios::in | ios::binary);
    if (!inP.is_open()) {
        OCP_WARNING("Can not open " + gridFile);
    }
    inP.read((OCP_CHAR*)(&nG), sizeof(nG));
    inP.read((OCP_CHAR*)(&nP), sizeof(nP));
    points_xyz.resize(nP * 3);
    inP.read((OCP_CHAR*)(&points_xyz[0]), sizeof(points_xyz[0]) * points_xyz.size());
    OCP_ULL len_cell_points;
    inP.read((OCP_CHAR*)(&len_cell_points), sizeof(len_cell_points));
    cell_points.resize(len_cell_points);
    inP.read((OCP_CHAR*)(&cell_points[0]), sizeof(cell_points[0]) * cell_points.size());
    cell_type.resize(nG);
    inP.read((OCP_CHAR*)(&cell_type[0]), sizeof(cell_type[0]) * cell_type.size());

    inP.close();
    if (remove(gridFile.c_str()) != 0) {
        OCP_WARNING("Failed to delete " + gridFile);
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/19/2022      Create file                          */
/*  Chensong Zhang      Feb/05/2023      Update output in vtk files           */
/*----------------------------------------------------------------------------*/
