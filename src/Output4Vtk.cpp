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

OCP_USI Output4Vtk::InitASCII(const string& dir,
                              const string& myFile,
                              const string& shortInfo) const
{

    // Input Grid info from disk, which was stored before
    vector<OCP_DBL> points_xyz;  ///< x,y,z coordinates
    vector<OCP_USI> cell_points; ///< numpoints, points index
    vector<USI>     cell_type;   ///< type of cell
    OCP_USI         len_cell_points; ///< length of cell_points

    const string pointsFile = dir + "points.out";
    ifstream inP(pointsFile, ios::in | ios::binary);
    if (!inP.is_open()) {
        OCP_WARNING("Can not open " + pointsFile);
    }
    inP.read((OCP_CHAR*)(&numGrid), sizeof(numGrid));
    inP.read((OCP_CHAR*)(&numPoints), sizeof(numPoints));
    points_xyz.resize(numPoints * 3);
    inP.read((OCP_CHAR*)(&points_xyz[0]), sizeof(points_xyz[0]) * points_xyz.size());
    inP.read((OCP_CHAR*)(&len_cell_points), sizeof(len_cell_points));
    cell_points.resize(len_cell_points);
    inP.read((OCP_CHAR*)(&cell_points[0]), sizeof(cell_points[0]) * cell_points.size());
    cell_type.resize(numGrid);
    inP.read((OCP_CHAR*)(&cell_type[0]), sizeof(cell_type[0]) * cell_type.size());

    inP.close();
    if (remove(pointsFile.c_str()) != 0) {
        OCP_WARNING("Failed to delete " + pointsFile);
    }

    // output Grid info

    ofstream outVtk(myFile);
    outVtk << VTK_HEADER << "\n";
    outVtk << shortInfo  << "\n";
    outVtk << VTK_ASCII  << "\n";
    outVtk << VTK_DATASET << " " << VTK_UNSTRUCTURED_GRID << "\n\n";
    // Output points
    outVtk << VTK_POINTS << " " << numPoints << " " << VTK_FLOAT << "\n";
    OCP_USI iterP = 0;
    for (OCP_USI n = 0; n < numPoints; n++) {
        outVtk << setw(6) << points_xyz[iterP]
            << setw(10) << points_xyz[iterP + 1]
            << setw(10) << points_xyz[iterP + 2] << "\n";
        iterP += 3;
    }
    
    // Output cells
    outVtk << "\n" << VTK_CELLS << " " << numGrid << " " << len_cell_points << "\n";
    OCP_USI iterC = 0;
    for (OCP_USI n = 0; n < numGrid; n++) {
        const USI nP = cell_points[iterC++];
        outVtk << nP;
        for (OCP_USI i = 0; i < nP; i++)
            outVtk << setw(8) << cell_points[iterC++];
        outVtk << "\n";
    }
    // OutPut cell types
    outVtk << "\n" << VTK_CELL_TYPES << " " << numGrid << "\n";
    for (OCP_USI n = 0; n < numGrid; n++)
        outVtk << cell_type[n] << "\n";
    outVtk.close();


    vector<OCP_DBL>().swap(points_xyz);
    vector<OCP_USI>().swap(cell_points);
    vector<USI>().swap(cell_type);

    return numGrid;
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/19/2022      Create file                          */
/*  Chensong Zhang      Feb/05/2023      Update output in vtk files           */
/*----------------------------------------------------------------------------*/
