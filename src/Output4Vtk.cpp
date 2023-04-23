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

void Output4Vtk::InitASCII(const string&          myFile,
                           const string&          shortInfo,
                           const vector<OCP_DBL>& points_xyz) const
{
    ofstream outVtk(myFile);
    outVtk << VTK_HEADER << "\n";
    outVtk << shortInfo  << "\n";
    outVtk << VTK_ASCII  << "\n";
    outVtk << VTK_DATASET << " " << VTK_UNSTRUCTURED_GRID << "\n\n";
    // Output points
    const OCP_USI numGrid = points_xyz.size() / (3 * 8);
    outVtk << VTK_POINTS << " " << numGrid * 8 << " " << VTK_FLOAT << "\n";
    OCP_USI iterP = 0;
    for (OCP_USI n = 0; n < numGrid * 8; n++) {
        outVtk << setw(6) << points_xyz[iterP]
            << setw(10) << points_xyz[iterP + 1]
            << setw(10) << points_xyz[iterP + 2] << "\n";
        iterP += 3;
    }
    
    // Output cells
    outVtk << "\n" << VTK_CELLS << " " << numGrid << " " << numGrid * 9 << "\n";
    OCP_USI iterC = 0;
    for (OCP_USI n = 0; n < numGrid; n++) {
        outVtk << 8;
        for (OCP_USI i = 0; i < 8; i++)
            outVtk << setw(8) << iterC++;
        outVtk << "\n";
    }
    // OutPut cell types
    outVtk << "\n" << VTK_CELL_TYPES << " " << numGrid << "\n";
    for (OCP_USI n = 0; n < numGrid; n++)
        outVtk << VTK_HEXAHEDRON << "\n";
    outVtk.close();
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/19/2022      Create file                          */
/*  Chensong Zhang      Feb/05/2023      Update output in vtk files           */
/*----------------------------------------------------------------------------*/
