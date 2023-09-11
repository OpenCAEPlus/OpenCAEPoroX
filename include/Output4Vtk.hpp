/*! \file    Output4Vtk.hpp
 *  \brief   Output reservoir information in vtk format
 *  \author  Shizhe Li
 *  \date    Oct/19/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OUTPUT4VTK_HEADER__
#define __OUTPUT4VTK_HEADER__

#include <fstream>
#include <string>
#include <vector>

using namespace std;

#include "OCPConst.hpp"

// Basic datatype
typedef OCP_DBL VTK_DBL;
typedef OCP_SIN VTK_SIN;
typedef OCP_USI VTK_USI;
typedef OCP_ULL VTK_ULL;

// Basic Keyword
const string  VTK_HEADER            = "# vtk DataFile Version 3.0";
const string  VTK_ASCII             = "ASCII";
const string  VTK_DATASET           = "DATASET";
const string  VTK_UNSTRUCTURED_GRID = "UNSTRUCTURED_GRID";
const string  VTK_POINTS            = "POINTS";
const string  VTK_CELLS             = "CELLS";
const string  VTK_CELL_TYPES        = "CELL_TYPES";
const string  VTK_CELL_DATA         = "CELL_DATA";
const string  VTK_POINT_DATA        = "POINT_DATA";
const VTK_USI VTK_MAX_TITLE_LENGTH  = 256;
const string  VTK_LOOKUP_TABLE      = "LOOKUP_TABLE";
const string  VTK_DEFAULT           = "default";
const string  VTK_SCALARS           = "SCALARS";

// Basic Cell Type
const VTK_USI VTK_POLY_LINE  = 4;
const VTK_USI VTK_TRIANGLE   = 5;
const VTK_USI VTK_QUAD       = 9;
const VTK_USI VTK_HEXAHEDRON = 12;


const string VTK_FLOAT        = "float";
const string VTK_UNSIGNED_INT = "unsigned_int";

class Output4Vtk
{
    friend class Out4VTK;

public:
    /// create a new file and write common information
    OCP_USI InitASCII(const string& dir, 
                      const string& myFile,
                      const string& shortInfo) const;
    template <typename T>
    void OutputCELL_DATA_SCALARS(ofstream&        outVtk,
                                 const string&    dataName,
                                 const string&    dataType,
                                 const vector<T>  tmpV,
                                 const OCP_USI&   bId,
                                 const OCP_USI&   nb,
                                 const USI&       digits) const;

public:
    static void OutputGridInfo(const string& dir, const OCP_USI& nG, const OCP_USI& nP, const vector<OCP_DBL>& points_xyz,
                               const vector<OCP_USI>& cell_points, const vector<USI>& cell_type);

protected:
    void InputGridInfo(const string& dir, OCP_USI& nG, OCP_USI& nP, vector<OCP_DBL>& points_xyz, vector<OCP_USI>& cell_points, vector<USI>& cell_type) const;

protected:
    static const string tmpFile;
};


template <typename T>
void Output4Vtk::OutputCELL_DATA_SCALARS(ofstream&        outVtk,
                                         const string&    dataName,
                                         const string&    dataType,
                                         const vector<T>  tmpV,                        
                                         const OCP_USI&   bId,
                                         const OCP_USI&   nb, 
                                         const USI&       digits) const
{
    outVtk << "\n" << VTK_SCALARS << " " << dataName << " " << dataType << " " << 1;
    outVtk << "\n" << VTK_LOOKUP_TABLE << " " << VTK_DEFAULT << "\n";
    outVtk << fixed << setprecision(digits);
    for (OCP_USI n = 0; n < nb; n++) {
        outVtk << tmpV[bId + n] << "\n";
    }
}

#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/19/2022      Create file                          */
/*  Chensong Zhang      Feb/05/2023      Update output in vtk files           */
/*----------------------------------------------------------------------------*/
