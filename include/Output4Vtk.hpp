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
const string  VTK_BINARY            = "BINARY";
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


template <typename T>
void SwapEnd(T* var, const OCP_ULL len)
{
    for (OCP_ULL n = 0; n < len; n++) {
        char* varArray = reinterpret_cast<char*>(&var[n]);
        for (long i = 0; i < static_cast<long>(sizeof(var[n]) / 2); i++)
            std::swap(varArray[sizeof(var[n]) - 1 - i], varArray[i]);
    }
}


class Output4Vtk
{
public:
    /// Setup
    void Setup(const OCP_BOOL& ascii) { ifASCII = ascii; }
    /// create a new file and write common information
    OCP_ULL Init(const string& dir, const string& myFile, const string& shortInfo) const;
    template <typename T>
    void OutputCELL_DATA_SCALARS(ofstream&        outVtk,
                                 const string&    dataName,
                                 const string&    dataType,
                                 const vector<T>& tmpV,
                                 const OCP_ULL&   bId,
                                 const OCP_ULL&   nb,
                                 const USI&       digits) const;

public:
    static void OutputGridInfo(const string& dir, const OCP_ULL& nG, const vector<OCP_SIN>& points_xyz,
                               const vector<OCP_ULL>& cell_points, const vector<USI>& cell_type);

protected:
    OCP_ULL InitASCII(const string& dir, const string& myFile, const string& shortInfo) const;
    OCP_ULL InitBINARY(const string& dir, const string& myFile, const string& shortInfo) const;
    void InputGridInfo(const string& dir, OCP_ULL& nG, OCP_ULL& nP, vector<OCP_SIN>& points_xyz, vector<OCP_ULL>& cell_points, vector<USI>& cell_type) const;

protected:
    OCP_BOOL            ifASCII{ OCP_FALSE };
    static const string tmpFile;
};


template <typename T>
void Output4Vtk::OutputCELL_DATA_SCALARS(ofstream&        outVtk,
                                         const string&    dataName,
                                         const string&    dataType,
                                         const vector<T>& tmpV,
                                         const OCP_ULL&   bId,
                                         const OCP_ULL&   nb,
                                         const USI&       digits) const
{
    outVtk << "\n" << VTK_SCALARS << " " << dataName << " " << dataType << " " << 1;
    outVtk << "\n" << VTK_LOOKUP_TABLE << " " << VTK_DEFAULT << "\n";

    if (ifASCII) {
        outVtk << fixed << setprecision(digits);
        for (OCP_ULL n = 0; n < nb; n++) {
            outVtk << tmpV[bId + n] << "\n";
        }
    }
    else {
        T* tmp = const_cast<T*>(&tmpV[bId]);
        SwapEnd(tmp, nb);
        outVtk.write((const char*)tmp, nb * sizeof(tmp[0]));
        outVtk << "\n";
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
