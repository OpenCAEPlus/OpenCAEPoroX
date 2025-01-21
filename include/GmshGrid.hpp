/*! \file    GmshGrid.hpp
 *  \brief   GmshGrid class declaration
 *  \author  Shizhe Li
 *  \date    Sep/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifdef OCP_USE_GMSH

#ifndef __GMSHGRID_HEADER__
#define __GMSHGRID_HEADER__


// Standard header files
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <set>

// Gmsh header files
#include <gmsh.h>

// OpenCAEPoroX header files
#include "UtilInput.hpp"
#include "UtilMesh.hpp"

using namespace std;

class Edge
{
public:
    Edge(const OCP_ULL& b, const OCP_ULL& e) {
        if (b < e) {
            bId = b;
            eId = e;
        }
        else {
            bId = e;
            eId = b;
        }
    }
    Edge(const OCP_ULL& b, const OCP_ULL& e, const OCP_ULL& tag_in, const OCP_ULL& face, const OCP_ULL& localEdge) : Edge(b, e) {
        tag = tag_in;
        faceIndex.push_back(face);
        faceIndex.push_back(localEdge);
    }
    Edge(const OCP_ULL& b, const OCP_ULL& e, const OCP_ULL& tag_in, const string& physical_in, const USI& phyIndex_in) : Edge(b, e) {
        tag      = tag_in;
        physical = physical_in;
        phyIndex = pow(2, phyIndex_in);
    }
    auto operator <(const Edge& e) const {
        if (bId < e.bId || (bId == e.bId && eId < e.eId)) return OCP_TRUE;
        else                                              return OCP_FALSE;
    }

public:
    /// index of begin node, bId < eId
    OCP_ULL                 bId;
    /// index of end node
    OCP_ULL                 eId;
    /// tag of edge(for debug)
    OCP_ULL                 tag;
    /// physical info(for debug)
    string                  physical;
    /// physical info binary encoding
    OCP_USI                 phyIndex{ 0 };
    /// index of connected face and local index of edge
    mutable vector<OCP_ULL> faceIndex;
    /// effective area from bId
    mutable vector<OCP_DBL> area;
};


class OCP_Polygon
{
public:
    /// constructor
    OCP_Polygon(const vector<OCP_ULL>& pIndex, const OCP_ULL& tag_in, const string& phyinfo, const OCP_USI& index);
    /// Calculate the center
    void CalCenter(const vector<OCP_DBL>& points);
    /// Calculate the area
    void CalArea(const vector<OCP_DBL>& points);
    /// judge if a point is in the element
    OCP_BOOL IfPointInElement(const Point3D& objP, const vector<OCP_DBL>& points);

public:   
    /// index of points(Store in order or reverse order)
    vector<OCP_ULL> p;
    /// tag of face (for debug)
    OCP_ULL         tag;
    /// physical info (for debug)
    string          physical;
    /// physical Index
    OCP_USI         phyIndex;
    /// boundary (for debug)
    string          boundary;
    /// index of boundary (binary encoding)
    OCP_USI         boundIndex{ 0 };
    /// effective area of boundary
    OCP_DBL         boundArea{ 0 };
    /// center
    Point3D         center;
    /// area
    OCP_DBL         area;
};

/// Facies property
class Facies
{
public:
    Facies(const string& faciesname) : name(faciesname) {};
public:
    /// name of Facies
    string  name;
    /// porosity of current facies
    OCP_DBL poro{ -1 };
    /// permeability of current facies
    OCP_DBL kx{ -1 };
    OCP_DBL ky{ -1 };
    OCP_DBL kz{ -1 };
};


class GMSHGrid
{
public:
    /// for 2-dimension now
    void InputGrid(const string& file);
    /// input property for each region
    void InputProperty(ifstream& ifs);
    /// get boundary name
    INT GetBoundaryName(vector<string>& names) const { 
        if (ifUse) {
            if (dimen == 2)      names = physicalNameSet[1];
            else if (dimen == 3) names = physicalNameSet[2];
            else                 OCP_ABORT("WRONG DIMENSION!");
            return names.size();
        }
        else {
            return 0;
        }
    }
    void PrintElementPoint(const OCP_ULL& n) const;

protected:
    void InputGrid2D(const string& file);
    void Setup();
    void CalAreaCenter2D();
    void SetupConnAreaAndBoundary2D();


public:
    /// If use the gmsh grid
    OCP_BOOL               ifUse{ OCP_FALSE };
    /// dimension
    USI                    dimen;
    /// coordinates of points 
    vector<OCP_DBL>        points;
    /// edges (for 2d now)
    set<Edge>              edges;
    /// elements (for 2d now)
    vector<OCP_Polygon>        elements;
    /// Facies
    vector<Facies>         facies;
    /// physical name set
    vector<vector<string>> physicalNameSet;
    /// thickness (for 2d now)
    OCP_DBL                thickness;
};


#endif /* end if __GMSHGRID_HEADER__ */


#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/
