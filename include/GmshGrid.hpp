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

#ifdef WITH_GMSH

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
    Edge(const OCP_USI& b, const OCP_USI& e) {
        if (b < e) {
            bId = b;
            eId = e;
        }
        else {
            bId = e;
            eId = b;
        }
    }
    Edge(const OCP_USI& b, const OCP_USI& e, const OCP_USI& tag_in) : Edge(b, e) {
        tag = tag_in;
    }
    Edge(const OCP_USI& b, const OCP_USI& e, const OCP_USI& tag_in, const OCP_USI& face, const OCP_USI& localEdge) : Edge(b, e) {
        tag = tag_in;
        faceIndex.push_back(face);
        faceIndex.push_back(localEdge);
    }
    Edge(const OCP_USI& b, const OCP_USI& e, const OCP_USI& tag_in, const string& phyinfo) : Edge(b, e) {
        tag      = tag_in;
        physical = phyinfo;
    }
    OCP_BOOL operator <(const Edge& e) const {
        if (bId < e.bId || (bId == e.bId && eId < e.eId)) return OCP_TRUE;
        else                                              return OCP_FALSE;
    }

public:
    /// index of begin node, bId < eId
    OCP_USI                 bId;
    /// index of end node
    OCP_USI                 eId;
    /// tag of edge(for debug)
    OCP_USI                 tag;
    /// physical info
    string                  physical;
    /// index of connected face and local index of edge
    mutable vector<OCP_USI> faceIndex;
    /// effective area from bId
    mutable vector<OCP_DBL> area;
};


class Polygon
{
public:
    /// constructor
    Polygon(const vector<OCP_USI>& pIndex, const OCP_USI& tag_in, const string& phyinfo);
    /// Calculate the center
    void CalCenter(const vector<OCP_DBL>& points);
    /// Calculate the area
    void CalArea(const vector<OCP_DBL>& points);
    /// judge if a point is in the element
    OCP_BOOL IfPointInElement(const Point3D& objP, const vector<OCP_DBL>& points);

public:   
    /// index of points(Store in order or reverse order)
    vector<OCP_USI> p;
    /// tag of face(for debug)
    OCP_USI         tag;
    /// physical info
    string          physical;
    /// location
    string          location;
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

protected:
    void InputGrid2D(const string& file);
    void Setup();
    void CalAreaCenter2D();
    void SetupConnAreaAndBoundary2D();

public:
    /// dimension
    USI             dimen;
    /// coordinates of points 
    vector<OCP_DBL> points;
    /// edges (for 2d now)
    set<Edge>       edges;
    /// elements (for 2d now)
    vector<Polygon> elements;
    /// Facies
    vector<Facies>  facies;
    /// element region index(facies)
    vector<INT>     faciesNum;
    /// map from facies (the order in which appear in GMSHPRO) to facies(the order in which appear in grids)
    vector<USI>     mapF2F;
    /// thickness (for 2d now)
    OCP_DBL         thickness;
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
