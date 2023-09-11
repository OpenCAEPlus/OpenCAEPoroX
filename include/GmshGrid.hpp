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
            bId = b - 1;
            eId = e - 1;
        }
        else {
            bId = e - 1;
            eId = b - 1;
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
    Polygon(const vector<OCP_USI>& pIndex, const OCP_USI& tag_in, const string& phyinfo) {
        p.clear();
        for (USI n = 0; n < pIndex.size(); n++) {
            p.push_back(pIndex[n] - 1);
        }    
        tag      = tag_in;
        physical = phyinfo;
    }

    void CalCenter(const vector<OCP_DBL>& points) {
        const USI np = p.size();
        center.resize(3);
        for (USI i = 0; i < np; i++) {
            center[0] += points[3 * p[i] + 0];
            center[1] += points[3 * p[i] + 1];
            center[2] += points[3 * p[i] + 2];
        }
        center[0] /= np;
        center[1] /= np;
        center[2] /= np;
    }
    void CalArea(const vector<OCP_DBL>& points) {
        if (p.size() == 3) {
            const OCP_DBL* p0 = &points[3 * p[0]];
            const OCP_DBL* p1 = &points[3 * p[1]];
            const OCP_DBL* p2 = &points[3 * p[2]];
            const USI      x  = 0;
            const USI      y  = 1;
            area = 0.5 * abs((p2[x] - p1[x]) * (p0[y] - p1[y]) - (p2[y] - p1[y]) * (p0[x] - p1[x]));
        }
        else {
            const OCP_DBL* p0 = &points[3 * p[0]];
            const OCP_DBL* p1 = &points[3 * p[1]];
            const OCP_DBL* p2 = &points[3 * p[2]];
            const OCP_DBL* p3 = &points[3 * p[3]];
            const USI      x  = 0;
            const USI      y  = 1;
            area = 0.5 * (abs((p2[x] - p1[x]) * (p0[y] - p1[y]) - (p2[y] - p1[y]) * (p0[x] - p1[x])) 
                 + abs((p0[x] - p3[x]) * (p2[y] - p3[y]) - (p0[y] - p3[y]) * (p2[x] - p3[x])));
        }
    }

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
    vector<OCP_DBL> center;
    /// area
    OCP_DBL         area;
};


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
    OCP_DBL k{ -1 };
};


class GMSHGrid
{
public:
    /// for 2-dimension now
    void InputGrid(const string& file);
    /// input property for each region
    void InputProperty(const string& file);

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
