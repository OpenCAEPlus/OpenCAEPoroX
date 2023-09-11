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
#include <gmsh.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <set>

// OpenCAEPoroX header files
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
    OCP_USI                  bId;
    /// index of end node
    OCP_USI                  eId;
    /// tag of edge(for debug)
    OCP_USI                  tag;
    /// physical info
    string                   physical;
    /// index of connected face and local index of edge
    mutable vector<OCP_USI>  faceIndex;
    /// effective area from bId
    mutable vector <OCP_DBL> area;
};


class Polygon
{
public:
    Polygon(const double* p0, const double* p1, const double* p2, const OCP_USI& tag_in, const string& phyinfo) {
        p.push_back(Point3D(*p0, *(p0 + 1), *(p0 + 2)));
        p.push_back(Point3D(*p1, *(p1 + 1), *(p1 + 2)));
        p.push_back(Point3D(*p2, *(p2 + 1), *(p2 + 2)));
        tag      = tag_in;
        physical = phyinfo;
        CalCenter();
        CalArea();
    }
    Polygon(const double* p0, const double* p1, const double* p2, const double* p3, const OCP_USI& tag_in, const string& phyinfo) {
        p.push_back(Point3D(*p0, *(p0 + 1), *(p0 + 2)));
        p.push_back(Point3D(*p1, *(p1 + 1), *(p1 + 2)));
        p.push_back(Point3D(*p2, *(p2 + 1), *(p2 + 2)));
        p.push_back(Point3D(*p3, *(p3 + 1), *(p3 + 2)));
        tag      = tag_in;
        physical = phyinfo;
        CalCenter();
        CalArea();
    }

protected:
    void CalCenter() {
        center.Reset();
        for (USI i = 0; i < p.size(); i++) {
            center += p[i];
        }
        center /= p.size();
    }
    void CalArea() {
        if (p.size() == 3) {
            Point3D v0, v1;
            v0 = p[2] - p[1];
            v1 = p[0] - p[1];
            area = 0.5 * abs(v0.x * v1.y - v0.y * v1.x);
        }
        else {
            Point3D v0, v1, v2, v3;
            v0 = p[2] - p[1];
            v1 = p[0] - p[1];
            v2 = p[0] - p[3];
            v3 = p[2] - p[3];
            area = 0.5 * (abs(v0.x * v1.y - v0.y * v1.x) + abs(v2.x * v3.y - v2.y * v3.x));
        }
    }

public:   
    /// points(Store in order or reverse order)
    vector<Point3D> p;
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



class GMSHGrid
{
public:
    /// for 2-dimension now
    void Input(const string& file);

protected:
    void Input2D(const string& file);
    void Setup();	
    void SetupConnAreaAndBoundary2D();


protected:
    USI             dimen;
    set<Edge>       edges;
    vector<Polygon> elements;
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
