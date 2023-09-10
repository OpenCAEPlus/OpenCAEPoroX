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
    Edge(const OCP_USI& b, const OCP_USI& e, const OCP_USI& tag_in, const OCP_USI& face) : Edge(b, e) {
        tag = tag_in;
        faceTag.push_back(face);
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
    OCP_USI         bId;
    /// index of end node
    OCP_USI         eId;
    /// tag of edge(for debug)
    OCP_USI         tag;
    /// tag of connected face
    mutable vector<OCP_USI> faceTag;
    /// physical info
    string          physical;
};


class Polygon
{
public:
    Polygon(const double* p0, const double* p1, const double* p2, const OCP_USI& tag_in) {
        p.push_back(Point3D(*p0, *(p0 + 1), *(p0 + 2)));
        p.push_back(Point3D(*p1, *(p1 + 1), *(p1 + 2)));
        p.push_back(Point3D(*p2, *(p2 + 1), *(p2 + 2)));
        tag = tag_in;
    }
    Polygon(const double* p0, const double* p1, const double* p2, const double* p3, const OCP_USI& tag_in) {
        p.push_back(Point3D(*p0, *(p0 + 1), *(p0 + 2)));
        p.push_back(Point3D(*p1, *(p1 + 1), *(p1 + 2)));
        p.push_back(Point3D(*p2, *(p2 + 1), *(p2 + 2)));
        p.push_back(Point3D(*p3, *(p3 + 1), *(p3 + 2)));
        tag = tag_in;
    }
public:
    /// points(Store in order or reverse order)
    vector<Point3D> p;
    /// tag of face(for debug)
    OCP_USI         tag;
    /// physical info
    string          physical;
};



class GMSHGrid
{
public:
    /// for 2-dimension now
	void Input(const string& file);


protected:
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
