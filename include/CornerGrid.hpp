/*! \file    CornerGrid.hpp
 *  \brief   Declaration of classes related to the corner grid
 *  \author  Shizhe Li
 *  \date    Nov/16/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __CORNERGRID_HEADER__
#define __CORNERGRID_HEADER__

// OpenCAEPoroX header files
#include "UtilMesh.hpp"

using namespace std;

// Constants used by corner-point grid code
const OCP_DBL TEENY        = 1E-3;  ///< Used for checking distance b/w center to face
const USI     MAX_NEIGHBOR = 80;    ///< Max number of  neighbors allowed


/// half-connections
class HalfConn
{
public:
    OCP_DBL    Ad_dd;
    Point3D    d;
    OCP_ULL    neigh;
    ConnDirect directionType; 
};

/// all connections
class ConnGrid
{
public:
    USI              nConn, maxConn;
    vector<HalfConn> halfConn;
    void             Allocate(const USI& max_neighbor);
    void             AddHalfConn(const OCP_ULL& n,
                                 const Point3D& area,
                                 const Point3D& d,
                                 const ConnDirect&     direction,
                                 const OCP_DBL& flag = 1);
};

/// General connection
class GeneralConnect
{
public:
    OCP_ULL    begin, end;
    ConnDirect directionType;
    OCP_DBL    Ad_dd_begin;
    OCP_DBL    Ad_dd_end;
};

/// corner grid
class OCP_COORD
{
    friend class PreParamGridWell;

public:
    void     Allocate(const USI& Nx, const USI& Ny, const USI& Nz);
    void     InputData(const vector<OCP_DBL>& coord, const vector<OCP_DBL>& zcorn);
    OCP_BOOL InputCOORDDATA(const vector<OCP_DBL>& coord);
    OCP_BOOL InputZCORNDATA(const vector<OCP_DBL>& zcorn);
    // New version
    void SetupCornerPoints();
    void SetAllFlags(const HexahedronFace& oFace, const HexahedronFace& Face);
    // functions
    OCP_DBL OCP_SIGN(const OCP_DBL& x) { return x >= 0 ? 1 : -1; }

private:
    USI           nx;
    USI           ny;
    USI           nz;
    OCP_DBL***    COORDDATA;
    OCP_DBL****   ZCORNDATA;
    Hexahedron*** cornerPoints;

    OCP_ULL         numGrid;
    OCP_ULL         numConn;
    OCP_ULL         numConnMax;
    vector<OCP_DBL> v;
    vector<OCP_DBL> depth;
    vector<OCP_DBL> dx;
    vector<OCP_DBL> dy;
    vector<OCP_DBL> dz;
    vector<Point3D> center;

    vector<GeneralConnect> connect;

    // Auxiliary variables
    OCP_BOOL       flagQuad;
    OCP_BOOL       upNNC, downNNC;
    OCP_BOOL       flagJump;
    HexahedronFace interFace;
    // after the Axes are determined, blocks will be placed along the y+, or along the
    // y-. if y+, then flagForward equals 1.0, else -1.0, this relates to calculation of
    // area normal vector
    OCP_DBL flagForward;
};

#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/16/2021      Create file                          */
/*  Chensong Zhang      Jan/16/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/