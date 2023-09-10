/*! \file    UtilMesh.hpp
 *  \brief   UtilMes class declaration
 *  \author  Shizhe Li
 *  \date    Sep/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __UTILMESH_HEADER__
#define __UTILMESH_HEADER__


// Standard header files
#include <math.h>
#include <stdlib.h>
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"

using namespace std;


/// A point in 2D.
class Point2D
{
public:
    OCP_DBL x;
    OCP_DBL y;

public:
    Point2D() = default;
    Point2D(OCP_DBL x0, OCP_DBL y0)
        : x(x0)
        , y(y0) {};
};

/// A point in 3D.
class Point3D
{
public:
    OCP_DBL x;
    OCP_DBL y;
    OCP_DBL z;

public:
    Point3D() = default;
    Point3D(OCP_DBL x0, OCP_DBL y0, OCP_DBL z0)
        : x(x0)
        , y(y0)
        , z(z0) {};

    Point3D& operator=(const Point3D& other);       ///< equal
    Point3D  operator+(const Point3D& other) const; ///< Addition
    Point3D  operator-(const Point3D& other) const; ///< Subtraction
    OCP_DBL  operator*(const Point3D& other) const; ///< Multiplication
    Point3D& operator+=(const Point3D& other);
    Point3D& operator*=(const OCP_DBL& a);
    Point3D& operator/=(const OCP_DBL& a);
    void     Reset()
    {
        x = 0;
        y = 0;
        z = 0;
    };
};

Point3D operator*(const Point3D& p, const OCP_DBL& a);      ///< Point * a
Point3D operator*(const OCP_DBL& a, const Point3D& p);      ///< a * Point
Point3D CrossProduct(const Point3D& p1, const Point3D& p2); ///< Cross product

/// A hexahedron cell.
class Hexahedron
{
public:
    /// Calculate volume of hexahedron
    OCP_DBL CalVolum() const;
    /// Calculate center of hexahedron
    Point3D CalCenter() const;

public:
    Point3D p0, p1, p2, p3, p4, p5, p6, p7;
};

/// A face of a hexahedron cell.
class HexahedronFace
{
public:
    /// Find the area normal vector of a face.
    Point3D CalAreaVector() const;
    /// Find the center of a face.
    Point3D CalCenter() const;
public:  
    Point3D p0, p1, p2, p3;
};

/// 3 by 3 matrix.
class Matrix3
{
public:
    OCP_DBL M[3][3];
    Point3D operator*(const Point3D& v) const;
};


/// Calculate the cross point of two lines
Point2D CalCrossingPoint(const Point2D Line1[2], const Point2D Line2[2]);

/// Calculate the area of interface of two HexahedronFace(the interface is not quadrilateral)
OCP_DBL CalAreaNotQuadr(const HexahedronFace& FACE1, const HexahedronFace& FACE2);





#endif /* end if __UTILMESH_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/