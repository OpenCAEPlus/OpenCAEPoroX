/*! \file    UtilMesh.cpp
 *  \brief   UtilMesh class declaration
 *  \author  Shizhe Li
 *  \date    Sep/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "UtilMesh.hpp"


Point3D& Point3D::operator=(const Point3D& other)
{
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
}

Point3D Point3D::operator+(const Point3D& other) const
{
    return Point3D(x + other.x, y + other.y, z + other.z);
}

Point3D Point3D::operator-(const Point3D& other) const
{
    return Point3D(x - other.x, y - other.y, z - other.z);
}

OCP_DBL Point3D::operator*(const Point3D& other) const
{
    return x * other.x + y * other.y + z * other.z;
}

Point3D& Point3D::operator+=(const Point3D& other)
{
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

Point3D& Point3D::operator*=(const OCP_DBL& a)
{
    x *= a;
    y *= a;
    z *= a;
    return *this;
}

Point3D& Point3D::operator/=(const OCP_DBL& a)
{
    x /= a;
    y /= a;
    z /= a;
    return *this;
}

Point3D operator*(const Point3D& p, const OCP_DBL& a)
{
    return Point3D(a * p.x, a * p.y, a * p.z);
}

Point3D operator*(const OCP_DBL& a, const Point3D& p)
{
    return Point3D(a * p.x, a * p.y, a * p.z);
}

Point3D CrossProduct(const Point3D& p1, const Point3D& p2)
{
    Point3D result;
    result.x = p1.y * p2.z - p1.z * p2.y;
    result.y = p1.z * p2.x - p1.x * p2.z;
    result.z = p1.x * p2.y - p1.y * p2.x;
    return result;
}

Point3D Matrix3::operator*(const Point3D& v) const
{
    Point3D result;
    result.x = M[0][0] * v.x + M[0][1] * v.y + M[0][2] * v.z;
    result.y = M[1][0] * v.x + M[1][1] * v.y + M[1][2] * v.z;
    result.z = M[2][0] * v.x + M[2][1] * v.y + M[2][2] * v.z;
    return result;
}


OCP_DBL Hexahedron::CalVolum() const
{
    OCP_DBL result =
        (p0.x * (p1.y * (-p2.z - p3.z + p4.z + p5.z) +
            p2.y * (p1.z - p3.z) +
            p3.y * (p1.z + p2.z - p4.z - p7.z) +
            p4.y * (-p1.z + p3.z - p5.z + p7.z) +
            p5.y * (-p1.z + p4.z) + p7.y * (p3.z - p4.z)) +
            p1.x * (p0.y * (+p2.z + p3.z - p4.z - p5.z) +
                p2.y * (-p0.z - p3.z + p5.z + p6.z) +
                p3.y * (-p0.z + p2.z) + p4.y * (p0.z - p5.z) +
                p5.y * (p0.z - p2.z + p4.z - p6.z) +
                p6.y * (-p2.z + p5.z)) +
            p2.x * (p0.y * (-p1.z + p3.z) +
                p1.y * (p0.z + p3.z - p5.z - p6.z) +
                p3.y * (-p0.z - p1.z + p6.z + p7.z) +
                p5.y * (p1.z - p6.z) +
                p6.y * (p1.z - p3.z + p5.z - p7.z) +
                p7.y * (-p3.z + p6.z)) +
            p3.x * (p0.y * (-p1.z - p2.z + p4.z + p7.z) +
                p1.y * (p0.z - p2.z) +
                p2.y * (p0.z + p1.z - p6.z - p7.z) +
                p4.y * (-p0.z + p7.z) + p6.y * (p2.z - p7.z) +
                p7.y * (-p0.z + p2.z - p4.z + p6.z)) +
            p4.x * (p0.y * (p1.z - p3.z + p5.z - p7.z) +
                p1.y * (-p0.z + p5.z) + p3.y * (p0.z - p7.z) +
                p5.y * (-p0.z - p1.z + p6.z + p7.z) +
                p6.y * (-p5.z + p7.z) +
                p7.y * (p0.z + p3.z - p5.z - p6.z)) +
            p5.x * (p0.y * (p1.z - p4.z) +
                p1.y * (-p0.z + p2.z - p4.z + p6.z) +
                p2.y * (-p1.z + p6.z) +
                p4.y * (p0.z + p1.z - p6.z - p7.z) +
                p6.y * (-p1.z - p2.z + p4.z + p7.z) +
                p7.y * (p4.z - p6.z)) +
            p6.x * (p1.y * (p2.z - p5.z) +
                p2.y * (-p1.z + p3.z - p5.z + p7.z) +
                p3.y * (-p2.z + p7.z) + p4.y * (p5.z - p7.z) +
                p5.y * (p1.z + p2.z - p4.z - p7.z) +
                p7.y * (-p2.z - p3.z + p4.z + p5.z)) +
            p7.x * (p0.y * (-p3.z + p4.z) + p2.y * (p3.z - p6.z) +
                p3.y * (p0.z - p2.z + p4.z - p6.z) +
                p4.y * (-p0.z - p3.z + p5.z + p6.z) +
                p5.y * (-p4.z + p6.z) +
                p6.y * (p2.z + p3.z - p4.z - p5.z))) /
        12;
    return result;
}

Point3D Hexahedron::CalCenter() const
{
    OCP_DBL r = 1.0 / 8.0;
    Point3D result = r * (p0 + p1 + p2 + p3 + p4 + p5 + p6 + p7);
    return result;
}


Point3D HexahedronFace::CalAreaVector() const
{
    Point3D v0, v1, v2, v3;
    v0 = p2 - p1;
    v1 = p0 - p1;
    v2 = p0 - p3;
    v3 = p2 - p3;
    Point3D result = 0.5 * (CrossProduct(v0, v1) + CrossProduct(v2, v3));
    return result;
}

Point3D HexahedronFace::CalCenter() const
{
    OCP_DBL r = 1.0 / 4.0;
    Point3D result = r * (p0 + p1 + p2 + p3);
    return result;
}

Point2D CalCrossingPoint(const Point2D Line1[2], const Point2D Line2[2])
{
    Point2D crosspoint;
    //
    //   LOCALS
    //
    OCP_DBL a11, a12, a21, a22, b1, b2, detA, detX, detY;
    //
    //   assume   x   =   crosspoint.x
    //            y   =   crosspoint.y
    //   calculate x and y with equations in the following
    //
    //    a11 a12     x       b1
    //   [        ] (   ) = (    )
    //    a21 a22     y       b2
    //
    a11 = Line1[1].y - Line1[0].y;
    a12 = Line1[0].x - Line1[1].x;
    a21 = Line2[1].y - Line2[0].y;
    a22 = Line2[0].x - Line2[1].x;
    b1 = a11 * Line1[0].x + a12 * Line1[0].y;
    b2 = a21 * Line2[0].x + a22 * Line2[0].y;
    detA = a11 * a22 - a12 * a21;

    if (fabs(detA) > 1E-10) {
        detX = b1 * a22 - b2 * a12;
        detY = a11 * b2 - a21 * b1;
        crosspoint.x = detX / detA;
        crosspoint.y = detY / detA;
    }
    else {
        crosspoint = Line1[0];
    }
    return crosspoint;
}

OCP_DBL CalAreaNotQuadr(const HexahedronFace& FACE1, const HexahedronFace& FACE2)
{
    // Attention! Only for non quadrilateral!!!  ---- Lishizhe
    //
    // This function calculate the common area of two quadrilaterals FACE1, FACE2.
    //
    // Order of points of Face follows
    //       1 --- 0        0 --- 1
    //       |     |    or  |     |
    //       2 --- 3        3 --- 2
    // p0, p1 are upper, p2, p3 are lower
    // y must be depth!!!
    //
    OCP_DBL CalAreaNotQuadr;
    //
    //   LOCALS
    //
    USI            iret;
    Point2D        crosspoint[4];
    Point2D        Line1[2], Line2[2];
    HexahedronFace FACEtmp1, FACEtmp2;
    Point3D        area, point1, point2, point3;
    //
    CalAreaNotQuadr = 0;
    iret = 0;
    //
    //   the crossing relations of 4 lines:
    //           Line1 : point0 and point1 of face1
    //           Line2 : point2 and point3 of face1
    //           Line3 : point0 and point1 of face2
    //           Line4 : point2 and point3 of face2
    //
    //   Line1 & Line3
    //
    Line1[0] = Point2D(FACE1.p0.x, FACE1.p0.y);
    Line1[1] = Point2D(FACE1.p1.x, FACE1.p1.y);
    Line2[0] = Point2D(FACE2.p0.x, FACE2.p0.y);
    Line2[1] = Point2D(FACE2.p1.x, FACE2.p1.y);
    crosspoint[0] = CalCrossingPoint(Line1, Line2);
    if ((crosspoint[0].x - Line1[0].x) * (crosspoint[0].x - Line1[1].x) < 0)
        iret = iret + 1;
    //
    //   Line2 & Line3
    //
    Line1[0] = Point2D(FACE1.p2.x, FACE1.p2.y);
    Line1[1] = Point2D(FACE1.p3.x, FACE1.p3.y);
    Line2[0] = Point2D(FACE2.p0.x, FACE2.p0.y);
    Line2[1] = Point2D(FACE2.p1.x, FACE2.p1.y);
    crosspoint[1] = CalCrossingPoint(Line1, Line2);
    if ((crosspoint[1].x - Line1[0].x) * (crosspoint[1].x - Line1[1].x) < 0)
        iret = iret + 2;
    //
    //   Line1 & Line4
    //
    Line1[0] = Point2D(FACE1.p0.x, FACE1.p0.y);
    Line1[1] = Point2D(FACE1.p1.x, FACE1.p1.y);
    Line2[0] = Point2D(FACE2.p2.x, FACE2.p2.y);
    Line2[1] = Point2D(FACE2.p3.x, FACE2.p3.y);
    crosspoint[2] = CalCrossingPoint(Line1, Line2);
    if ((crosspoint[2].x - Line1[0].x) * (crosspoint[2].x - Line1[1].x) < 0)
        iret = iret + 4;
    //
    //   Line2 & Line4
    //
    Line1[0] = Point2D(FACE1.p2.x, FACE1.p2.y);
    Line1[1] = Point2D(FACE1.p3.x, FACE1.p3.y);
    Line2[0] = Point2D(FACE2.p2.x, FACE2.p2.y);
    Line2[1] = Point2D(FACE2.p3.x, FACE2.p3.y);
    crosspoint[3] = CalCrossingPoint(Line1, Line2);
    if ((crosspoint[3].x - Line1[0].x) * (crosspoint[3].x - Line1[1].x) < 0)
        iret = iret + 8;
    //
    //   consider 12 cases of crossing relation combinations
    //
    switch (iret) {
    case 1:
        //
        //  Line1 & Line3 only
        //
        FACEtmp1.p1 = Point3D(crosspoint[0].x, crosspoint[0].y, 0);
        FACEtmp2.p0 = Point3D(crosspoint[0].x, crosspoint[0].y, 0);

        if (FACE1.p0.y > FACE2.p0.y) {
            FACEtmp1.p0 = FACE1.p0;
        }
        else {
            FACEtmp1.p0 = FACE2.p0;
        }

        if (FACE1.p1.y > FACE2.p1.y) {
            FACEtmp2.p1 = FACE1.p1;
        }
        else {
            FACEtmp2.p1 = FACE2.p1;
        }

        if (FACE1.p3.y > FACE2.p3.y) {
            FACEtmp1.p3 = FACE2.p3;
        }
        else {
            FACEtmp1.p3 = FACE1.p3;
        }

        if (FACE1.p2.y > FACE2.p2.y) {
            FACEtmp2.p2 = FACE2.p2;
        }
        else {
            FACEtmp2.p2 = FACE1.p2;
        }

        FACEtmp1.p2 = Point3D(0.5 * (FACEtmp1.p3.x + FACEtmp2.p2.x),
            0.5 * (FACEtmp1.p3.y + FACEtmp2.p2.y), 0);
        FACEtmp2.p3 = FACEtmp1.p2;

        area = FACEtmp1.CalAreaVector();
        CalAreaNotQuadr = fabs(area.z);
        area = FACEtmp2.CalAreaVector();
        CalAreaNotQuadr = CalAreaNotQuadr + fabs(area.z);
        break;
    case 2:
        //
        //  Line2 & Line3 only
        //
        if (FACE1.p3.y > FACE2.p0.y) {
            point1 = FACE1.p3;
            point2 = FACE2.p0;
        }
        else {
            point1 = FACE1.p2;
            point2 = FACE2.p1;
        }
        point3 = Point3D(crosspoint[1].x, crosspoint[1].y, 0);
        area = CrossProduct(point1 - point3, point2 - point3);
        CalAreaNotQuadr = fabs(area.z) * 0.5;
        break;
    case 3:
        //
        //  Line1 & Line3
        //  Line2 & Line3
        //
        FACEtmp1.p0 = Point3D(crosspoint[0].x, crosspoint[0].y, 0);
        FACEtmp1.p1 = Point3D(crosspoint[1].x, crosspoint[1].y, 0);
        if (FACE1.p0.y < FACE2.p0.y) {
            FACEtmp1.p2 = FACE1.p2;
            FACEtmp1.p3 = FACE1.p1;
        }
        else {
            FACEtmp1.p2 = FACE1.p3;
            FACEtmp1.p3 = FACE1.p0;
        }
        area = FACEtmp1.CalAreaVector();
        CalAreaNotQuadr = fabs(area.z);
        break;
    case 4:
        //
        //  Line1 & Line4 only
        //
        if (FACE1.p0.y < FACE2.p3.y) {
            point1 = FACE1.p0;
            point2 = FACE2.p3;
        }
        else {
            point1 = FACE1.p1;
            point2 = FACE2.p2;
        }
        point3 = Point3D(crosspoint[2].x, crosspoint[2].y, 0);
        area = CrossProduct(point1 - point3, point2 - point3);
        CalAreaNotQuadr = fabs(area.z) * 0.5;
        break;
    case 5:
        //
        //  Line1 & Line3
        //  Line1 & Line4
        //
        FACEtmp1.p0 = Point3D(crosspoint[0].x, crosspoint[0].y, 0);
        FACEtmp1.p1 = Point3D(crosspoint[2].x, crosspoint[2].y, 0);
        if (FACE2.p3.y > FACE1.p0.y) {
            FACEtmp1.p2 = FACE2.p3;
            FACEtmp1.p3 = FACE2.p0;
        }
        else {
            FACEtmp1.p2 = FACE2.p2;
            FACEtmp1.p3 = FACE2.p1;
        }
        area = FACEtmp1.CalAreaVector();
        CalAreaNotQuadr = fabs(area.z);
        break;
    case 8:
        //
        //  Line2 & Line4 only
        //
        FACEtmp1.p2 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);
        FACEtmp2.p3 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);

        if (FACE1.p0.y > FACE2.p0.y) {
            FACEtmp1.p0 = FACE1.p0;
        }
        else {
            FACEtmp1.p0 = FACE2.p0;
        }

        if (FACE1.p1.y > FACE2.p1.y) {
            FACEtmp2.p1 = FACE1.p1;
        }
        else {
            FACEtmp2.p1 = FACE2.p1;
        }

        if (FACE1.p3.y > FACE2.p3.y) {
            FACEtmp1.p3 = FACE2.p3;
        }
        else {
            FACEtmp1.p3 = FACE1.p3;
        }

        if (FACE1.p2.y > FACE2.p2.y) {
            FACEtmp2.p2 = FACE2.p2;
        }
        else {
            FACEtmp2.p2 = FACE1.p2;
        }

        FACEtmp1.p1 = Point3D(0.5 * (FACEtmp1.p0.x + FACEtmp2.p1.x),
            0.5 * (FACEtmp1.p0.y + FACEtmp2.p1.y), 0);
        FACEtmp2.p0 = FACEtmp1.p1;

        area = FACEtmp1.CalAreaVector();
        CalAreaNotQuadr = fabs(area.z);
        area = FACEtmp2.CalAreaVector();
        CalAreaNotQuadr = CalAreaNotQuadr + fabs(area.z);
        break;
    case 9:
        //
        //  Line1 & Line3
        //  Line2 & Line4
        //
        FACEtmp1.p1 = Point3D(crosspoint[0].x, crosspoint[0].y, 0);
        FACEtmp2.p0 = Point3D(crosspoint[0].x, crosspoint[0].y, 0);
        FACEtmp1.p2 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);
        FACEtmp2.p3 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);

        if (FACE1.p0.y > FACE2.p0.y) {
            FACEtmp1.p0 = FACE1.p0;
        }
        else {
            FACEtmp1.p0 = FACE2.p0;
        }

        if (FACE1.p1.y > FACE2.p1.y) {
            FACEtmp2.p1 = FACE1.p1;
        }
        else {
            FACEtmp2.p1 = FACE2.p1;
        }

        if (FACE1.p3.y > FACE2.p3.y) {
            FACEtmp1.p3 = FACE2.p3;
        }
        else {
            FACEtmp1.p3 = FACE1.p3;
        }

        if (FACE1.p2.y > FACE2.p2.y) {
            FACEtmp2.p2 = FACE2.p2;
        }
        else {
            FACEtmp2.p2 = FACE1.p2;
        }

        area = FACEtmp1.CalAreaVector();
        CalAreaNotQuadr = fabs(area.z);
        area = FACEtmp2.CalAreaVector();
        CalAreaNotQuadr = CalAreaNotQuadr + fabs(area.z);
        break;
    case 10:
        //
        //  Line2 & Line3
        //  Line2 & Line4
        //
        FACEtmp1.p0 = Point3D(crosspoint[1].x, crosspoint[1].y, 0);
        FACEtmp1.p1 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);
        if (FACE1.p2.y > FACE2.p1.y) {
            FACEtmp1.p2 = FACE2.p1;
            FACEtmp1.p3 = FACE2.p2;
        }
        else {
            FACEtmp1.p2 = FACE2.p3;
            FACEtmp1.p3 = FACE2.p0;
        }
        area = FACEtmp1.CalAreaVector();
        CalAreaNotQuadr = fabs(area.z);
        break;
    case 11:
        //
        //  Line1 & Line3
        //  Line2 & Line3
        //  Line2 & Line4
        //
        FACEtmp1.p0 = Point3D(crosspoint[0].x, crosspoint[0].y, 0);
        FACEtmp1.p2 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);
        FACEtmp1.p3 = Point3D(crosspoint[1].x, crosspoint[1].y, 0);
        FACEtmp2.p3 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);

        if (FACE1.p0.y < FACE2.p0.y) {
            FACEtmp2.p1 = FACE1.p1;
            FACEtmp2.p2 = FACE2.p2;
        }
        else {
            FACEtmp2.p1 = FACE1.p0;
            FACEtmp2.p2 = FACE2.p3;
        }

        FACEtmp1.p1 = Point3D(0.5 * (FACEtmp1.p0.x + FACEtmp2.p1.x),
            0.5 * (FACEtmp1.p0.y + FACEtmp2.p1.y), 0);
        FACEtmp2.p0 = FACEtmp1.p1;

        area = FACEtmp1.CalAreaVector();
        CalAreaNotQuadr = fabs(area.z);
        area = FACEtmp2.CalAreaVector();
        CalAreaNotQuadr = CalAreaNotQuadr + fabs(area.z);
        break;
    case 12:
        //
        //  Line1 & Line4
        //  Line2 & Line4
        //
        FACEtmp1.p0 = Point3D(crosspoint[2].x, crosspoint[2].y, 0);
        FACEtmp1.p1 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);
        if (FACE1.p2.y > FACE2.p2.y) {
            FACEtmp1.p2 = FACE1.p3;
            FACEtmp1.p3 = FACE1.p0;
        }
        else {
            FACEtmp1.p2 = FACE1.p2;
            FACEtmp1.p3 = FACE1.p1;
        }
        area = FACEtmp1.CalAreaVector();
        CalAreaNotQuadr = fabs(area.z);
        break;
    case 13:
        //
        //  Line1 & Line3
        //  Line1 & Line4
        //  Line2 & Line4
        //
        FACEtmp1.p2 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);
        FACEtmp2.p1 = Point3D(crosspoint[0].x, crosspoint[0].y, 0);
        FACEtmp2.p2 = Point3D(crosspoint[2].x, crosspoint[2].y, 0);
        FACEtmp2.p3 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);

        if (FACE1.p2.y > FACE2.p2.y) {
            FACEtmp1.p0 = FACE2.p0;
            FACEtmp1.p3 = FACE1.p3;
        }
        else {
            FACEtmp1.p0 = FACE2.p1;
            FACEtmp1.p3 = FACE1.p2;
        }

        FACEtmp1.p1 = Point3D(0.5 * (FACEtmp1.p0.x + FACEtmp2.p1.x),
            0.5 * (FACEtmp1.p0.y + FACEtmp2.p1.y), 0);
        FACEtmp2.p0 = FACEtmp1.p1;

        area = FACEtmp1.CalAreaVector();
        CalAreaNotQuadr = fabs(area.z);
        area = FACEtmp2.CalAreaVector();
        CalAreaNotQuadr = CalAreaNotQuadr + fabs(area.z);
        break;
    case 15:
        //
        //  Line1 & Line3
        //  Line2 & Line3
        //  Line1 & Line4
        //  Line2 & Line4
        //
        FACEtmp1.p0 = Point3D(crosspoint[0].x, crosspoint[0].y, 0);
        FACEtmp1.p1 = Point3D(crosspoint[2].x, crosspoint[2].y, 0);
        FACEtmp1.p2 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);
        FACEtmp1.p3 = Point3D(crosspoint[1].x, crosspoint[1].y, 0);
        area = FACEtmp1.CalAreaVector();
        CalAreaNotQuadr = fabs(area.z);
        break;
    default:
        CalAreaNotQuadr = 0;
        break;
    }
    return CalAreaNotQuadr;
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/