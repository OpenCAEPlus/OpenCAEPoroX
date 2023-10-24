/*! \file    CornerGrid.cpp
 *  \brief   Declaration of classes related to the corner grid
 *  \author  Shizhe Li
 *  \date    Nov/19/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "CornerGrid.hpp"


void ConnGrid::Allocate(const USI& max_neighbor)
{
    nConn   = 0;
    maxConn = max_neighbor;
    halfConn.resize(maxConn);
}

void ConnGrid::AddHalfConn(const OCP_ULL& n,
                           const Point3D& area,
                           const Point3D& d,
                           const ConnDirect&     direction,
                           const OCP_DBL& flag)
{
    if (nConn >= maxConn) {
        maxConn *= 2;
        halfConn.resize(maxConn);
        // get larger space
        if (maxConn > MAX_NEIGHBOR) {
            // OCP_ABORT("Too many Neighbors!");
        }
    }

    const OCP_DBL areaE = area * d / (d * d) * flag;

    if (!isfinite(areaE)) {
        if (d * d < TINY) {
            OCP_WARNING("Hexahedron reduces to a face!");
        }
        else {
            OCP_ABORT("Effective Area is NAN");
        }
    }

    halfConn[nConn].Ad_dd         = areaE;
    halfConn[nConn].d             = d;
    halfConn[nConn].neigh         = n;
    halfConn[nConn].directionType = direction;
    nConn++;
}

void OCP_COORD::Allocate(const USI& Nx, const USI& Ny, const USI& Nz)
{
    nx      = Nx;
    ny      = Ny;
    nz      = Nz;
    numGrid = nx * ny * nz;

    COORDDATA = new OCP_DBL**[3];
    for (USI i = 0; i < 3; i++) {
        COORDDATA[i] = new OCP_DBL*[2];
        for (USI j = 0; j < 2; j++) {
            COORDDATA[i][j] = new OCP_DBL[(nx + 1) * (ny + 1)];
        }
    }

    ZCORNDATA = new OCP_DBL***[nx];
    for (USI i = 0; i < nx; i++) {
        ZCORNDATA[i] = new OCP_DBL**[ny];
        for (USI j = 0; j < ny; j++) {
            ZCORNDATA[i][j] = new OCP_DBL*[nz];
            for (USI k = 0; k < nz; k++) {
                ZCORNDATA[i][j][k] = new OCP_DBL[8];
            }
        }
    }

    cornerPoints = new Hexahedron**[nx];
    for (USI i = 0; i < nx; i++) {
        cornerPoints[i] = new Hexahedron*[ny];
        for (USI j = 0; j < ny; j++) {
            cornerPoints[i][j] = new Hexahedron[nz];
        }
    }

    v.resize(numGrid);
    depth.resize(numGrid);
    dx.resize(numGrid);
    dy.resize(numGrid);
    dz.resize(numGrid);
    center.resize(numGrid);
}

void OCP_COORD::InputData(const vector<OCP_DBL>& coord, const vector<OCP_DBL>& zcorn)
{
    if (coord.empty() || !InputCOORDDATA(coord)) {
        OCP_ABORT("ERROR COORD!");
    }
    if (zcorn.empty() || !InputZCORNDATA(zcorn)) {
        OCP_ABORT("ERROR ZCORN!");
    }
}

OCP_BOOL OCP_COORD::InputCOORDDATA(const vector<OCP_DBL>& coord)
{
    // See Eclipse -- COORD
    OCP_BOOL flag = OCP_FALSE;
    OCP_ULL  iter = 0;

    for (USI J = 0; J < ny + 1; J++) {
        for (USI I = 0; I < nx + 1; I++) {
            // top
            for (USI i = 0; i < 3; i++) {
                COORDDATA[i][0][J * (nx + 1) + I] = coord[iter];
                iter++;
            }
            // bottom
            for (USI i = 0; i < 3; i++) {
                COORDDATA[i][1][J * (nx + 1) + I] = coord[iter];
                iter++;
            }
        }
    }

    flag = OCP_TRUE;
    return flag;
}

OCP_BOOL OCP_COORD::InputZCORNDATA(const vector<OCP_DBL>& zcorn)
{
    // See Eclipse -- ZCORN
    OCP_BOOL flag = OCP_FALSE;
    OCP_ULL  iter = 0;

    for (USI K = 0; K < nz; K++) {
        for (USI J = 0; J < ny; J++) {
            for (USI I = 0; I < nx; I++) {
                ZCORNDATA[I][J][K][0] = zcorn[iter];
                iter++;
                ZCORNDATA[I][J][K][1] = zcorn[iter];
                iter++;
            }
            for (USI I = 0; I < nx; I++) {
                ZCORNDATA[I][J][K][3] = zcorn[iter];
                iter++;
                ZCORNDATA[I][J][K][2] = zcorn[iter];
                iter++;
            }
        }
        for (USI J = 0; J < ny; J++) {
            for (USI I = 0; I < nx; I++) {
                ZCORNDATA[I][J][K][4] = zcorn[iter];
                iter++;
                ZCORNDATA[I][J][K][5] = zcorn[iter];
                iter++;
            }
            for (USI I = 0; I < nx; I++) {
                ZCORNDATA[I][J][K][7] = zcorn[iter];
                iter++;
                ZCORNDATA[I][J][K][6] = zcorn[iter];
                iter++;
            }
        }
    }

    flag = OCP_TRUE;
    return flag;
}

void OCP_COORD::SetAllFlags(const HexahedronFace& oFace, const HexahedronFace& Face)
{
    upNNC   = OCP_FALSE;
    downNNC = OCP_FALSE;

    // if the face reduce to a line
    if (sqrt((Face.p0 - Face.p1) * (Face.p0 - Face.p1)) <= TEENY && 
        sqrt((Face.p3 - Face.p2) * (Face.p3 - Face.p2)) <= TEENY) {
        interFace = Face;
        flagJump  = OCP_TRUE;     
    }
    else {
        // if the i th point of oFace is deeper than the one of Face, then flagpi = 1;
        // if the i th point of oFace is higher than the one of Face, then flagpi = -1;
        // if the i th point of oFace is very close to the one of Face, then flagpi = 0;
        OCP_INT flagp0, flagp1, flagp2, flagp3;

        interFace = Face;
        if (oFace.p0.z > Face.p0.z + TEENY) {
            interFace.p0 = oFace.p0;
            flagp0       = 1;
            upNNC        = OCP_TRUE;
        }
        else if (oFace.p0.z < Face.p0.z - TEENY)  flagp0 = -1;
        else                                      flagp0 = 0;
    
        if (oFace.p1.z < Face.p1.z - TEENY) {
            interFace.p1 = oFace.p1;
            flagp1       = -1;
            downNNC      = OCP_TRUE;
        }
        else if (oFace.p1.z > Face.p1.z + TEENY)  flagp1 = 1;
        else                                      flagp1 = 0;
       
        if (oFace.p2.z < Face.p2.z - TEENY) {
            interFace.p2 = oFace.p2;
            flagp2       = -1;
            downNNC      = OCP_TRUE;
        }
        else if (oFace.p2.z > Face.p2.z + TEENY)  flagp2 = 1;
        else                                      flagp2 = 0;

        if (oFace.p3.z > Face.p3.z + TEENY) {
            interFace.p3 = oFace.p3;
            flagp3       = 1;
            upNNC        = OCP_TRUE;
        }
        else if (oFace.p3.z < Face.p3.z - TEENY)  flagp3 = -1;
        else                                      flagp3 = 0;

        // check if interface is empty set
        // check if interface is quadrilateral
        // check if the one contains the other one

        if (((oFace.p1.z <= Face.p0.z) && (oFace.p2.z <= Face.p3.z)) ||
            ((oFace.p0.z >= Face.p1.z) && (oFace.p3.z >= Face.p2.z))) {
            flagJump = OCP_TRUE;
        }
        else {
            flagJump = OCP_FALSE;
            if ((flagp0 * flagp3 >= 0) && (oFace.p0.z <= Face.p1.z) &&
                (oFace.p3.z <= Face.p2.z) && (flagp1 * flagp2 >= 0) &&
                (oFace.p1.z >= Face.p0.z) && (oFace.p2.z >= Face.p3.z)) {
                flagQuad = OCP_TRUE;
            }
            else {
                flagQuad = OCP_FALSE;
            }
        }
    }
}

void OCP_COORD::SetupCornerPoints()
{
    OCP_ULL cindex, oindex; // current block index and the other block index
    OCP_ULL nxny = nx * ny;

    // allocate memoery for connections
    vector<ConnGrid> blockconn(numGrid);
    for (OCP_ULL iloop = 0; iloop < numGrid; iloop++) {
        blockconn[iloop].Allocate(10);
    }

    // setup each block including coordinates of points, center, depth, and volume
    OCP_DBL xtop, ytop, ztop, xbottom, ybottom, zbottom, xvalue, yvalue, zvalue;
    for (USI k = 0; k < nz; k++) {
        for (USI j = 0; j < ny; j++) {
            for (USI i = 0; i < nx; i++) {
                //
                // corner point 0 and 4
                //
                xtop    = COORDDATA[0][0][j * (nx + 1) + i];
                ytop    = COORDDATA[1][0][j * (nx + 1) + i];
                ztop    = COORDDATA[2][0][j * (nx + 1) + i];
                xbottom = COORDDATA[0][1][j * (nx + 1) + i];
                ybottom = COORDDATA[1][1][j * (nx + 1) + i];
                zbottom = COORDDATA[2][1][j * (nx + 1) + i];

                if (fabs(xtop - xbottom) < TINY && fabs(ytop - ybottom) < TINY) {
                    xvalue = xtop;
                    yvalue = ytop;
                    cornerPoints[i][j][k].p0 = Point3D(xvalue, yvalue, ZCORNDATA[i][j][k][0]);
                    cornerPoints[i][j][k].p4 = Point3D(xvalue, yvalue, ZCORNDATA[i][j][k][4]);
                }
                else {
                    zvalue = ZCORNDATA[i][j][k][0];
                    xvalue = xbottom - (zbottom - zvalue) / (zbottom - ztop) * (xbottom - xtop);
                    yvalue = ybottom - (zbottom - zvalue) / (zbottom - ztop) * (ybottom - ytop);
                    cornerPoints[i][j][k].p0 = Point3D(xvalue, yvalue, zvalue);

                    zvalue = ZCORNDATA[i][j][k][4];
                    xvalue = xbottom - (zbottom - zvalue) / (zbottom - ztop) * (xbottom - xtop);
                    yvalue = ybottom - (zbottom - zvalue) / (zbottom - ztop) * (ybottom - ytop);
                    cornerPoints[i][j][k].p4 = Point3D(xvalue, yvalue, zvalue);
                }

                //
                //    corner point 1 and 5
                //

                xtop    = COORDDATA[0][0][j * (nx + 1) + i + 1];
                ytop    = COORDDATA[1][0][j * (nx + 1) + i + 1];
                ztop    = COORDDATA[2][0][j * (nx + 1) + i + 1];
                xbottom = COORDDATA[0][1][j * (nx + 1) + i + 1];
                ybottom = COORDDATA[1][1][j * (nx + 1) + i + 1];
                zbottom = COORDDATA[2][1][j * (nx + 1) + i + 1];

                if (fabs(xtop - xbottom) < TINY && fabs(ytop - ybottom) < TINY) {
                    xvalue = xtop;
                    yvalue = ytop;
                    cornerPoints[i][j][k].p1 = Point3D(xvalue, yvalue, ZCORNDATA[i][j][k][1]);
                    cornerPoints[i][j][k].p5 = Point3D(xvalue, yvalue, ZCORNDATA[i][j][k][5]);
                }
                else {
                    zvalue = ZCORNDATA[i][j][k][1];
                    xvalue = xbottom - (zbottom - zvalue) / (zbottom - ztop) * (xbottom - xtop);
                    yvalue = ybottom - (zbottom - zvalue) / (zbottom - ztop) * (ybottom - ytop);
                    cornerPoints[i][j][k].p1 = Point3D(xvalue, yvalue, zvalue);

                    zvalue = ZCORNDATA[i][j][k][5];
                    xvalue = xbottom - (zbottom - zvalue) / (zbottom - ztop) * (xbottom - xtop);
                    yvalue = ybottom - (zbottom - zvalue) / (zbottom - ztop) * (ybottom - ytop);
                    cornerPoints[i][j][k].p5 = Point3D(xvalue, yvalue, zvalue);
                }


                //
                //    corner point 2 and 6
                //
                xtop    = COORDDATA[0][0][(j + 1) * (nx + 1) + i + 1];
                ytop    = COORDDATA[1][0][(j + 1) * (nx + 1) + i + 1];
                ztop    = COORDDATA[2][0][(j + 1) * (nx + 1) + i + 1];
                xbottom = COORDDATA[0][1][(j + 1) * (nx + 1) + i + 1];
                ybottom = COORDDATA[1][1][(j + 1) * (nx + 1) + i + 1];
                zbottom = COORDDATA[2][1][(j + 1) * (nx + 1) + i + 1];

                if (fabs(xtop - xbottom) < TINY && fabs(ytop - ybottom) < TINY) {
                    xvalue = xtop;
                    yvalue = ytop;
                    cornerPoints[i][j][k].p2 = Point3D(xvalue, yvalue, ZCORNDATA[i][j][k][2]);
                    cornerPoints[i][j][k].p6 = Point3D(xvalue, yvalue, ZCORNDATA[i][j][k][6]);
                }
                else {
                    zvalue = ZCORNDATA[i][j][k][2];
                    xvalue = xbottom - (zbottom - zvalue) / (zbottom - ztop) * (xbottom - xtop);
                    yvalue = ybottom - (zbottom - zvalue) / (zbottom - ztop) * (ybottom - ytop);
                    cornerPoints[i][j][k].p2 = Point3D(xvalue, yvalue, zvalue);

                    zvalue = ZCORNDATA[i][j][k][6];
                    xvalue = xbottom - (zbottom - zvalue) / (zbottom - ztop) * (xbottom - xtop);
                    yvalue = ybottom - (zbottom - zvalue) / (zbottom - ztop) * (ybottom - ytop);
                    cornerPoints[i][j][k].p6 = Point3D(xvalue, yvalue, zvalue);
                }

                //
                //    corner point 3 and 7
                //
                xtop    = COORDDATA[0][0][(j + 1) * (nx + 1) + i];
                ytop    = COORDDATA[1][0][(j + 1) * (nx + 1) + i];
                ztop    = COORDDATA[2][0][(j + 1) * (nx + 1) + i];
                xbottom = COORDDATA[0][1][(j + 1) * (nx + 1) + i];
                ybottom = COORDDATA[1][1][(j + 1) * (nx + 1) + i];
                zbottom = COORDDATA[2][1][(j + 1) * (nx + 1) + i];

                if (fabs(xtop - xbottom) < TINY && fabs(ytop - ybottom) < TINY) {
                    xvalue = xtop;
                    yvalue = ytop;
                    cornerPoints[i][j][k].p3 = Point3D(xvalue, yvalue, ZCORNDATA[i][j][k][3]);
                    cornerPoints[i][j][k].p7 = Point3D(xvalue, yvalue, ZCORNDATA[i][j][k][7]);
                }
                else {
                    zvalue = ZCORNDATA[i][j][k][3];
                    xvalue = xbottom - (zbottom - zvalue) / (zbottom - ztop) * (xbottom - xtop);
                    yvalue = ybottom - (zbottom - zvalue) / (zbottom - ztop) * (ybottom - ytop);
                    cornerPoints[i][j][k].p3 = Point3D(xvalue, yvalue, zvalue);

                    zvalue = ZCORNDATA[i][j][k][7];
                    xvalue = xbottom - (zbottom - zvalue) / (zbottom - ztop) * (xbottom - xtop);
                    yvalue = ybottom - (zbottom - zvalue) / (zbottom - ztop) * (ybottom - ytop);
                    cornerPoints[i][j][k].p7 = Point3D(xvalue, yvalue, zvalue);
                }

                //    calculate volumes and pore volumes
                cindex = k * nxny + j * nx + i;
                //
                // NOTE: if there are several points not well ordered, the calculated
                // volume will be negative.
                //
                v[cindex]      = fabs(cornerPoints[i][j][k].CalVolum());
                center[cindex] = cornerPoints[i][j][k].CalCenter();
                depth[cindex]  = center[cindex].z;
            }
        }
    }

    // find neighbor and calculate transmissibility
    OCP_ULL num_conn = 0;         // record the num of connection, a->b & b->a are both included
    Point3D Pcenter, Pface, Pc2f; // center of Hexahedron
    HexahedronFace Face, oFace;   // current face, the other face
    HexahedronFace FaceP, oFaceP; // Projection of Face and the other face
    Point3D        areaV;         // area vector of interface
    OCP_DBL        areaP;         // area of projection of interface
    OCP_INT        iznnc;
    Point3D        dxpoint, dypoint, dzpoint;

    /////////////////////////////////////////////////////////////////////
    // Attention that The coordinate axis follows the right-hand rule ! //
    /////////////////////////////////////////////////////////////////////
    //
    //      o----> x
    //     /|
    //    y z
    // For a face, p0 and p3 are the points of upper edge of quadrilateral,
    // p1 and p2 are the points of lower edge
    //       p0 ---- p3
    //        |       |
    //        |       |
    //       p1 ---- p2

    // Determine flagForward
    if (COORDDATA[1][0][nx + 1] > COORDDATA[1][0][0])
        flagForward = 1.0;
    else
        flagForward = -1.0;

    ConnDirect direction;

    for (USI k = 0; k < nz; k++) {
        for (USI j = 0; j < ny; j++) {
            for (USI i = 0; i < nx; i++) {
                // begin from each block
                const Hexahedron& block = cornerPoints[i][j][k];
                cindex                  = k * nxny + j * nx + i;
                Pcenter                 = center[cindex];

                // cout << "============= " << cindex << " =============" << endl;
                //
                // (x-) direction
                //

                direction = ConnDirect::xm;

                Face.p0 = block.p0;
                Face.p1 = block.p4;
                Face.p2 = block.p7;
                Face.p3 = block.p3;
                Pface   = Face.CalCenter();
                Pc2f    = Pface - Pcenter;
                dxpoint = Pc2f;

                if (i == 0) {
                    // nothing to do
                } else {

                    const Hexahedron& leftblock = cornerPoints[i - 1][j][k];
                    oindex                      = k * nxny + j * nx + i - 1;

                    oFace.p0 = leftblock.p1;
                    oFace.p1 = leftblock.p5;
                    oFace.p2 = leftblock.p6;
                    oFace.p3 = leftblock.p2;

                    SetAllFlags(oFace, Face);

                    // calculate the interface of two face
                    if (flagJump) {
                        // nothing to do
                    } else {
                        if (flagQuad) {
                            areaV = interFace.CalAreaVector();
                        } else {
                            FaceP.p0  = Point3D(Face.p3.y, Face.p3.z, 0);
                            FaceP.p1  = Point3D(Face.p0.y, Face.p0.z, 0);
                            FaceP.p2  = Point3D(Face.p1.y, Face.p1.z, 0);
                            FaceP.p3  = Point3D(Face.p2.y, Face.p2.z, 0);
                            oFaceP.p0 = Point3D(oFace.p3.y, oFace.p3.z, 0);
                            oFaceP.p1 = Point3D(oFace.p0.y, oFace.p0.z, 0);
                            oFaceP.p2 = Point3D(oFace.p1.y, oFace.p1.z, 0);
                            oFaceP.p3 = Point3D(oFace.p2.y, oFace.p2.z, 0);
                            areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                            // attention the direction of vector
                            areaV = Face.CalAreaVector();
                            // correct
                            if (fabs(areaV.x) < 1E-6) {
                                OCP_WARNING("x is too small");
                            } else {
                                areaV.y = areaV.y / fabs(areaV.x) * areaP;
                                areaV.z = areaV.z / fabs(areaV.x) * areaP;
                                areaV.x = OCP_SIGN(areaV.x) * areaP;
                            }
                        }
                        blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, direction,
                                                      flagForward);
                        num_conn++;
                    }

                    // then find all NNC for current block
                    // check if upNNC and downNNC exist

                    direction = ConnDirect::x;

                    iznnc = -1;
                    while (upNNC) {
                        // if (-iznnc > k) break;
                        if (-iznnc - static_cast<OCP_INT>(k) > 0) break;
                        // find object block
                        const Hexahedron& leftblock = cornerPoints[i - 1][j][k + iznnc];
                        oindex   = (k + iznnc) * nxny + j * nx + i - 1;
                        oFace.p0 = leftblock.p1;
                        oFace.p1 = leftblock.p5;
                        oFace.p2 = leftblock.p6;
                        oFace.p3 = leftblock.p2;

                        SetAllFlags(oFace, Face);

                        // calculate the interface of two face
                        if (flagJump) {
                            // nothing to do
                        } else {
                            if (flagQuad) {
                                areaV = interFace.CalAreaVector();
                            } else {
                                FaceP.p0  = Point3D(Face.p3.y, Face.p3.z, 0);
                                FaceP.p1  = Point3D(Face.p0.y, Face.p0.z, 0);
                                FaceP.p2  = Point3D(Face.p1.y, Face.p1.z, 0);
                                FaceP.p3  = Point3D(Face.p2.y, Face.p2.z, 0);
                                oFaceP.p0 = Point3D(oFace.p3.y, oFace.p3.z, 0);
                                oFaceP.p1 = Point3D(oFace.p0.y, oFace.p0.z, 0);
                                oFaceP.p2 = Point3D(oFace.p1.y, oFace.p1.z, 0);
                                oFaceP.p3 = Point3D(oFace.p2.y, oFace.p2.z, 0);
                                areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                                // attention the direction of vector
                                areaV = Face.CalAreaVector();
                                // correct
                                if (fabs(areaV.x) < 1E-6) {
                                    OCP_WARNING("x is too small");
                                } else {
                                    areaV.y = areaV.y / fabs(areaV.x) * areaP;
                                    areaV.z = areaV.z / fabs(areaV.x) * areaP;
                                    areaV.x = OCP_SIGN(areaV.x) * areaP;
                                }
                            }
                            blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, direction,
                                                          flagForward);
                            num_conn++;
                        }
                        iznnc--;
                    }

                    iznnc = 1;
                    while (downNNC) {
                        if (k + iznnc > nz - 1) break;
                        // find object block
                        const Hexahedron& leftblock = cornerPoints[i - 1][j][k + iznnc];
                        oindex   = (k + iznnc) * nxny + j * nx + i - 1;
                        oFace.p0 = leftblock.p1;
                        oFace.p1 = leftblock.p5;
                        oFace.p2 = leftblock.p6;
                        oFace.p3 = leftblock.p2;

                        SetAllFlags(oFace, Face);

                        // calculate the interface of two face
                        if (flagJump) {
                            // nothing to do
                        } else {
                            if (flagQuad) {
                                areaV = interFace.CalAreaVector();
                            } else {
                                FaceP.p0  = Point3D(Face.p3.y, Face.p3.z, 0);
                                FaceP.p1  = Point3D(Face.p0.y, Face.p0.z, 0);
                                FaceP.p2  = Point3D(Face.p1.y, Face.p1.z, 0);
                                FaceP.p3  = Point3D(Face.p2.y, Face.p2.z, 0);
                                oFaceP.p0 = Point3D(oFace.p3.y, oFace.p3.z, 0);
                                oFaceP.p1 = Point3D(oFace.p0.y, oFace.p0.z, 0);
                                oFaceP.p2 = Point3D(oFace.p1.y, oFace.p1.z, 0);
                                oFaceP.p3 = Point3D(oFace.p2.y, oFace.p2.z, 0);
                                areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                                // attention the direction of vector
                                areaV = Face.CalAreaVector();
                                // correct
                                if (fabs(areaV.x) < 1E-6) {
                                    OCP_WARNING("x is too small");
                                } else {
                                    areaV.y = areaV.y / fabs(areaV.x) * areaP;
                                    areaV.z = areaV.z / fabs(areaV.x) * areaP;
                                    areaV.x = OCP_SIGN(areaV.x) * areaP;
                                }
                            }
                            blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, direction,
                                                          flagForward);
                            num_conn++;
                        }
                        iznnc++;
                    }
                }

                //
                // (x+) direction
                //
                direction = ConnDirect::xp;

                Face.p0 = block.p2;
                Face.p1 = block.p6;
                Face.p2 = block.p5;
                Face.p3 = block.p1;
                Pface   = Face.CalCenter();
                Pc2f    = Pface - Pcenter;
                dxpoint = Pc2f - dxpoint;

                if (i == nx - 1) {
                    // nothing to do
                } else {

                    const Hexahedron& rightblock = cornerPoints[i + 1][j][k];
                    oindex                       = k * nxny + j * nx + i + 1;

                    oFace.p0 = rightblock.p3;
                    oFace.p1 = rightblock.p7;
                    oFace.p2 = rightblock.p4;
                    oFace.p3 = rightblock.p0;

                    SetAllFlags(oFace, Face);

                    // calculate the interface of two face
                    if (flagJump) {
                        // nothing to do
                    } else {
                        if (flagQuad) {
                            areaV = interFace.CalAreaVector();
                        } else {
                            FaceP.p0  = Point3D(Face.p3.y, Face.p3.z, 0);
                            FaceP.p1  = Point3D(Face.p0.y, Face.p0.z, 0);
                            FaceP.p2  = Point3D(Face.p1.y, Face.p1.z, 0);
                            FaceP.p3  = Point3D(Face.p2.y, Face.p2.z, 0);
                            oFaceP.p0 = Point3D(oFace.p3.y, oFace.p3.z, 0);
                            oFaceP.p1 = Point3D(oFace.p0.y, oFace.p0.z, 0);
                            oFaceP.p2 = Point3D(oFace.p1.y, oFace.p1.z, 0);
                            oFaceP.p3 = Point3D(oFace.p2.y, oFace.p2.z, 0);
                            areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                            // attention the direction of vector
                            areaV = Face.CalAreaVector();
                            // correct
                            if (fabs(areaV.x) < 1E-6) {
                                OCP_WARNING("x is too small");
                            } else {
                                areaV.y = areaV.y / fabs(areaV.x) * areaP;
                                areaV.z = areaV.z / fabs(areaV.x) * areaP;
                                areaV.x = OCP_SIGN(areaV.x) * areaP;
                            }
                        }
                        blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, direction,
                                                      flagForward);
                        num_conn++;
                    }

                    // then find all NNC for current block
                    direction = ConnDirect::x;

                    iznnc = -1;
                    while (upNNC) {
                        // if (-iznnc > k) break;
                        if (-iznnc - static_cast<OCP_INT>(k) > 0) break;
                        // find object block
                        const Hexahedron& rightblock =
                            cornerPoints[i + 1][j][k + iznnc];
                        oindex   = (k + iznnc) * nxny + j * nx + i + 1;
                        oFace.p0 = rightblock.p3;
                        oFace.p1 = rightblock.p7;
                        oFace.p2 = rightblock.p4;
                        oFace.p3 = rightblock.p0;

                        SetAllFlags(oFace, Face);

                        // calculate the interface of two face
                        if (flagJump) {
                            // nothing to do
                        } else {
                            if (flagQuad) {
                                areaV = interFace.CalAreaVector();
                            } else {
                                FaceP.p0  = Point3D(Face.p3.y, Face.p3.z, 0);
                                FaceP.p1  = Point3D(Face.p0.y, Face.p0.z, 0);
                                FaceP.p2  = Point3D(Face.p1.y, Face.p1.z, 0);
                                FaceP.p3  = Point3D(Face.p2.y, Face.p2.z, 0);
                                oFaceP.p0 = Point3D(oFace.p3.y, oFace.p3.z, 0);
                                oFaceP.p1 = Point3D(oFace.p0.y, oFace.p0.z, 0);
                                oFaceP.p2 = Point3D(oFace.p1.y, oFace.p1.z, 0);
                                oFaceP.p3 = Point3D(oFace.p2.y, oFace.p2.z, 0);
                                areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                                // attention the direction of vector
                                areaV = Face.CalAreaVector();
                                // correct
                                if (fabs(areaV.x) < 1E-6) {
                                    OCP_WARNING("x is too small");
                                } else {
                                    areaV.y = areaV.y / fabs(areaV.x) * areaP;
                                    areaV.z = areaV.z / fabs(areaV.x) * areaP;
                                    areaV.x = OCP_SIGN(areaV.x) * areaP;
                                }
                            }
                            blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, direction,
                                                          flagForward);
                            num_conn++;
                        }
                        iznnc--;
                    }

                    iznnc = 1;
                    while (downNNC) {
                        if (k + iznnc > nz - 1) break;
                        // find object block
                        const Hexahedron& rightblock =
                            cornerPoints[i + 1][j][k + iznnc];
                        oindex   = (k + iznnc) * nxny + j * nx + i + 1;
                        oFace.p0 = rightblock.p3;
                        oFace.p1 = rightblock.p7;
                        oFace.p2 = rightblock.p4;
                        oFace.p3 = rightblock.p0;

                        SetAllFlags(oFace, Face);

                        // calculate the interface of two face
                        if (flagJump) {
                            // nothing to do
                        } else {
                            if (flagQuad) {
                                areaV = interFace.CalAreaVector();
                            } else {
                                FaceP.p0  = Point3D(Face.p3.y, Face.p3.z, 0);
                                FaceP.p1  = Point3D(Face.p0.y, Face.p0.z, 0);
                                FaceP.p2  = Point3D(Face.p1.y, Face.p1.z, 0);
                                FaceP.p3  = Point3D(Face.p2.y, Face.p2.z, 0);
                                oFaceP.p0 = Point3D(oFace.p3.y, oFace.p3.z, 0);
                                oFaceP.p1 = Point3D(oFace.p0.y, oFace.p0.z, 0);
                                oFaceP.p2 = Point3D(oFace.p1.y, oFace.p1.z, 0);
                                oFaceP.p3 = Point3D(oFace.p2.y, oFace.p2.z, 0);
                                areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                                // attention the direction of vector
                                areaV = Face.CalAreaVector();
                                // correct
                                if (fabs(areaV.x) < 1E-6) {
                                    OCP_WARNING("x is too small");
                                } else {
                                    areaV.y = areaV.y / fabs(areaV.x) * areaP;
                                    areaV.z = areaV.z / fabs(areaV.x) * areaP;
                                    areaV.x = OCP_SIGN(areaV.x) * areaP;
                                }
                            }
                            blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, direction,
                                                          flagForward);
                            num_conn++;
                        }
                        iznnc++;
                    }
                }

                //
                // (y-) direction
                //
                direction = ConnDirect::ym;

                Face.p0 = block.p1;
                Face.p1 = block.p5;
                Face.p2 = block.p4;
                Face.p3 = block.p0;
                Pface   = Face.CalCenter();
                Pc2f    = Pface - Pcenter;
                dypoint = Pc2f;

                if (j == 0) {
                    // nothing to do
                } else {

                    const Hexahedron& backblock = cornerPoints[i][j - 1][k];
                    oindex                      = k * nxny + (j - 1) * nx + i;

                    oFace.p0 = backblock.p2;
                    oFace.p1 = backblock.p6;
                    oFace.p2 = backblock.p7;
                    oFace.p3 = backblock.p3;

                    SetAllFlags(oFace, Face);

                    // calculate the interface of two face
                    if (flagJump) {
                        // nothing to do
                    } else {
                        if (flagQuad) {
                            areaV = interFace.CalAreaVector();
                        } else {
                            FaceP.p0  = Point3D(Face.p0.x, Face.p0.z, 0);
                            FaceP.p1  = Point3D(Face.p3.x, Face.p3.z, 0);
                            FaceP.p2  = Point3D(Face.p2.x, Face.p2.z, 0);
                            FaceP.p3  = Point3D(Face.p1.x, Face.p1.z, 0);
                            oFaceP.p0 = Point3D(oFace.p0.x, oFace.p0.z, 0);
                            oFaceP.p1 = Point3D(oFace.p3.x, oFace.p3.z, 0);
                            oFaceP.p2 = Point3D(oFace.p2.x, oFace.p2.z, 0);
                            oFaceP.p3 = Point3D(oFace.p1.x, oFace.p1.z, 0);
                            areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                            // attention the direction of vector
                            areaV = Face.CalAreaVector();
                            // correct
                            if (fabs(areaV.y) < 1E-6) {
                                OCP_WARNING("y is too small");
                            } else {
                                areaV.x = areaV.x / fabs(areaV.y) * areaP;
                                areaV.z = areaV.z / fabs(areaV.y) * areaP;
                                areaV.y = OCP_SIGN(areaV.y) * areaP;
                            }
                        }
                        blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, direction,
                                                      flagForward);
                        num_conn++;
                    }

                    // then find all NNC for current block
                    direction = ConnDirect::y;

                    iznnc = -1;
                    while (upNNC) {
                        // if (-iznnc > k) break;
                        if (-iznnc - static_cast<OCP_INT>(k) > 0) break;
                        // find object block
                        const Hexahedron& backblock = cornerPoints[i][j - 1][k + iznnc];
                        oindex   = (k + iznnc) * nxny + (j - 1) * nx + i;
                        oFace.p0 = backblock.p2;
                        oFace.p1 = backblock.p6;
                        oFace.p2 = backblock.p7;
                        oFace.p3 = backblock.p3;

                        SetAllFlags(oFace, Face);

                        // calculate the interface of two face
                        if (flagJump) {
                            // nothing to do
                        } else {
                            if (flagQuad) {
                                areaV = interFace.CalAreaVector();
                            } else {
                                FaceP.p0  = Point3D(Face.p0.x, Face.p0.z, 0);
                                FaceP.p1  = Point3D(Face.p3.x, Face.p3.z, 0);
                                FaceP.p2  = Point3D(Face.p2.x, Face.p2.z, 0);
                                FaceP.p3  = Point3D(Face.p1.x, Face.p1.z, 0);
                                oFaceP.p0 = Point3D(oFace.p0.x, oFace.p0.z, 0);
                                oFaceP.p1 = Point3D(oFace.p3.x, oFace.p3.z, 0);
                                oFaceP.p2 = Point3D(oFace.p2.x, oFace.p2.z, 0);
                                oFaceP.p3 = Point3D(oFace.p1.x, oFace.p1.z, 0);
                                areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                                // attention the direction of vector
                                areaV = Face.CalAreaVector();
                                // correct
                                if (fabs(areaV.y) < 1E-6) {
                                    OCP_WARNING("y is too small");
                                } else {
                                    areaV.x = areaV.x / fabs(areaV.y) * areaP;
                                    areaV.z = areaV.z / fabs(areaV.y) * areaP;
                                    areaV.y = OCP_SIGN(areaV.y) * areaP;
                                }
                            }
                            blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, direction,
                                                          flagForward);
                            num_conn++;
                        }
                        iznnc--;
                    }

                    iznnc = 1;
                    while (downNNC) {
                        if (k + iznnc > nz - 1) break;
                        // find object block
                        const Hexahedron& backblock = cornerPoints[i][j - 1][k + iznnc];
                        oindex   = (k + iznnc) * nxny + (j - 1) * nx + i;
                        oFace.p0 = backblock.p2;
                        oFace.p1 = backblock.p6;
                        oFace.p2 = backblock.p7;
                        oFace.p3 = backblock.p3;

                        SetAllFlags(oFace, Face);

                        // calculate the interface of two face
                        if (flagJump) {
                            // nothing to do
                        } else {
                            if (flagQuad) {
                                areaV = interFace.CalAreaVector();
                            } else {
                                FaceP.p0  = Point3D(Face.p0.x, Face.p0.z, 0);
                                FaceP.p1  = Point3D(Face.p3.x, Face.p3.z, 0);
                                FaceP.p2  = Point3D(Face.p2.x, Face.p2.z, 0);
                                FaceP.p3  = Point3D(Face.p1.x, Face.p1.z, 0);
                                oFaceP.p0 = Point3D(oFace.p0.x, oFace.p0.z, 0);
                                oFaceP.p1 = Point3D(oFace.p3.x, oFace.p3.z, 0);
                                oFaceP.p2 = Point3D(oFace.p2.x, oFace.p2.z, 0);
                                oFaceP.p3 = Point3D(oFace.p1.x, oFace.p1.z, 0);
                                areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                                // attention the direction of vector
                                areaV = Face.CalAreaVector();
                                // correct
                                if (fabs(areaV.y) < 1E-6) {
                                    OCP_WARNING("y is too small");
                                } else {
                                    areaV.x = areaV.x / fabs(areaV.y) * areaP;
                                    areaV.z = areaV.z / fabs(areaV.y) * areaP;
                                    areaV.y = OCP_SIGN(areaV.y) * areaP;
                                }
                            }
                            blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, direction,
                                                          flagForward);
                            num_conn++;
                        }
                        iznnc++;
                    }
                }

                //
                // (y+) direction
                //
                direction = ConnDirect::yp;

                Face.p0 = block.p3;
                Face.p1 = block.p7;
                Face.p2 = block.p6;
                Face.p3 = block.p2;
                Pface   = Face.CalCenter();
                Pc2f    = Pface - Pcenter;
                dypoint = Pc2f - dypoint;

                if (j == ny - 1) {
                    // nothing to do
                } else {

                    const Hexahedron& frontblock = cornerPoints[i][j + 1][k];
                    oindex                       = k * nxny + (j + 1) * nx + i;

                    oFace.p0 = frontblock.p0;
                    oFace.p1 = frontblock.p4;
                    oFace.p2 = frontblock.p5;
                    oFace.p3 = frontblock.p1;

                    SetAllFlags(oFace, Face);

                    // calculate the interface of two face
                    if (flagJump) {
                        // nothing to do
                    } else {
                        if (flagQuad) {
                            areaV = interFace.CalAreaVector();
                        } else {
                            FaceP.p0  = Point3D(Face.p0.x, Face.p0.z, 0);
                            FaceP.p1  = Point3D(Face.p3.x, Face.p3.z, 0);
                            FaceP.p2  = Point3D(Face.p2.x, Face.p2.z, 0);
                            FaceP.p3  = Point3D(Face.p1.x, Face.p1.z, 0);
                            oFaceP.p0 = Point3D(oFace.p0.x, oFace.p0.z, 0);
                            oFaceP.p1 = Point3D(oFace.p3.x, oFace.p3.z, 0);
                            oFaceP.p2 = Point3D(oFace.p2.x, oFace.p2.z, 0);
                            oFaceP.p3 = Point3D(oFace.p1.x, oFace.p1.z, 0);
                            areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                            // attention the direction of vector
                            areaV = Face.CalAreaVector();
                            // correct
                            if (fabs(areaV.y) < 1E-6) {
                                OCP_WARNING("y is too small");
                            } else {
                                areaV.x = areaV.x / fabs(areaV.y) * areaP;
                                areaV.z = areaV.z / fabs(areaV.y) * areaP;
                                areaV.y = OCP_SIGN(areaV.y) * areaP;
                            }
                        }
                        blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, direction,
                                                      flagForward);
                        num_conn++;
                    }

                    // then find all NNC for current block
                    direction = ConnDirect::y;

                    iznnc = -1;
                    while (upNNC) {
                        // if (-iznnc > k) break;
                        if (-iznnc - static_cast<OCP_INT>(k) > 0) break;
                        // find object block
                        const Hexahedron& frontblock =
                            cornerPoints[i][j + 1][k + iznnc];
                        oindex   = (k + iznnc) * nxny + (j + 1) * nx + i;
                        oFace.p0 = frontblock.p0;
                        oFace.p1 = frontblock.p4;
                        oFace.p2 = frontblock.p5;
                        oFace.p3 = frontblock.p1;

                        SetAllFlags(oFace, Face);

                        // calculate the interface of two face
                        if (flagJump) {
                            // nothing to do
                        } else {
                            if (flagQuad) {
                                areaV = interFace.CalAreaVector();
                            } else {
                                FaceP.p0  = Point3D(Face.p0.x, Face.p0.z, 0);
                                FaceP.p1  = Point3D(Face.p3.x, Face.p3.z, 0);
                                FaceP.p2  = Point3D(Face.p2.x, Face.p2.z, 0);
                                FaceP.p3  = Point3D(Face.p1.x, Face.p1.z, 0);
                                oFaceP.p0 = Point3D(oFace.p0.x, oFace.p0.z, 0);
                                oFaceP.p1 = Point3D(oFace.p3.x, oFace.p3.z, 0);
                                oFaceP.p2 = Point3D(oFace.p2.x, oFace.p2.z, 0);
                                oFaceP.p3 = Point3D(oFace.p1.x, oFace.p1.z, 0);
                                areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                                // attention the direction of vector
                                areaV = Face.CalAreaVector();
                                // correct
                                if (fabs(areaV.y) < 1E-6) {
                                    OCP_WARNING("y is too small");
                                } else {
                                    areaV.x = areaV.x / fabs(areaV.y) * areaP;
                                    areaV.z = areaV.z / fabs(areaV.y) * areaP;
                                    areaV.y = OCP_SIGN(areaV.y) * areaP;
                                }
                            }
                            blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, direction,
                                                          flagForward);
                            num_conn++;
                        }
                        iznnc--;
                    }

                    iznnc = 1;
                    while (downNNC) {
                        if (k + iznnc > nz - 1) break;
                        // find object block
                        const Hexahedron& frontblock =
                            cornerPoints[i][j + 1][k + iznnc];
                        oindex   = (k + iznnc) * nxny + (j + 1) * nx + i;
                        oFace.p0 = frontblock.p0;
                        oFace.p1 = frontblock.p4;
                        oFace.p2 = frontblock.p5;
                        oFace.p3 = frontblock.p1;

                        SetAllFlags(oFace, Face);

                        // calculate the interface of two face
                        if (flagJump) {
                            // nothing to do
                        } else {
                            if (flagQuad) {
                                areaV = interFace.CalAreaVector();
                            } else {
                                FaceP.p0  = Point3D(Face.p0.x, Face.p0.z, 0);
                                FaceP.p1  = Point3D(Face.p3.x, Face.p3.z, 0);
                                FaceP.p2  = Point3D(Face.p2.x, Face.p2.z, 0);
                                FaceP.p3  = Point3D(Face.p1.x, Face.p1.z, 0);
                                oFaceP.p0 = Point3D(oFace.p0.x, oFace.p0.z, 0);
                                oFaceP.p1 = Point3D(oFace.p3.x, oFace.p3.z, 0);
                                oFaceP.p2 = Point3D(oFace.p2.x, oFace.p2.z, 0);
                                oFaceP.p3 = Point3D(oFace.p1.x, oFace.p1.z, 0);
                                areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                                // attention the direction of vector
                                areaV = Face.CalAreaVector();
                                // correct
                                if (fabs(areaV.y) < 1E-6) {
                                    OCP_WARNING("y is too small");
                                } else {
                                    areaV.x = areaV.x / fabs(areaV.y) * areaP;
                                    areaV.z = areaV.z / fabs(areaV.y) * areaP;
                                    areaV.y = OCP_SIGN(areaV.y) * areaP;
                                }
                            }
                            blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, direction,
                                                          flagForward);
                            num_conn++;
                        }
                        iznnc++;
                    }
                }

                //
                // (z-) direction
                //
                direction = ConnDirect::zm;

                Face.p0 = block.p0;
                Face.p1 = block.p3;
                Face.p2 = block.p2;
                Face.p3 = block.p1;
                Pface   = Face.CalCenter();
                Pc2f    = Pface - Pcenter;
                dzpoint = Pc2f;
                if (k == 0) {
                    // nothing to do
                } else {
                    // upblock
                    oindex = (k - 1) * nxny + j * nx + i;

                    interFace = Face;
                    areaV   = interFace.CalAreaVector();
                    blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, direction, flagForward);
                    num_conn++;                   
                }

                //
                // (z+) direction
                //
                direction = ConnDirect::zp;

                Face.p0 = block.p5;
                Face.p1 = block.p6;
                Face.p2 = block.p7;
                Face.p3 = block.p4;
                Pface   = Face.CalCenter();
                Pc2f    = Pface - Pcenter;
                dzpoint = Pc2f - dzpoint;

                if (k == nz - 1) {
                    // nothing to do
                } else {
                    // downblock
                    oindex = (k + 1) * nxny + j * nx + i;

                    interFace = Face;
                    areaV   = interFace.CalAreaVector();
                    blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, direction, flagForward);
                    num_conn++;
                }

                // calculate dx,dy,dz
                dx[cindex] = sqrt(dxpoint.x * dxpoint.x + dxpoint.y * dxpoint.y +
                                  dxpoint.z * dxpoint.z);
                dy[cindex] = sqrt(dypoint.x * dypoint.x + dypoint.y * dypoint.y +
                                  dypoint.z * dypoint.z);
                dz[cindex] = sqrt(dzpoint.x * dzpoint.x + dzpoint.y * dzpoint.y +
                                  dzpoint.z * dzpoint.z);

                OCP_ASSERT(isfinite(dx[cindex]), "Wrong dx!");
                OCP_ASSERT(isfinite(dy[cindex]), "Wrong dy!");
                OCP_ASSERT(isfinite(dz[cindex]), "Wrong dz!");

            }
        }
    }

    OCP_ASSERT(num_conn % 2 == 0, "Wrong Conn!");
    numConnMax = num_conn;
    connect.resize(numConnMax);
    //
    //    calculate the x,y,z direction transmissibilities of each block and save them
    //
    // make the connections
    OCP_ULL iter_conn = 0;
    for (OCP_ULL n = 0; n < numGrid; n++) {
        for (USI j = 0; j < blockconn[n].nConn; j++) {
            OCP_ULL nn = blockconn[n].halfConn[j].neigh;
            USI jj;
            for (jj = 0; jj < blockconn[nn].nConn; jj++) {
                if (blockconn[nn].halfConn[jj].neigh == n) {
                    break;
                }
            }
            if (jj == blockconn[nn].nConn) {
                continue;
            }

            if (blockconn[n].halfConn[j].Ad_dd <= 0 ||
                blockconn[nn].halfConn[jj].Ad_dd <= 0) {
                // false connection
                continue;
            }

            //
            // now, blockconn[n].halfConn[j]
            //     blockconn[nn].halfConn[jj]
            //     are a pair of connections

            connect[iter_conn].begin         = n;
            connect[iter_conn].Ad_dd_begin   = blockconn[n].halfConn[j].Ad_dd;
            connect[iter_conn].end           = nn;
            connect[iter_conn].Ad_dd_end     = blockconn[nn].halfConn[jj].Ad_dd;
            connect[iter_conn].directionType = blockconn[n].halfConn[j].directionType;
            iter_conn++;
        }
    }
    numConn = iter_conn;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/19/2021      Create file                          */
/*----------------------------------------------------------------------------*/