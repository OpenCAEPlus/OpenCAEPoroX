/*! \file    UtilMath.cpp
 *  \brief   UtilMath class declaration
 *  \author  Shizhe Li
 *  \date    Jul/23/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "UtilMath.hpp"


USI CubicRoot(const OCP_DBL& a, const OCP_DBL& b, const OCP_DBL& c, const OCP_BOOL& NTflag, vector<OCP_DBL>& z)
{

    OCP_DBL Q = (a * a - 3 * b) / 9;
    OCP_DBL R = (2 * a * a * a - 9 * a * b + 27 * c) / 54;

    OCP_DBL Q3 = Q * Q * Q;
    OCP_DBL M = R * R - Q3;

    if (M <= 0) {
        // 3 real roots
        OCP_DBL theta = acos(R / sqrt(Q3));
        z[0] = -2 * sqrt(Q) * cos(theta / 3) - a / 3;
        z[1] = -2 * sqrt(Q) * cos((theta + 2 * PI) / 3) - a / 3;
        z[2] = -2 * sqrt(Q) * cos((theta - 2 * PI) / 3) - a / 3;

        if (NTflag) {
            NTcubicroot(z[0], a, b, c);
            NTcubicroot(z[1], a, b, c);
            NTcubicroot(z[2], a, b, c);
        }

        sort(z.begin(), z.end());

        return 3;
    }
    else {
        OCP_DBL tmp1 = -R + sqrt(M);
        OCP_DBL tmp2 = R + sqrt(M);
        OCP_DBL S = signD(tmp1) * pow(fabs(tmp1), 1.0 / 3);
        OCP_DBL T = -signD(tmp2) * pow(fabs(tmp2), 1.0 / 3);
        z[0] = S + T - a / 3;

        if (NTflag) {
            NTcubicroot(z[0], a, b, c);
        }
        return 1;
    }
}


void NTcubicroot(OCP_DBL& root, const OCP_DBL& a, const OCP_DBL& b, const OCP_DBL& c)
{
    OCP_DBL e = root * (root * (root + a) + b) + c;
    OCP_DBL df;
    OCP_DBL iter = 0;
    OCP_DBL optroot = root;
    OCP_DBL opte = fabs(e);

    while (fabs(e) > 1E-8) {

        df = root * (3 * root + 2 * a) + b;
        root = root - e / df;
        iter++;
        if (iter > 10) {
            // std::cout << "WARNING: INEXACT ROOT FOR CUBIC EQUATIONS" << std::endl;
            break;
        }
        e = root * (root * (root + a) + b) + c;
        if (fabs(e) <= opte) {
            opte = fabs(e);
            optroot = root;
        }
    }
    root = optroot;
}


/// Return the sign of double d
OCP_DBL signD(const OCP_DBL& d)
{
    if (d > 0) {
        return 1.0;
    }
    else if (d < 0) {
        return -1.0;
    }
    else {
        return 0.0;
    }
}

OCP_DBL delta(const USI& i, const USI& j)
{
    if (i == j) {
        return 1.0;
    }
    return 0.0;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/23/2023      Create file                          */
/*----------------------------------------------------------------------------*/