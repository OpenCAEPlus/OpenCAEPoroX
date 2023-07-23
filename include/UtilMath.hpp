/*! \file    UtilMath.hpp
 *  \brief   UtilMath class declaration
 *  \author  Shizhe Li
 *  \date    Jul/23/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __UTILMATH_HEADER__
#define __UTILMATH_HEADER__

#include "OCPConst.hpp"
#include <vector>
#include <algorithm>

using namespace std;

USI CubicRoot(const OCP_DBL& a, const OCP_DBL& b, const OCP_DBL& c, const OCP_BOOL& NTflag, vector<OCP_DBL>& z);

void NTcubicroot(OCP_DBL& root, const OCP_DBL& a, const OCP_DBL& b, const OCP_DBL& c);

OCP_DBL signD(const OCP_DBL& d);

OCP_DBL delta(const USI& i, const USI& j);

#endif /* end if __UTILMATH_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/23/2023      Create file                          */
/*----------------------------------------------------------------------------*/