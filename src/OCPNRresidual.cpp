/*! \file    OCPNRresidual.cpp
 *  \brief   data structure used in non-linear iterations
 *  \author  Shizhe Li
 *  \date    Oct/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


#include "OCPNRresidual.hpp"


void OCPNRresidual::SetupIsoT(const OCP_USI& nb, const OCP_USI& nw, const OCP_USI& nc)
{
    OCP_USI reslen = (nb + nw) * (nc + 1);
    resAbs.resize(reslen);
    resRelV.resize(nb);
    resRelN.resize(nb);
}


void OCPNRresidual::SetupT(const OCP_USI& nb, const OCP_USI& nw, const OCP_USI& nc)
{
    OCP_USI reslen = (nb + nw) * (nc + 2);
    resAbs.resize(reslen);
    resRelV.resize(nb);
    resRelN.resize(nb);
    resRelE.resize(nb);
}


void OCPNRresidual::SetZero()
{
    fill(resAbs.begin(), resAbs.end(), 0);
    fill(resRelV.begin(), resRelV.end(), 0);
    fill(resRelN.begin(), resRelN.end(), 0);
    fill(resRelE.begin(), resRelE.end(), 0);
    maxRelRes_V = 0;
    maxRelRes_N = 0;
    maxRelRes_E = 0;
    maxWellRelRes_mol = 0;
    maxId_V = 0;
    maxId_N = 0;
    maxId_E = 0;
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/