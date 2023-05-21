/*! \file    BulkConn.cpp
 *  \brief   BulkConn class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// Standard header files
#include <cassert>
#include <cmath>
#include <ctime>

// OpenCAEPoroX header files
#include "BulkConn.hpp"

/////////////////////////////////////////////////////////////////////
// General
/////////////////////////////////////////////////////////////////////

void BulkConn::SetupIsoT(const Bulk& bk)
{
    CalAkd(bk);
    SetConnTypeIsoT(bk);
}


void BulkConn::SetupT(const Bulk& bk)
{
    CalAkd(bk);
    SetConnTypeT(bk);
}


void BulkConn::CalAkd(const Bulk& bk)
{
    OCP_USI bId, eId;
    OCP_DBL areaB, areaE;
    OCP_DBL T1, T2;
    for (OCP_USI c = 0; c < numConn; c++) {
        bId   = iteratorConn[c].bId;
        eId   = iteratorConn[c].eId;
        areaB = iteratorConn[c].areaB;
        areaE = iteratorConn[c].areaE;
        switch (iteratorConn[c].direction) {
            case 1:
                T1 = bk.ntg[bId] * bk.rockKx[bId] * areaB;
                T2 = bk.ntg[eId] * bk.rockKx[eId] * areaE;
                break;
            case 2:
                T1 = bk.ntg[bId] * bk.rockKy[bId] * areaB;
                T2 = bk.ntg[eId] * bk.rockKy[eId] * areaE;
                break;
            case 3:
                T1 = bk.rockKz[bId] * areaB;
                T2 = bk.rockKz[eId] * areaE;
                break;
            default:
                OCP_ABORT("Wrong Direction!");
        }
        iteratorConn[c].area = 1 / (1 / T1 + 1 / T2);
    }
}


void BulkConn::SetConnTypeIsoT(const Bulk& bk)
{           
    flux.push_back(new OCPFlux_IsoT(bk));
    for (auto& iter : iteratorConn) {
        iter.type = 0;
    }
}


void BulkConn::SetConnTypeT(const Bulk& bk)
{   
    flux.push_back(new OCPFlux_T(bk));

    for (auto& iter : iteratorConn) {
        iter.type = 0;
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/