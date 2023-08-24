/*! \file    BulkConnArea.cpp
 *  \brief   BulkConnArea class definition
 *  \author  Shizhe Li
 *  \date    Aug/24/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "BulkConnArea.hpp"


void BulkConnAreaMethod01::CalArea(BulkConnPair& bp, const Bulk& bk)
{
    const BulkVarSet& bcv = bk.GetVarSet();

    OCP_DBL T1, T2;
    const auto& bId   = bp.bId;
    const auto& eId   = bp.eId;
    const auto& areaB = bp.areaB;
    const auto& areaE = bp.areaE;

    switch (bp.direction) 
    {
    case 1:
        T1 = bcv.ntg[bId] * bcv.rockKx[bId] * areaB;
        T2 = bcv.ntg[eId] * bcv.rockKx[eId] * areaE;
        break;
    case 2:
        T1 = bcv.ntg[bId] * bcv.rockKy[bId] * areaB;
        T2 = bcv.ntg[eId] * bcv.rockKy[eId] * areaE;
        break;
    case 3:
        T1 = bcv.rockKz[bId] * areaB;
        T2 = bcv.rockKz[eId] * areaE;
        break;
    default:
        OCP_ABORT("Wrong Direction!");
    }
    bp.area = 1 / (1 / T1 + 1 / T2);
}


void BulkConnArea::Setup()
{
    bcaM = new BulkConnAreaMethod01();
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/24/2023      Create file                          */
/*----------------------------------------------------------------------------*/