/*! \file    BulkConnTrans.cpp
 *  \brief   BulkConnTrans class definition
 *  \author  Shizhe Li
 *  \date    Aug/24/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "BulkConnTrans.hpp"


void BulkConnTransMethod01::CalTrans(BulkConnPair& bp, const Bulk& bk)
{
    const BulkVarSet& bcv = bk.GetVarSet();

    OCP_DBL T1, T2;
    const auto& bId   = bp.bId;
    const auto& eId   = bp.eId;
    const auto& areaB = bp.areaB;
    const auto& areaE = bp.areaE;

    switch (bp.direction) 
    {
    case ConnDirect::x:
        T1 = bcv.ntg[bId] * bcv.rockKx[bId] * areaB;
        T2 = bcv.ntg[eId] * bcv.rockKx[eId] * areaE;
        bp.trans = 1 / (1 / T1 + 1 / T2);
        break;
    case ConnDirect::y:
        T1 = bcv.ntg[bId] * bcv.rockKy[bId] * areaB;
        T2 = bcv.ntg[eId] * bcv.rockKy[eId] * areaE;
        bp.trans = 1 / (1 / T1 + 1 / T2);
        break;
    case ConnDirect::z:
        T1 = bcv.rockKz[bId] * areaB;
        T2 = bcv.rockKz[eId] * areaE;
        bp.trans = 1 / (1 / T1 + 1 / T2);
        break;
    case ConnDirect::mf:
        bp.trans = bcv.rockKx[bId] * bcv.v[bId] * bcv.sigma[bId];
        break;
    default:
        OCP_ABORT("Wrong BulkConnType!");
    }
    
}


void BulkConnTrans::Setup()
{
    bcaM = new BulkConnTransMethod01();
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/24/2023      Create file                          */
/*----------------------------------------------------------------------------*/