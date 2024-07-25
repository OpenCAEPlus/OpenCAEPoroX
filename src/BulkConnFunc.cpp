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

#include "BulkConnFunc.hpp"


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
    case ConnDirect::xp:
    case ConnDirect::xm:
        T1 = bcv.ntg[bId] * bcv.rockKx[bId] * areaB;
        T2 = bcv.ntg[eId] * bcv.rockKx[eId] * areaE;
        bp.trans = bp.transMult / (1 / T1 + 1 / T2);
        break;
    case ConnDirect::y:
    case ConnDirect::yp:
    case ConnDirect::ym:
        T1 = bcv.ntg[bId] * bcv.rockKy[bId] * areaB;
        T2 = bcv.ntg[eId] * bcv.rockKy[eId] * areaE;
        bp.trans = bp.transMult / (1 / T1 + 1 / T2);
        break;
    case ConnDirect::z:
    case ConnDirect::zp:
    case ConnDirect::zm:
        T1 = bcv.rockKz[bId] * areaB;
        T2 = bcv.rockKz[eId] * areaE;
        bp.trans = bp.transMult / (1 / T1 + 1 / T2);
        break;
    case ConnDirect::mf:
    case ConnDirect::fm:
        bp.trans = bp.transMult * bcv.rockKx[bId] * bcv.v[bId] * bcv.sigma[bId];
        break;
    case ConnDirect::usg:
        // homogeneous now
        T1 = bcv.ntg[bId] * bcv.rockKx[bId] * areaB;
        T2 = bcv.ntg[eId] * bcv.rockKx[eId] * areaE;
        bp.trans = bp.transMult / (1 / T1 + 1 / T2);
        break;
    default:
        OCP_ABORT("Wrong BulkConnType!");
    }

    if (!isfinite(bp.trans)) {
        OCP_ABORT("Transmissbility is NAN!");
    }

    bp.trans *= CONV_DARCY;
}


void BulkConnTrans::Setup()
{
    bctM = new BulkConnTransMethod01();
}



void BulkConnDiffuMethod01::CalDiffu(BulkConnPair& bp, const Bulk& bk)
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
    case ConnDirect::xp:
    case ConnDirect::xm:
        //T1 = bcv.ntg[bId] * bcv.poro[bId] * areaB;
        //T2 = bcv.ntg[eId] * bcv.poro[eId] * areaE;
        T1 = bcv.poro[bId] * areaB;
        T2 = bcv.poro[eId] * areaE;
        bp.diffu = 1 / (1 / T1 + 1 / T2);
        break;
    case ConnDirect::y:
    case ConnDirect::yp:
    case ConnDirect::ym:
        //T1 = bcv.ntg[bId] * bcv.poro[bId] * areaB;
        //T2 = bcv.ntg[eId] * bcv.poro[eId] * areaE;
        T1 = bcv.poro[bId] * areaB;
        T2 = bcv.poro[eId] * areaE;
        bp.diffu = 1 / (1 / T1 + 1 / T2);
        break;
    case ConnDirect::z:
    case ConnDirect::zp:
    case ConnDirect::zm:
        T1 = bcv.poro[bId] * areaB;
        T2 = bcv.poro[eId] * areaE;
        bp.diffu = 1 / (1 / T1 + 1 / T2);
        break;
    case ConnDirect::mf:
    case ConnDirect::fm:
        bp.diffu = bcv.poro[bId] * bcv.v[bId] * bcv.sigma[bId];
        break;
    case ConnDirect::usg:
        // homogeneous now
        //T1 = bcv.ntg[bId] * bcv.poro[bId] * areaB;
        //T2 = bcv.ntg[eId] * bcv.poro[eId] * areaE;
        T1 = bcv.poro[bId] * areaB;
        T2 = bcv.poro[eId] * areaE;
        bp.diffu = 1 / (1 / T1 + 1 / T2);
        break;
    default:
        OCP_ABORT("Wrong BulkConnType!");
    }

    if (!isfinite(bp.diffu)) {
        OCP_ABORT("Transmissbility is NAN!");
    }
}


void BulkConnDiffu::Setup()
{
    bcdM = new BulkConnDiffuMethod01();
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/24/2023      Create file                          */
/*----------------------------------------------------------------------------*/