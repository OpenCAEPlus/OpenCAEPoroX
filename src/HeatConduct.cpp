/*! \file    HeatConduct.cpp
 *  \brief   HeatConduct class definition
 *  \author  Shizhe Li
 *  \date    Aug/26/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "HeatConduct.hpp"

HeatConductMethod01::HeatConductMethod01(const ParamReservoir& rs_param, HeatConductVarSet& hlvs)
{
    thconp.resize(hlvs.np);
    if (hlvs.o >= 0) thconp[hlvs.o] = rs_param.thcono;
    if (hlvs.g >= 0) thconp[hlvs.g] = rs_param.thcong;
    if (hlvs.w >= 0) thconp[hlvs.w] = rs_param.thconw;

    thconr = rs_param.thconr;


    hlvs.kt.resize(hlvs.nb);
    hlvs.ktP.resize(hlvs.nb);
    hlvs.ktT.resize(hlvs.nb);
    hlvs.ktS.resize(hlvs.nb * hlvs.np);

    hlvs.lkt.resize(hlvs.nb);
    hlvs.lktP.resize(hlvs.nb);
    hlvs.lktT.resize(hlvs.nb);
    hlvs.lktS.resize(hlvs.nb * hlvs.np);

}


void HeatConductMethod01::CalHeatConduct(const OCP_USI& bId, HeatConductVarSet& hlvs, const BulkVarSet& bvs) const
{
    const OCP_USI& np = hlvs.np;

    if (bvs.cType[bId] == BulkContent::rf) {
        // fluid bulk
        OCP_DBL tmp = 0;
        for (USI j = 0; j < np; j++) {
            tmp                    += bvs.S[bId * np + j] * thconp[j];
            hlvs.ktS[bId * np + j]  = bvs.poro[bId] * thconp[j];
        }
        hlvs.kt[bId]  = bvs.poro[bId] * tmp + (1 - bvs.poro[bId]) * thconr;
        hlvs.ktP[bId] = bvs.poroP[bId] * (tmp - thconr);
        hlvs.ktT[bId] = bvs.poroT[bId] * (tmp - thconr);
    }
    else {
        // non fluid bulk
        hlvs.kt[bId]  = thconr;
        hlvs.ktP[bId] = 0;
        hlvs.ktT[bId] = 0;
        for (USI j = 0; j < np; j++) {
            hlvs.ktS[bId * np + j] = 0;
        }
    }
}


void HeatConduct::Setup(const ParamReservoir& rs_param, const BulkVarSet& bvs)
{
    ifUse = rs_param.ifThcon;
    if (ifUse) {
        vs.SetNb(bvs.nb, bvs.np, bvs.oIndex, bvs.gIndex, bvs.wIndex);
        hcM.push_back(new HeatConductMethod01(rs_param, vs));
    }
}


void HeatConduct::CalHeatConduct(const BulkVarSet& bvs)
{ 
    if (ifUse) {
        for (OCP_USI n = 0; n < bvs.nb; n++) {
            hcM[0]->CalHeatConduct(n, vs, bvs);
        }
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/26/2023      Create file                          */
/*----------------------------------------------------------------------------*/