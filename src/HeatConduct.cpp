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

HeatConductMethod01::HeatConductMethod01(const ParamReservoir& rs_param, HeatConductVarSet& hcvs)
{
    thconP.resize(hcvs.np);
    if (hcvs.o >= 0) thconP[hcvs.o] = rs_param.thcono;
    if (hcvs.g >= 0) thconP[hcvs.g] = rs_param.thcong;
    if (hcvs.w >= 0) thconP[hcvs.w] = rs_param.thconw;

    thconR = rs_param.thconr;


    hcvs.kt.resize(hcvs.nb);
    hcvs.ktP.resize(hcvs.nb);
    hcvs.ktT.resize(hcvs.nb);
    hcvs.ktS.resize(hcvs.nb * hcvs.np);

    hcvs.lkt.resize(hcvs.nb);
    hcvs.lktP.resize(hcvs.nb);
    hcvs.lktT.resize(hcvs.nb);
    hcvs.lktS.resize(hcvs.nb * hcvs.np);

}


void HeatConductMethod01::CalConductCoeff(const OCP_USI& bId, HeatConductVarSet& hcvs, const BulkVarSet& bvs) const
{
    const OCP_USI& np = hcvs.np;

    if (bvs.cType[bId] == BulkContent::rf) {
        // fluid bulk
        OCP_DBL tmp = 0;
        for (USI j = 0; j < np; j++) {
            tmp                    += bvs.S[bId * np + j] * thconP[j];
            hcvs.ktS[bId * np + j]  = bvs.poro[bId] * thconP[j];
        }
        hcvs.kt[bId]  = bvs.poro[bId] * tmp + (1 - bvs.poro[bId]) * thconR;
        hcvs.ktP[bId] = bvs.poroP[bId] * (tmp - thconR);
        hcvs.ktT[bId] = bvs.poroT[bId] * (tmp - thconR);
    }
    else {
        // non fluid bulk
        hcvs.kt[bId]  = thconR;
        hcvs.ktP[bId] = 0;
        hcvs.ktT[bId] = 0;
        for (USI j = 0; j < np; j++) {
            hcvs.ktS[bId * np + j] = 0;
        }
    }
}


OCP_DBL HeatConductMethod01::CalFlux(const HeatConductVarSet& hcvs, const BulkConnPair& bp, const BulkVarSet& bvs) const
{
    const OCP_USI bId = bp.BId();
    const OCP_USI eId = bp.EId();
	const OCP_DBL T1  = hcvs.kt[bId] * bp.AreaB();
	const OCP_DBL T2  = hcvs.kt[eId] * bp.AreaE();
	return (bvs.T[bId] - bvs.T[eId]) / (1 / T1 + 1 / T2);
}


void HeatConductMethod01::AssembleFIM(const BulkConnPair& bp, const HeatConductVarSet& hcvs, const BulkVarSet& bvs, FluxVarSet& fvs) const
{
    const USI& np     = fvs.np;
    const USI& ncol1  = fvs.ncol1;
    const USI& ncol2  = fvs.ncol2;
    auto&      dFdXpB = fvs.dFdXpB;
    auto&      dFdXpE = fvs.dFdXpE;
    auto&      dFdXsB = fvs.dFdXsB;
    auto&      dFdXsE = fvs.dFdXsE;

    const OCP_USI bId   = bp.BId();
    const OCP_USI eId   = bp.EId();
    const OCP_DBL areaB = bp.AreaB();
    const OCP_DBL areaE = bp.AreaE();

	const OCP_DBL T1 = hcvs.kt[bId] * areaB;
	const OCP_DBL T2 = hcvs.kt[eId] * areaE;
	const OCP_DBL Adkt = 1 / (1 / T1 + 1 / T2);
	const OCP_DBL tmpB = pow(Adkt, 2) / pow(T1, 2) * areaB;
	const OCP_DBL tmpE = pow(Adkt, 2) / pow(T2, 2) * areaE;
	const OCP_DBL dT = bvs.T[bId] - bvs.T[eId];
	// Thermal Conduction
	// dP
	dFdXpB[ncol1 * ncol1 - ncol1] += tmpB * hcvs.ktP[bId] * dT;
	dFdXpE[ncol1 * ncol1 - ncol1] += tmpE * hcvs.ktP[eId] * dT;
	// dT
	dFdXpB[ncol1 * ncol1 - 1] += Adkt + tmpB * hcvs.ktT[bId] * dT;
	dFdXpE[ncol1 * ncol1 - 1] += -Adkt + tmpE * hcvs.ktT[eId] * dT;
	// dS
	for (OCP_USI j = 0; j < np; j++) {
		dFdXsB[(ncol1 - 1) * ncol2 + j] += tmpB * hcvs.ktS[bId * np + j] * dT;
		dFdXsE[(ncol1 - 1) * ncol2 + j] += tmpE * hcvs.ktS[eId * np + j] * dT;
	}
}


void HeatConduct::Setup(const ParamReservoir& rs_param, const BulkVarSet& bvs)
{
    ifUse = rs_param.ifThcon;
    if (ifUse) {
        if (rs_param.thermal == OCP_FALSE) {
            ifUse = OCP_FALSE;
            OCP_WARNING("HEATCONDUCT is IGNORED in ISOTHERMAL MODEL!");
            return;
        }
        vs.SetUp(bvs.nb, bvs.np, bvs.o, bvs.g, bvs.w);
        hcM.push_back(new HeatConductMethod01(rs_param, vs));
    }
}


void HeatConduct::CalConductCoeff(const BulkVarSet& bvs)
{ 
    if (ifUse) {
        for (OCP_USI n = 0; n < bvs.nb; n++) {
            hcM[0]->CalConductCoeff(n, vs, bvs);
        }
    }
}


void HeatConduct::CalFlux(const BulkConnPair& bp, const BulkVarSet& bvs)
{
    if (ifUse) {
        vs.conduct_H = hcM[0]->CalFlux(vs, bp, bvs);
    }
    else {
        vs.conduct_H = 0;
    }
}


void HeatConduct::AssembleFIM(const BulkConnPair& bp, const BulkVarSet& bvs, FluxVarSet& fvs) const 
{
    if (ifUse) {
        hcM[0]->AssembleFIM(bp, vs, bvs, fvs);
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/26/2023      Create file                          */
/*----------------------------------------------------------------------------*/