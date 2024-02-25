/*! \file    BulkAccumuModule.cpp
 *  \brief   BulkAccumuModule class definition
 *  \author  Shizhe Li
 *  \date    Aug/30/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "BulkAccumuModule.hpp"


BulkAccumuTerm01::BulkAccumuTerm01(const BulkVarSet& bvs, const OptionalModules* opt)
{
	dim = bvs.nc + 1;
	dFdXp.resize(dim * dim);
	res.resize(dim);
	optM = opt;
}


const vector<OCP_DBL>& BulkAccumuTerm01::CalResFIM(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt) const
{
    fill(res.begin(), res.end(), 0.0);

	res[0] = bvs.rockVp[bId] - bvs.vf[bId];
	for (USI i = 0; i < bvs.nc; i++) {
		res[i + 1] = bvs.Ni[bId * bvs.nc + i] - bvs.lNi[bId * bvs.nc + i];
	}

    // for spe11 now
    if (optM->boundary.boundaryFlow.IfUse(bId)) {
        const OCP_DBL dP = (bvs.Pj[bId * 2 + 0] - GRAVITY_FACTOR * bvs.rho[bId * 2 + 0] * bvs.depth[bId]) -
            (PRESSURE_STD - GRAVITY_FACTOR * bvs.rho[bId * 2 + 0] * 0);
        const OCP_DBL tmp = dt * CONV_DARCY * optM->boundary.boundArea[bId] * bvs.rockKx[bId] * bvs.kr[bId * 2 + 0] / bvs.mu[bId * 2 + 0] * dP;
        res[1 + bvs.w] += tmp;
        
        // cout << scientific << setprecision(10) << dP << "   " << bvs.Ni[bId * 2] << "   " << tmp << endl;
    }

	return res;
}


const vector<OCP_DBL>& BulkAccumuTerm01::CaldFdXpFIM(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt) const
{
    fill(dFdXp.begin(), dFdXp.end(), 0.0);

	dFdXp[0] = bvs.v[bId] * bvs.poroP[bId] - bvs.vfP[bId];
	for (USI i = 0; i < bvs.nc; i++) {
		dFdXp[i + 1] = -bvs.vfi[bId * bvs.nc + i];
	}
	for (USI i = 1; i < dim; i++) {
		dFdXp[i * dim + i] = 1;
	}

    if (optM->boundary.boundaryFlow.IfUse(bId)) {
        const OCP_DBL dP = (bvs.Pj[bId * 2 + 0] - GRAVITY_FACTOR * bvs.rho[bId * 2 + 0] * bvs.depth[bId]) -
            (PRESSURE_STD - GRAVITY_FACTOR * bvs.rho[bId * 2 + 0] * 0);
        OCP_DBL tmp = dt * CONV_DARCY * optM->boundary.boundArea[bId]
            * bvs.kr[bId * 2 + 0] / bvs.mu[bId * 2 + 0] * (1.0 - GRAVITY_FACTOR * bvs.rhoP[bId * 2 + 0] * bvs.depth[bId]
                + GRAVITY_FACTOR * bvs.rhoP[bId * 2 + 0] * bvs.depth[bId] * 0);
        tmp += -dt * CONV_DARCY * optM->boundary.boundArea[bId] * bvs.kr[bId * 2 + 0] * dP * bvs.muP[bId * 2 + 0] / (bvs.mu[bId * 2 + 0] * bvs.mu[bId * 2 + 0]);
        dFdXp[(1 + bvs.w) * dim + 0] += tmp;
    }

	return dFdXp;
}


void BulkAccumuTerm01::CalValRhsIMPEC(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt, OCP_DBL& valout, OCP_DBL& rhsout) const
{
	valout = bvs.v[bId] * bvs.poroP[bId] - bvs.vfP[bId];
	rhsout = (bvs.v[bId] * bvs.poroP[bId] - bvs.vfP[bId]) * bvs.P[bId] + dt * (bvs.vf[bId] - bvs.rockVp[bId]);
}


BulkAccumuTerm02::BulkAccumuTerm02(const BulkVarSet& bvs, const OptionalModules* opt)
{
	dim = bvs.nc + 2;
	dFdXp.resize(dim * dim);
	res.resize(dim);
	optM = opt;
}


const vector<OCP_DBL>& BulkAccumuTerm02::CalResFIM(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt) const
{
    fill(res.begin(), res.end(), 0.0);

    const USI& nc = bvs.nc;
    if (bvs.cType[bId] == BulkContent::rf) {
        // Fluid bulk
        // Volume Conservation
        res[0] = bvs.rockVp[bId] - bvs.vf[bId];
        // Mass Conservation
        for (USI i = 0; i < nc; i++) {
            res[1 + i] = bvs.Ni[bId * nc + i] - bvs.lNi[bId * nc + i];
        }
        // Energy Conservation
        res[nc + 1] = (bvs.vf[bId] * bvs.Uf[bId] + bvs.vr[bId] * bvs.Hr[bId]) -
                      (bvs.lvf[bId] * bvs.lUf[bId] + bvs.lvr[bId] * bvs.lHr[bId]);
    }
    else {
        // only rock
        res[nc + 1] = bvs.vr[bId] * bvs.Hr[bId] - bvs.lvr[bId] * bvs.lHr[bId];
    }
    
    // Heat Loss
    if (optM->boundary.heatLoss.IfUse(bId)) {
        // dT
        res[nc + 1] += dt * optM->boundary.heatLoss.GetHl(bId);
    }
    return res;
}


const vector<OCP_DBL>& BulkAccumuTerm02::CaldFdXpFIM(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt) const
{

    fill(dFdXp.begin(), dFdXp.end(), 0.0);

    const USI& nc = bvs.nc;
    if (bvs.cType[bId] == BulkContent::rf) {
        // Volume consevation
        // dP
        dFdXp[0] = bvs.v[bId] * bvs.poroP[bId] - bvs.vfP[bId];
        // dNi
        for (USI i = 0; i < nc; i++) {
            dFdXp[i + 1] = -bvs.vfi[bId * nc + i];
        }
        // dT
        dFdXp[nc + 1] = bvs.v[bId] * bvs.poroT[bId] - bvs.vfT[bId];
        // Mass consevation
        // dNi
        for (USI i = 1; i < nc + 1; i++) {
            dFdXp[i * dim + i] = 1;
        }
        // Energy consevation
        // dP
        dFdXp[dim * (dim - 1)] = bvs.vfP[bId] * bvs.Uf[bId] + bvs.vf[bId] * bvs.UfP[bId] + bvs.vrP[bId] * bvs.Hr[bId];
        // dNi
        for (USI i = 0; i < nc; i++) {
            dFdXp[dim * (dim - 1) + i + 1] = bvs.vfi[bId * nc + i] * bvs.Uf[bId] + bvs.vf[bId] * bvs.Ufi[bId * nc + i];
        }
        // dT
        dFdXp[dim * dim - 1] = bvs.vfT[bId] * bvs.Uf[bId] + bvs.vf[bId] * bvs.UfT[bId] + bvs.vrT[bId] * bvs.Hr[bId] + bvs.vr[bId] * bvs.HrT[bId];
    }
    else {
        // only rock
        // Energy consevation
        // dT
        dFdXp[dim * dim - 1] = bvs.vrT[bId] * bvs.Hr[bId] + bvs.vr[bId] * bvs.HrT[bId];
    }

    // Heat Loss iterm
    if (optM->boundary.heatLoss.IfUse(bId)) {
        // dT
        dFdXp[dim * dim - 1] += dt * optM->boundary.heatLoss.GetHlT(bId);
    }

    return dFdXp;
}


/// Setup accumulation module
void BulkAccumuModule::Setup(const ParamReservoir& param, const BulkVarSet& bvs, const OptionalModules& opt)
{
    if (param.thermal) bacT = new BulkAccumuTerm02(bvs, &opt);
    else               bacT = new BulkAccumuTerm01(bvs, &opt);
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/30/2023      Create file                          */
/*----------------------------------------------------------------------------*/