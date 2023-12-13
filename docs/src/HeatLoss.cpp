/*! \file    HeatLoss.cpp
 *  \brief   HeatLoss class definition
 *  \author  Shizhe Li
 *  \date    Aug/24/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "HeatLoss.hpp"


HeatLossMethod01::HeatLossMethod01(const OCP_DBL& bK_in, const OCP_DBL& bC_in, HeatLossVarSet& hlvs)
{
	bK = bK_in;
	bD = bK_in / bC_in;

	hlvs.I.resize(hlvs.nb);
	hlvs.hl.resize(hlvs.nb);
	hlvs.hlT.resize(hlvs.nb);
	hlvs.lI.resize(hlvs.nb);
	hlvs.lhl.resize(hlvs.nb);
	hlvs.lhlT.resize(hlvs.nb);
}


void HeatLossMethod01::CalHeatLoss(const OCP_USI& bId, HeatLossVarSet& hlvs, const BulkVarSet& bvs, const OCP_DBL& t, const OCP_DBL& dt)
{

	const OCP_DBL lambda = bD;
	const OCP_DBL kappa  = bK;
	
	const OCP_DBL dT    = bvs.T[bId] - bvs.lT[bId];
	const OCP_DBL theta = bvs.T[bId] - bvs.initT[bId];
	const OCP_DBL d     = sqrt(lambda * t) / 2;
	const OCP_DBL tmp   = 3 * pow(d, 2) + lambda * dt;
	const OCP_DBL p     = (theta * (lambda * dt / d) + hlvs.lI[bId] - dT * (pow(d, 3) / (lambda * dt))) / tmp;
	const OCP_DBL pT    = (lambda * dt / d - pow(d, 3) / (lambda * dt)) / tmp;
	const OCP_DBL q     = (2 * p * d - theta + pow(d, 2) * dT / (lambda * dt)) / (2 * pow(d, 2));
	
	hlvs.I[bId]   = theta * d + p * pow(d, 2) + 2 * q * pow(d, 3);
	hlvs.hl[bId]  = kappa * (2 * (bvs.T[bId] - bvs.initT[bId]) / sqrt(lambda * t) - p) * (bvs.dx[bId] * bvs.dy[bId]);
	hlvs.hlT[bId] = kappa * (2 / sqrt(lambda * t) - pT) * (bvs.dx[bId] * bvs.dy[bId]);
}


void HeatLoss::Setup(const ParamReservoir& rs_param, const BulkVarSet& bvs, const vector<USI>& boundIndex)
{
	ifUse = rs_param.hLoss.ifHLoss;
	if (ifUse) {
		if (rs_param.thermal == OCP_FALSE) {
			ifUse = OCP_FALSE;
			OCP_WARNING("HEATLOSS is IGNORED in ISOTHERMAL MODEL!");
			return;
		}
		vs.SetNb(bvs.nbI);

		INT obIndex = -1;
		if (rs_param.hLoss.obUse) {
			// use overburden heatloss
			hlM.push_back(new HeatLossMethod01(rs_param.hLoss.obK, rs_param.hLoss.obC, vs));
			obIndex = hlM.size() - 1;
		}
		INT ubIndex = -1;
		if (rs_param.hLoss.ubUse) {
			// use underburden heatloss
			hlM.push_back(new HeatLossMethod01(rs_param.hLoss.ubK, rs_param.hLoss.ubC, vs));
			ubIndex = hlM.size() - 1;
		}
		mIndex.resize(vs.nb, -1);
		for (OCP_USI n = 0; n < vs.nb; n++) {
			if (boundIndex[n] == 1 && rs_param.hLoss.obUse) {
				mIndex[n] = obIndex;
			}
			else if (boundIndex[n] == 2 && rs_param.hLoss.ubUse) {
				mIndex[n] = ubIndex;
			}
		}
	}
}


void HeatLoss::CalHeatLoss(const BulkVarSet& bvs, const OCP_DBL& t, const OCP_DBL& dt)
{
	if (ifUse) {
		for (OCP_USI n = 0; n < vs.nb; n++) {
			if (mIndex[n] >= 0) {
				hlM[mIndex[n]]->CalHeatLoss(n, vs, bvs, t, dt);
			}			
		}
	}
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/24/2023      Create file                          */
/*----------------------------------------------------------------------------*/