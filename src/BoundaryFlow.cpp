/*! \file    BoundaryFlow.cpp
 *  \brief   BoundaryFlow class declaration
 *  \author  Shizhe Li
 *  \date    Sep/27/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


#include "BoundaryFlow.hpp"


void BoundaryFlow::Setup(const ParamReservoir& rs_param, const BulkVarSet& bvs, const vector<string>& boundName, const vector<USI>& boundIndex)
{

	if (rs_param.BDparam.size() > 0) {
		ifUse = OCP_TRUE;
		mIndex.resize(bvs.nbI);

		for (const auto& bP : rs_param.BDparam) {
			if (bP.constP) {
				bfM.push_back(new BoundaryFlowMethod01(bP, vs));
			}
		}

		const USI bLen = boundName.size();
		vector<USI> boundCode(bLen, 0);
		for (OCP_USI n = 0; n < bvs.nbI; n++) {
			fill(boundCode.begin(), boundCode.end(), 0.0);
			for (OCP_INT i = bLen - 1; i >= 0; i--) {
				boundCode[i] = (((boundIndex[n] >> i) & 1));
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/27/2023      Create file                          */
/*----------------------------------------------------------------------------*/
