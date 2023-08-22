/*! \file    ROCKModule.cpp
 *  \brief   ROCKModule class definition
 *  \author  Shizhe Li
 *  \date    Aug/21/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "ROCKModule.hpp"

void ROCKModule::Setup(const ParamReservoir& rs_param, const BulkVarSet& bvs, OptionalFeatures& opts)
{
	NTROCK = rs_param.NTROOC;

	if (rs_param.thermal) {
		for (USI i = 0; i < NTROCK; i++) {
			if (rs_param.rockSet[i].type == "LINEAR") {
				ROCKs.push_back(new OCPRockT_Linear(rs_param.rockSet[i]));
			}
			else {
				ROCKs.push_back(new OCPRockT_Exp(rs_param.rockSet[i]));
			}
		}
	}
	else {
		for (USI i = 0; i < NTROCK; i++) {
			ROCKs.push_back(new OCPRockIsoT_Linear(rs_param.rockSet[i]));
		}
	}

	if (ROCKNUM.empty()) ROCKNUM.resize(bvs.nb, 0);
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/21/2023      Create file                          */
/*----------------------------------------------------------------------------*/