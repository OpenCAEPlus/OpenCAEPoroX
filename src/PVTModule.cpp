/*! \file    PVTModule.cpp
 *  \brief   PVTModule class definition
 *  \author  Shizhe Li
 *  \date    Aug/19/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "PVTModule.hpp"

void PVTModule::Setup(const ParamReservoir& rs_param, BulkVarSet& bvs, OptionalFeatures& opts)
{
	NTPVT = rs_param.NTPVT;

	if (rs_param.thermal) {
		for (USI i = 0; i < NTPVT; i++)
			PVTs.push_back(new MixtureUnitThermal_OW(rs_param, i, opts));
	}
	else if (rs_param.blackOil) {
        if (rs_param.water && rs_param.oil && !rs_param.gas) {
            for (USI i = 0; i < NTPVT; i++)
                PVTs.push_back(new MixtureUnitBlkOil_OW(rs_param, i, opts));
        }
        else if (rs_param.water && rs_param.oil && rs_param.gas) {
            for (USI i = 0; i < NTPVT; i++)
                PVTs.push_back(new MixtureUnitBlkOil_OGW(rs_param, i, opts));
		}
		else {
			OCP_ABORT("Inavailable Mixture Type!");
		}
	}
	else if (rs_param.comps) {
		for (USI i = 0; i < NTPVT; i++) 
			PVTs.push_back(new MixtureUnitComp(rs_param, i, opts));
	}

    mixType     = PVTs[0]->GetMixtureType();

	if (PVTNUM.empty()) PVTNUM.resize(bvs.nb, 0);

	bvs.np      = PVTs[0]->GetVs()->np;
	bvs.nc      = PVTs[0]->GetVs()->nc;
	bvs.oIndex  = PVTs[0]->GetVs()->o;
	bvs.gIndex  = PVTs[0]->GetVs()->g;
	bvs.wIndex  = PVTs[0]->GetVs()->w;


}



 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Aug/19/2023      Create file                          */
 /*----------------------------------------------------------------------------*/