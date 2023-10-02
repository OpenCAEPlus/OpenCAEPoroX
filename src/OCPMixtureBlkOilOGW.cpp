/*! \file    OCPMixtureBlkOilOGW.cpp
 *  \brief   OCPMixtureBlkOilOGW class declaration
 *  \author  Shizhe Li
 *  \date    Jul/19/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPMixtureBlkOilOGW.hpp"


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOGW 
/////////////////////////////////////////////////////


void OCPMixtureBlkOilOGW::Setup(const ParamReservoir& rs_param, const USI& i)
{
	if (rs_param.PVCO_T.data.size() > 0 &&
		rs_param.PVDG_T.data.size() > 0 &&
		rs_param.PVTW_T.data.size() > 0) {
		pmMethod = new OCPMixtureKOGWMethod01(rs_param, i, vs);
	}
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/19/2023      Create file                          */
/*----------------------------------------------------------------------------*/