/*! \file    OCPMixtureBlkOilOW.cpp
 *  \brief   OCPMixtureBlkOilOW class declaration
 *  \author  Shizhe Li
 *  \date    Jul/13/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPMixtureBlkOilOW.hpp"


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOW 
/////////////////////////////////////////////////////

void OCPMixtureBlkOilOW::Setup(const ParamReservoir& rs_param, const USI& i)
{  
    vs.Init(2, 2);
    GetStdRhoOW(rs_param);
	if (rs_param.PVTW_T.data.size() > 0 && rs_param.PVDO_T.data.size() > 0) {
		pmMethod = new OCPMixtureBlkOilOWMethod01(rs_param.PVDO_T.data[i], 
                                                  rs_param.PVTW_T.data[i], 
                                                  stdRhoO, stdRhoW, &vs);
	}
}


void OCPMixtureBlkOilOW::GetStdRhoOW(const ParamReservoir& rs_param)
{
    if (rs_param.density.activity) {
        stdRhoO = rs_param.density.data[0];
        stdRhoW = rs_param.density.data[1];
    }
    else {
        stdRhoO = (141.5 * RHOW_STD) / (rs_param.gravity.data[0] + 131.5);
        stdRhoW = RHOW_STD * rs_param.gravity.data[1];
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/13/2023      Create file                          */
/*----------------------------------------------------------------------------*/