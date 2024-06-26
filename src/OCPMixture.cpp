/*! \file    OCPMixture.cpp
 *  \brief   OCPMixture class declaration
 *  \author  Shizhe Li
 *  \date    Oct/02/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPMixture.hpp"

using namespace std;



////////////////////////////////////////////////////////////////
// OCPMixtureComp 
////////////////////////////////////////////////////////////////


OCPMixtureComp::OCPMixtureComp(const ParamReservoir& rs_param, const USI& i, BulkOptionalModules& opts)
{
    if (rs_param.PVTW_T.data[i].size() > 0) {
        pmMethod = new OCPMixtureMethodComp01(rs_param, i, vs);
    }

    // Skip stability analysis
    skipPSA = &opts.skipPSA;
    skipMethodIndex = skipPSA->Setup(opts.nb, pmMethod);
}



/////////////////////////////////////////////////////
// OCPMixtureK 
/////////////////////////////////////////////////////


OCPMixtureK::OCPMixtureK(const ParamReservoir& rs_param, const USI& i, BulkOptionalModules& opts)
{
    if (rs_param.PVCO_T.data.size() > 0 &&
        rs_param.PVDG_T.data.size() > 0 &&
        rs_param.PVTW_T.data.size() > 0) {
        pmMethod = new OCPMixtureMethodK_OGW01(rs_param, i, vs);
    }
    else if (rs_param.PVTW_T.data.size() > 0 &&
        (rs_param.PVDO_T.data.size() > 0 || rs_param.PVCDO_T.data.size() > 0)) {
        pmMethod = new OCPMixtureMethodK_OW01(rs_param, i, vs);
    }
    else if (rs_param.PVTCO2.data.size() > 0 && rs_param.PVTH2O.data.size() > 0) {
        pmMethod = new OCPMixtureMethodK_GW01(rs_param, i, vs);
    }
    else if (rs_param.PVDG_T.data.size() > 0 && rs_param.PVTW_T.data.size() > 0) {
        pmMethod = new OCPMixtureMethodK_GW02(rs_param, i, vs);
    }
    else if (rs_param.thermal) {
        pmMethod = new OCPMixtureMethodK_OW01T(rs_param.comsParam, i, vs);
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/02/2023      Create file                          */
/*----------------------------------------------------------------------------*/