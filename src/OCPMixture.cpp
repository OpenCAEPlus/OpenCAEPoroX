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


OCPMixtureComp::OCPMixtureComp(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts)
{
    if (rs_param.PVTW_T.data[i].size() > 0) {
        pmMethod = new OCPMixtureCompMethod01(rs_param, i, vs);
    }

    // Skip stability analysis
    skipPSA = &opts.skipPSA;
    skipMethodIndex = skipPSA->Setup(opts.nb, pmMethod);
}



////////////////////////////////////////////////////////////////
// OCPMixtureK 
////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////
// OCPMixtureBlkOilOGW 
/////////////////////////////////////////////////////



OCPMixtureBlkOilOGW::OCPMixtureBlkOilOGW(const ParamReservoir& rs_param, const USI& i)
{
    if (rs_param.PVCO_T.data.size() > 0 &&
        rs_param.PVDG_T.data.size() > 0 &&
        rs_param.PVTW_T.data.size() > 0) {
        pmMethod = new OCPMixtureKOGWMethod01(rs_param, i, vs);
    }
}


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOW 
/////////////////////////////////////////////////////


OCPMixtureBlkOilOW::OCPMixtureBlkOilOW(const ParamReservoir& rs_param, const USI& i)
{
    if (rs_param.PVTW_T.data.size() > 0 &&
        (rs_param.PVDO_T.data.size() > 0 || rs_param.PVCDO_T.data.size() > 0)) {
        pmMethod = new OCPMixtureKOWMethod01(rs_param, i, vs);
    }
}


/////////////////////////////////////////////////////
// OCPMixtureBlkOilGW 
/////////////////////////////////////////////////////


OCPMixtureBlkOilGW::OCPMixtureBlkOilGW(const ParamReservoir& rs_param, const USI& i)
{
    if (rs_param.PVTCO2.data.size() > 0 && rs_param.PVTH2O.data.size() > 0) {
        pmMethod = new OCPMixtureKGWMethod01(rs_param, i, vs);
    }
}


/////////////////////////////////////////////////////
// OCPMixtureUnitThermalOW 
/////////////////////////////////////////////////////


OCPMixtureUnitThermalOW::OCPMixtureUnitThermalOW(const ParamReservoir& rs_param, const USI& i)
{
    pmMethod = new OCPMixtureKOWMethod01T(rs_param.comsParam, i, vs);
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/02/2023      Create file                          */
/*----------------------------------------------------------------------------*/