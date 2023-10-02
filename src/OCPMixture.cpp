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



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/02/2023      Create file                          */
/*----------------------------------------------------------------------------*/