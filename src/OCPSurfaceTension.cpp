/*! \file    OCPSurfaceTension.cpp
 *  \brief   OCPSurfaceTension class declaration
 *  \author  Shizhe Li
 *  \date    Jul/03/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPSurfaceTension.hpp"


 /////////////////////////////////////////////////////////////////////
 // SurTenMethod01
 /////////////////////////////////////////////////////////////////////


OCP_DBL SurTenMethod01::CalSurfaceTension() const
{
    if (*NP == 1)
        return OCP_HUGE;
    else {
        const OCP_DBL Bv = *bv * CONV7;
        const OCP_DBL Bl = *bl * CONV7;
        OCP_DBL surTen = 0;
        for (USI i = 0; i < NC; i++)
            surTen += parachor[i] * (Bv * xv[i] - Bl * xl[i]);
        return pow(surTen, 4.0);
    }
}


/////////////////////////////////////////////////////////////////////
// SurfaceTension
/////////////////////////////////////////////////////////////////////

USI SurfaceTension::Setup(const SurTenMethodParams& params)
{
    if (params.params01.IfUse()) {
        stMethod.push_back(new SurTenMethod01(params.params01));
    }
    else {
        OCP_ABORT("NO Matched Surface Tension Model!");
    }

    return stMethod.size() - 1;
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/03/2022      Create file                          */
/*----------------------------------------------------------------------------*/