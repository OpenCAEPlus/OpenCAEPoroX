/*! \file    OCPSurfaceTension.cpp
 *  \brief   OCPSurfaceTension class declaration
 *  \author  Shizhe Li
 *  \date    Jul/03/2023
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


void SurTenMethod01::CalSurfaceTension(const OCP_USI& bId, SurTenVarSet& stvs, const OCPMixtureVarSet& mvs) const
{
    if (!mvs.phaseExist[0] || !mvs.phaseExist[1])
        stvs.surTen[bId] = OCP_HUGE;
    else {
        const OCPMixtureVarSet& vs = mvs;
        const OCP_DBL           Bl = vs.xi[0] * CONV6;
        const OCP_DBL           Bv = vs.xi[1] * CONV6;    
        const OCP_DBL*          xl = vs.GetXj(0);
        const OCP_DBL*          xv = vs.GetXj(1);
        OCP_DBL surTen = 0;
        for (USI i = 0; i < NC; i++)
            surTen += parachor[i] * (Bv * xv[i] - Bl * xl[i]);
        stvs.surTen[bId] = pow(surTen, 4.0);
    }
}


/////////////////////////////////////////////////////////////////////
// SurfaceTension
/////////////////////////////////////////////////////////////////////

USI SurfaceTension::Setup(const ParamReservoir& rs_param, const USI& i, const OCP_USI& nb)
{

    if (rs_param.comsParam.parachor.data.size() > 0) {
        ifUse = OCP_TRUE;
        vs.SetNb(nb);
        stMethod.push_back(new SurTenMethod01(rs_param.comsParam.parachor.data[i], vs));      
    }

    return stMethod.size() - 1;
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/03/2023      Create file                          */
/*----------------------------------------------------------------------------*/