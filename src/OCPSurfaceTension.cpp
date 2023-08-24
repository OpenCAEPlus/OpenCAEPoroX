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


void SurTenMethod01::CalSurfaceTension(const OCP_USI& bId, SurTenVarSet& stvs) const
{
    if (mixture->GetNPPE() == 1)
        stvs.surTen[bId] = OCP_HUGE;
    else {
        const OCPMixtureVarSet& vs = mixture->GetVarSet();
        const OCP_DBL           Bl = vs.xi[0] * CONV7;
        const OCP_DBL           Bv = vs.xi[1] * CONV7;    
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

USI SurfaceTension::Setup(const ParamReservoir& rs_param, const USI& i, const OCP_USI& nb, OCPMixtureComp* mix)
{

    if (rs_param.comsParam.parachor.data.size() > 0) {
        ifUse = OCP_TRUE;
        vs.SetNb(nb);
        stMethod.push_back(new SurTenMethod01(rs_param.comsParam.parachor.data[i], mix, vs));      
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