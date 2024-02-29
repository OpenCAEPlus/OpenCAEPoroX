/*! \file    OCPDiffusion.cpp
 *  \brief   OCPDiffusion class declaration
 *  \author  Shizhe Li
 *  \date    Feb/28/2024
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPDiffusion.hpp"


OCPDiffusionMethod01::OCPDiffusionMethod01(const ParamReservoir& rs_param, OCPDiffusionVarSet& dvs)
{
    // temp
    diffuCP.resize(2);
    diffuCP[dvs.g].resize(dvs.nc, 1.6E-1);
    diffuCP[dvs.w].resize(dvs.nc, 1.0E-5);
}


void OCPDiffusionMethod01::CalFlux(const BulkConnPair& bp, const OCPDiffusionVarSet& dvs, const BulkVarSet& bvs, FluxVarSet& fvs) const
{

}


void OCPDiffusionMethod01::AssembleFIM(const BulkConnPair& bp, const OCPDiffusionVarSet& dvs, const BulkVarSet& bvs, FluxVarSet& fvs) const
{

}


void OCPDiffusion::Setup(const ParamReservoir& rs_param, const BulkVarSet& bvs)
{
    ifUse = OCP_FALSE;
    if (ifUse) {
        if (rs_param.thermal == OCP_FALSE) {
            ifUse = OCP_FALSE;
            OCP_WARNING("HEATCONDUCT is IGNORED in ISOTHERMAL MODEL!");
            return;
        }
        vs.SetUp(bvs.np, bvs.nc, bvs.o, bvs.g, bvs.w);
        dM.push_back(new OCPDiffusionMethod01(rs_param, vs));
    }
}

void OCPDiffusion::CalFlux(const BulkConnPair& bp, const BulkVarSet& bvs, FluxVarSet& fvs)
{
    if (ifUse) {
        dM[0]->CalFlux(bp, vs, bvs, fvs);
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Feb/28/2024      Create file                          */
/*----------------------------------------------------------------------------*/