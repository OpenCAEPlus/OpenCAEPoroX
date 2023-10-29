/*! \file    FlowUnit.cpp
 *  \brief   FlowUnit class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "FlowUnit.hpp"


///////////////////////////////////////////////
// FlowUnit
///////////////////////////////////////////////


FlowUnit::FlowUnit(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts)
{
    flow = new OCPFlow(rs_param, i);
    // for miscible
    misCurve      = &opts.misCur;
    mcMethodIndex = misCurve->Setup(flow, &opts.misFac);
    // for scalePcow
    scalePcow     = &opts.scalePcow;
    spMethodIndex = scalePcow->Setup(flow);
}


void FlowUnit::SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin) const
{
    scalePcow->SetScaleVal(bId, spMethodIndex, Swinout, Pcowin);
}


void FlowUnit::CalKrPc(const OCP_USI& bId, const OCP_DBL* S) const
{
    flow->CalKrPc(S);
    misCurve->CorrectCurve(bId, mcMethodIndex);
    scalePcow->Scale(bId, spMethodIndex);
}


void FlowUnit::CalKrPcFIM(const OCP_USI& bId, const OCP_DBL* S) const
{
    flow->CalKrPcDer(S);
    misCurve->CorrectCurveDer(bId, mcMethodIndex);
    scalePcow->ScaleDer(bId, spMethodIndex);
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/