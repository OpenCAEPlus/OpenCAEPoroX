/*! \file    OCPScalePcow.cpp
 *  \brief   OCPScalePcow class declaration
 *  \author  Shizhe Li
 *  \date    Jul/01/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPScalePcow.hpp"

/////////////////////////////////////////////////////////////////////
// Scale The Water-Oil Capillary Pressure Curves (From SWATINIT)
/////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////
// ScalePcowMethod01
////////////////////////////////////////////////////////////////


ScalePcowMethod01::ScalePcowMethod01(OCPFlow* flowin, ScalePcowVarSet& spvs) 
{

    switch (flowin->FlowType())
    {
    case OCPFlowType::OGW:
    case OCPFlowType::OW:
        break;
    default:
        OCP_ABORT("Wrong FlowType for ScalePcow!");
        break;
    }

    flow    = flowin;
    Swco    = flow->GetSwco();
    maxPcow = flow->GetMaxPcow();
    minPcow = flow->GetMinPcow();

    spvs.scaleVal.resize(spvs.swatInit.size());
}


void ScalePcowMethod01::SetScaleVal(const OCP_USI& bId, ScalePcowVarSet& spvs, OCP_DBL& Swinout, const OCP_DBL& Pcowin) const
{
    const auto& swatInit = spvs.swatInit[bId];
    auto&       sv       = spvs.scaleVal[bId];

    if (swatInit <= Swco) {
        Swinout = Swco;
        const OCP_DBL PcowInit = flow->CalPcowBySw(Swinout);
        sv = (Pcowin / PcowInit * maxPcow - minPcow) / (maxPcow - minPcow);
    }
    else {
        Swinout = swatInit;
        if (Pcowin > 0) {
            const OCP_DBL PcowInit = flow->CalPcowBySw(Swinout);
            if (PcowInit > 0) {
                sv = (Pcowin / PcowInit * maxPcow - minPcow) / (maxPcow - minPcow);
            }
        }
        else {
            sv = 1.0;
        }
    }   
}


void ScalePcowMethod01::ScaleDer(const OCP_USI& bId, const ScalePcowVarSet& spvs) const
{
    const auto&    sv = spvs.scaleVal[bId];
    OCPFlowVarSet& vs = flow->GetVarSet();
    const INT&     w  = vs.w;
    vs.Pc[w]          = -((-vs.Pc[w] - minPcow) * sv + minPcow);
    vs.dPcdS[vs.ww]   *= sv;
}


void ScalePcowMethod01::Scale(const OCP_USI& bId, const ScalePcowVarSet& spvs) const
{
    const auto&    sv = spvs.scaleVal[bId];
    OCPFlowVarSet& vs = flow->GetVarSet();
    const INT&    w  = vs.w;
    vs.Pc[w] = -((-vs.Pc[w] - minPcow) * sv + minPcow);
}


////////////////////////////////////////////////////////////////
// ScalePcow
////////////////////////////////////////////////////////////////


USI ScalePcow::Setup(OCPFlow* flowin)
{

    if (vs.swatInit.size() > 0) {
        ifScale = OCP_TRUE;

        scalePcowMethod.push_back(new ScalePcowMethod01(flowin, vs));
        return scalePcowMethod.size() - 1;
    }
    else {
        ifScale = OCP_FALSE;
    }
    return 0;
}


void ScalePcow::SetScaleVal(const OCP_USI& bId, const USI& mIndex, OCP_DBL& Swinout, const OCP_DBL& Pcowin)
{
    if (ifScale) {
        scalePcowMethod[mIndex]->SetScaleVal(bId, vs, Swinout, Pcowin);
    }  
}


void ScalePcow::ScaleDer(const OCP_USI& bId, const USI& mIndex) const
{
    if (ifScale) {
        scalePcowMethod[mIndex]->ScaleDer(bId, vs);
    }
}

void ScalePcow::Scale(const OCP_USI& bId, const USI& mIndex) const
{
    if (ifScale) {
        scalePcowMethod[mIndex]->Scale(bId, vs);
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/01/2023      Create file                          */
/*----------------------------------------------------------------------------*/