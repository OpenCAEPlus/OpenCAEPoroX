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


ScalePcowMethod01::ScalePcowMethod01(OCP3PhaseFlow* pf3) {
    PF3     = pf3;
    Swco    = PF3->GetSwco();
    maxPcow = PF3->GetMaxPcow();
    minPcow = PF3->GetMinPcow();
}


OCP_DBL ScalePcowMethod01::SetScaleVal(const OCP_DBL& swatInit, OCP_DBL& Swinout, const OCP_DBL& Pcowin) const
{
    if (swatInit <= Swco) {
        Swinout = Swco;
        const OCP_DBL PcowInit = PF3->CalPcowBySw(Swinout);
        return (Pcowin / PcowInit * maxPcow - minPcow) / (maxPcow - minPcow);
    }
    else {
        Swinout = swatInit;
        if (Pcowin > 0) {
            const OCP_DBL PcowInit = PF3->CalPcowBySw(Swinout);
            if (PcowInit > 0) {
                return (Pcowin / PcowInit * maxPcow - minPcow) / (maxPcow - minPcow);
            }
        }
    }
    return 1.0;
}


void ScalePcowMethod01::ScaleDer(const OCP_DBL& sv) const
{
    OCP3PFVarSet vs = PF3->GetVarSet();
    vs.Pcwo      = -((-vs.Pcwo - minPcow) * sv + minPcow);
    vs.dPcwodSw *= sv;
}


void ScalePcowMethod01::Scale(const OCP_DBL& sv) const
{
    OCP3PFVarSet vs = PF3->GetVarSet();
    vs.Pcwo = -((-vs.Pcwo - minPcow) * sv + minPcow);
}


////////////////////////////////////////////////////////////////
// ScalePcow
////////////////////////////////////////////////////////////////


USI ScalePcow::Setup(OCP3PhaseFlow* pf3)
{
    if (ifScale) {
        OCP_ASSERT(swatInit.size() > 0, "SWATINIT is MISSING!");
        scaleVal.resize(swatInit.size(), 1.0);

        scalePcowMethod.push_back(new ScalePcowMethod01(pf3));
        return scalePcowMethod.size() - 1;
    }
    else {
        return 0;
    }
}


void ScalePcow::SetScaleVal(const OCP_USI& bId, const USI& mIndex, OCP_DBL& Swinout, const OCP_DBL& Pcowin)
{
    if (ifScale) {
        scaleVal[bId] = scalePcowMethod[mIndex]->SetScaleVal(swatInit[bId], Swinout, Pcowin);
    }  
}


void ScalePcow::ScaleDer(const OCP_USI& bId, const USI& mIndex) const
{
    if (ifScale) {
        scalePcowMethod[mIndex]->ScaleDer(scaleVal[bId]);
    }
}

void ScalePcow::Scale(const OCP_USI& bId, const USI& mIndex) const
{
    if (ifScale) {
        scalePcowMethod[mIndex]->Scale(scaleVal[bId]);
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/01/2023      Create file                          */
/*----------------------------------------------------------------------------*/