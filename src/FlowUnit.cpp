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
// FlowUnit_OW
///////////////////////////////////////////////


void FlowUnit_OW::SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin)
{
    scalePcow->SetScaleVal(bId, spMethodIndex, Swinout, Pcowin);
}


void FlowUnit_OW::CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId)
{
    bulkId = bId;
    OWF.CalKrPc(S_in[oIndex], S_in[wIndex]);
    scalePcow->Scale(bId, spMethodIndex);
    AssinValue();
}


void FlowUnit_OW::CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId)
{
    bulkId = bId;
    OWF.CalKrPcDer(S_in[oIndex], S_in[wIndex]);
    scalePcow->ScaleDer(bId, spMethodIndex);
    AssinValueDer();
}


void FlowUnit_OW::AssinValueDer()
{
    const OCPFlowVarSet& vs = OWF.GetVarSet();

    kr    = vs.kr;
    pc    = vs.Pc;
    dKrdS = vs.dKrdS;
    dPcdS = vs.dPcdS;
}


void FlowUnit_OW::AssinValue()
{
    const OCPFlowVarSet& vs = OWF.GetVarSet();

    kr = vs.kr;
    pc = vs.Pc;
}


///////////////////////////////////////////////
// FlowUnit_OG
///////////////////////////////////////////////


void FlowUnit_OG::CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId)
{
    bulkId = bId;
    OGF.CalKrPc(S_in[oIndex], S_in[gIndex]);
    AssinValue();
}

void FlowUnit_OG::CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId)
{
    bulkId = bId;
    OGF.CalKrPcDer(S_in[oIndex], S_in[gIndex]);
    AssinValueDer();
}


void FlowUnit_OG::AssinValueDer()
{
    const OCPFlowVarSet& vs = OGF.GetVarSet();

    kr    = vs.kr;
    pc    = vs.Pc;
    dKrdS = vs.dKrdS;
    dPcdS = vs.dPcdS;
}


void FlowUnit_OG::AssinValue()
{
    const OCPFlowVarSet& vs = OGF.GetVarSet();

    kr = vs.kr;
    pc = vs.Pc;
}


///////////////////////////////////////////////
// FlowUnit_GW
///////////////////////////////////////////////


void FlowUnit_GW::CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId)
{
    bulkId = bId;
    GWF.CalKrPc(S_in[gIndex], S_in[wIndex]);
    AssinValue();
}


void FlowUnit_GW::CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId)
{
    bulkId = bId;
    GWF.CalKrPcDer(S_in[gIndex], S_in[wIndex]);
    AssinValueDer();
}

void FlowUnit_GW::AssinValueDer()
{
    const OCPFlowVarSet& vs = GWF.GetVarSet();

    kr    = vs.kr;
    pc    = vs.Pc; 
    dKrdS = vs.dKrdS;
    dPcdS = vs.dPcdS;
}


void FlowUnit_GW::AssinValue()
{
    const OCPFlowVarSet& vs = GWF.GetVarSet();

    kr = vs.kr;
    pc = vs.Pc;
}

///////////////////////////////////////////////
// FlowUnit_OGW
///////////////////////////////////////////////


void FlowUnit_OGW::SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin)
{
    scalePcow->SetScaleVal(bId, spMethodIndex, Swinout, Pcowin);
}


void FlowUnit_OGW::CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId)
{
    bulkId = bId;
    OGWF.CalKrPc(S_in[oIndex], S_in[gIndex], S_in[wIndex]);
    misCurve->CorrectCurve(bId, mcMethodIndex);
    scalePcow->Scale(bId, spMethodIndex);
    AssinValue();
}


void FlowUnit_OGW::CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId)
{
    bulkId = bId;
    OGWF.CalKrPcDer(S_in[oIndex], S_in[gIndex], S_in[wIndex]);
    misCurve->CorrectCurveDer(bId, mcMethodIndex);
    scalePcow->ScaleDer(bId, spMethodIndex);
    AssinValueDer(); 
}


void FlowUnit_OGW::AssinValueDer()
{
    const OCPFlowVarSet& vs = OGWF.GetVarSet();
    kr    = vs.kr;
    pc    = vs.Pc;
    dKrdS = vs.dKrdS;
    dPcdS = vs.dPcdS;
}

void FlowUnit_OGW::AssinValue()
{
    const OCPFlowVarSet& vs = OGWF.GetVarSet();
    kr = vs.kr;
    pc = vs.Pc;
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/