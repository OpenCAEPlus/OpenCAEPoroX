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


void FlowUnit_OW::SetupOptionalFeatures(OptionalFeatures& optFeatures)
{
    // for scalePcow
    scalePcow     = &optFeatures.scalePcow;
    spMethodIndex = scalePcow->Setup(&OWF);
}


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

    kr[oIndex]  = vs.kro;
    kr[wIndex]  = vs.krw;
    pc[oIndex]  = 0;
    pc[wIndex]  = vs.Pcwo;

    dKrdS[oo] = vs.dKrodSo;
    dKrdS[ow] = vs.dKrodSw;
    dKrdS[wo] = vs.dKrwdSo;
    dKrdS[ww] = vs.dKrwdSw;

    dPcdS[oo] = 0;
    dPcdS[ow] = 0;
    dPcdS[wo] = vs.dPcwodSo;
    dPcdS[ww] = vs.dPcwodSw;
}


void FlowUnit_OW::AssinValue()
{
    const OCPFlowVarSet& vs = OWF.GetVarSet();

    kr[oIndex] = vs.kro;
    kr[wIndex] = vs.krw;
    pc[oIndex] = 0;
    pc[wIndex] = vs.Pcwo;
}


///////////////////////////////////////////////
// FlowUnit_OG
///////////////////////////////////////////////


FlowUnit_OG::FlowUnit_OG(const ParamReservoir& rs_param, const USI& i)
{
    Allocate(2);
    OGF.Setup(rs_param, i);
}

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

    kr[oIndex] = vs.kro;
    kr[gIndex] = vs.krg;
    pc[oIndex] = 0;
    pc[gIndex] = vs.Pcgo;

    dKrdS[oo] = vs.dKrodSo;
    dKrdS[og] = vs.dKrodSg;
    dKrdS[go] = vs.dKrgdSo;
    dKrdS[gg] = vs.dKrgdSg;

    dPcdS[oo] = 0;
    dPcdS[og] = 0;
    dPcdS[go] = vs.dPcgodSo;
    dPcdS[gg] = vs.dPcgodSg;
}


void FlowUnit_OG::AssinValue()
{
    const OCPFlowVarSet& vs = OGF.GetVarSet();

    kr[oIndex] = vs.kro;
    kr[gIndex] = vs.krg;
    pc[oIndex] = 0;
    pc[gIndex] = vs.Pcgo;
}

///////////////////////////////////////////////
// FlowUnit_OGW
///////////////////////////////////////////////


void FlowUnit_OGW::SetupOptionalFeatures(OptionalFeatures& optFeatures)
{
    // for miscible
    miscible       = &optFeatures.miscible;
    mcMethodIndex  = miscible->Setup(&OGWF);
    // for scalePcow
    scalePcow      = &optFeatures.scalePcow;
    spMethodIndex  = scalePcow->Setup(&OGWF);
}

void FlowUnit_OGW::SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin)
{
    scalePcow->SetScaleVal(bId, spMethodIndex, Swinout, Pcowin);
}


void FlowUnit_OGW::CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId)
{
    bulkId = bId;
    OGWF.CalKrPc(S_in[oIndex], S_in[gIndex], S_in[wIndex]);
    miscible->CorrectCurve(bId, mcMethodIndex);
    scalePcow->Scale(bId, spMethodIndex);
    AssinValue();
}


void FlowUnit_OGW::CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId)
{
    bulkId = bId;
    OGWF.CalKrPcDer(S_in[oIndex], S_in[gIndex], S_in[wIndex]);
    miscible->CorrectCurveDer(bId, mcMethodIndex);
    scalePcow->ScaleDer(bId, spMethodIndex);
    AssinValueDer(); 
}


void FlowUnit_OGW::AssinValueDer()
{
    const OCPFlowVarSet& vs = OGWF.GetVarSet();
    kr[oIndex] = vs.kro;
    kr[gIndex] = vs.krg;
    kr[wIndex] = vs.krw;
    pc[oIndex] = 0;;
    pc[gIndex] = vs.Pcgo;
    pc[wIndex] = vs.Pcwo;

    dKrdS[oo] = vs.dKrodSo;
    dKrdS[og] = vs.dKrodSg;
    dKrdS[ow] = vs.dKrodSw;
    dKrdS[go] = vs.dKrgdSo;
    dKrdS[gg] = vs.dKrgdSg;
    dKrdS[gw] = vs.dKrgdSw;
    dKrdS[wo] = vs.dKrwdSo;
    dKrdS[wg] = vs.dKrwdSg;
    dKrdS[ww] = vs.dKrwdSw;

    dPcdS[oo] = 0;
    dPcdS[og] = 0;
    dPcdS[ow] = 0;
    dPcdS[go] = vs.dPcgodSo;
    dPcdS[gg] = vs.dPcgodSg;
    dPcdS[gw] = vs.dPcgodSw;
    dPcdS[wo] = vs.dPcwodSo;
    dPcdS[wg] = vs.dPcwodSg;
    dPcdS[ww] = vs.dPcwodSw;
}

void FlowUnit_OGW::AssinValue()
{
    const OCPFlowVarSet& vs = OGWF.GetVarSet();
    kr[oIndex] = vs.kro;
    kr[gIndex] = vs.krg;
    kr[wIndex] = vs.krw;
    pc[oIndex] = 0;
    pc[gIndex] = vs.Pcgo;
    pc[wIndex] = vs.Pcwo;
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/