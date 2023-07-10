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
// FlowUnit_W
///////////////////////////////////////////////

void FlowUnit_W::CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId)
{
    // constant value
    //kr_out[0] = 1;
    //pc_out[0] = 0;
}

void FlowUnit_W::CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId)
{
    // constant value
    //kr_out[0] = 1;
    //pc_out[0] = 0;
    //dkrdS[0]  = 0;
    //dPcjdS[0] = 0;
}

///////////////////////////////////////////////
// FlowUnit_OW
///////////////////////////////////////////////

FlowUnit_OW::FlowUnit_OW(const ParamReservoir& rs_param, const USI& i)
{
    Allocate(2);

    SWOF.Setup(rs_param.SWOF_T.data[i]);
    Swco = SWOF.GetSwco();

    data.resize(4, 0);
    cdata.resize(4, 0);
}

void FlowUnit_OW::CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId)
{
    const OCP_DBL Sw = S_in[1];

    // three phase black oil model using stone 2
    OCP_DBL krw, krow, Pcwo;
    SWOF.CalKrwKrowPcwo(Sw, krw, krow, Pcwo);

    kr[0] = krow;
    kr[1] = krw;
    //pc[0] = 0;
    pc[1] = Pcwo;
}

void FlowUnit_OW::CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId)
{
    const OCP_DBL Sw = S_in[1];
    OCP_DBL krw, krow, Pcwo, dKrwdSw, dKrowdSw, dPcwodSw;
    SWOF.CalKrwKrowPcwoDer(Sw, krw, krow, Pcwo, dKrwdSw, dKrowdSw, dPcwodSw);

    kr[0] = krow;
    kr[1] = krw;
    //pc[0] = 0;
    pc[1] = Pcwo;

    //dKrdS[0] = 0;
    dKrdS[1] = dKrowdSw;
    //dKrdS[2] = 0;
    dKrdS[3] = dKrwdSw;

    //dPcdS[0] = 0;
    //dPcdS[1] = 0;
    //dPcdS[2] = 0;
    dPcdS[3] = dPcwodSw;
}

///////////////////////////////////////////////
// FlowUnit_OG
///////////////////////////////////////////////

FlowUnit_OG::FlowUnit_OG(const ParamReservoir& rs_param, const USI& i)
{
    Allocate(2);

    SGOF.Setup(rs_param.SGOF_T.data[i]);

    data.resize(4, 0);
    cdata.resize(4, 0);
}

void FlowUnit_OG::CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId)
{
    const OCP_DBL Sg = S_in[1];

    OCP_DBL krg, krog, Pcgo;
    SGOF.CalKrgKrogPcgo(Sg, krg, krog, Pcgo);

    kr[0] = krog;
    kr[1] = krg;
    //pc[0] = 0;
    pc[1] = Pcgo;
}

void FlowUnit_OG::CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId)
{
    OCP_ABORT("Not Completed Now!");
}

///////////////////////////////////////////////
// FlowUnit_OGW
///////////////////////////////////////////////


void FlowUnit_OGW::SetupOptionalFeatures(OptionalFeatures& optFeatures)
{
    // for miscible
    miscible       = &optFeatures.miscible;
    mcMethodIndex  = miscible->Setup(&PF3);
    // for scalePcow
    scalePcow      = &optFeatures.scalePcow;
    spMethodIndex  = scalePcow->Setup(&PF3);
}

void FlowUnit_OGW::SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin)
{
    scalePcow->SetScaleVal(bId, spMethodIndex, Swinout, Pcowin);
}


void FlowUnit_OGW::CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId)
{
    bulkId = bId;
    PF3.CalKrPc(S_in[oIndex], S_in[gIndex], S_in[wIndex]);
    miscible->CorrectCurve(bId, mcMethodIndex);
    scalePcow->Scale(bId, spMethodIndex);
    AssinValue();
}


void FlowUnit_OGW::CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId)
{
    bulkId = bId;
    PF3.CalKrPcDer(S_in[oIndex], S_in[gIndex], S_in[wIndex]);
    miscible->CorrectCurve(bId, mcMethodIndex);
    scalePcow->ScaleDer(bId, spMethodIndex);
    AssinValueDer(); 
}


void FlowUnit_OGW::AssinValueDer()
{
    const OCP3PFVarSet& vs = PF3.GetVarSet();
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
    const OCP3PFVarSet& vs = PF3.GetVarSet();
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