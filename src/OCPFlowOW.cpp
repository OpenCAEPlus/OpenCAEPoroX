/*! \file    OCPFlowOW.cpp
 *  \brief   OCPFlowOW class declaration
 *  \author  Shizhe Li
 *  \date    Jul/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPFlowOW.hpp"


/////////////////////////////////////////////////////
// OCPOWFMethod01
/////////////////////////////////////////////////////

OCPOWFMethod01::OCPOWFMethod01(const vector<vector<OCP_DBL>>& SWOFin, OCPFlowVarSet& vs)
{
    vs.Init(OCPFlowType::OW, 2, 2);

    SWOF.Setup(SWOFin);

    vs.Swco = SWOF.GetSwco();
}


void OCPOWFMethod01::CalKrPc(OCPFlowVarSet& vs)
{
    const INT& o = vs.o;
    const INT& w = vs.w;
    SWOF.CalKrwKrowPcwo(vs.S[w], vs.kr[w], vs.kr[o], vs.Pc[w]);
}


void OCPOWFMethod01::CalKrPcDer(OCPFlowVarSet& vs)
{
    const INT& o = vs.o;
    const INT& w = vs.w;
    SWOF.CalKrwKrowPcwoDer(vs.S[w], vs.kr[w], vs.kr[o], vs.Pc[w], vs.dKrdS[vs.ww], vs.dKrdS[vs.ow], vs.dPcwdSw);
}



void OCPFlowOW::Setup(const ParamReservoir& rs_param, const USI& i)
{
    if (rs_param.SWOF_T.data.size() > 0) {
        pfMethod = new OCPOWFMethod01(rs_param.SWOF_T.data[i], vs);
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/