/*! \file    OCPFlowOG.cpp
 *  \brief   OCPFlowOG class declaration
 *  \author  Shizhe Li
 *  \date    Jul/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPFlowOG.hpp"


 /////////////////////////////////////////////////////
 // OCPOGFMethod01
 /////////////////////////////////////////////////////

OCPOGFMethod01::OCPOGFMethod01(const vector<vector<OCP_DBL>>& SWOFin, OCPFlowVarSet& vs)
{
    vs.Init(OCPFlowType::OG, 2, 2);
    SGOF.Setup(SWOFin);
}


void OCPOGFMethod01::CalKrPc(OCPFlowVarSet& vs)
{
    const INT& o = vs.o;
    const INT& g = vs.g;
    SGOF.CalKrgKrogPcgo(vs.S[g], vs.kr[g], vs.kr[o], vs.Pcg);
}


void OCPOGFMethod01::CalKrPcDer(OCPFlowVarSet& vs)
{
    const INT& o = vs.o;
    const INT& g = vs.g;
    SGOF.CalKrgKrogPcgoDer(vs.S[g], vs.kr[g], vs.kr[o], vs.Pcg, vs.dKrgdSg, vs.dKrodSg, vs.dPcgdSg);
}


void OCPFlowOG::Setup(const ParamReservoir& rs_param, const USI& i)
{
    if (rs_param.SGOF_T.data.size() > 0) {
        pfMethod = new OCPOGFMethod01(rs_param.SGOF_T.data[i], vs);
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/