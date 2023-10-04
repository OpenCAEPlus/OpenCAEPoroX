/*! \file    OCPFlowGW.cpp
 *  \brief   OCPFlowGW class declaration
 *  \author  Shizhe Li
 *  \date    Sep/30/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


#include "OCPFlowGW.hpp"


/////////////////////////////////////////////////////
// OCPGWFMethod01
/////////////////////////////////////////////////////


void OCPGWFMethod01::CalKrPc(OCPFlowVarSet& vs)
{
	const INT& g = vs.g;
	const INT& w = vs.w;

	bc.CalKrPcN(vs.S[g], vs.kr[g], vs.Pc[g]);
	bc.CalKrPcW(vs.S[w], vs.kr[w], vs.Pc[w]);
}


void OCPGWFMethod01::CalKrPcDer(OCPFlowVarSet& vs)
{
	const INT& g = vs.g;
	const INT& w = vs.w;

	bc.CalKrPcDerN(vs.S[g], vs.kr[g], vs.Pc[g], vs.dKrdS[vs.gg], vs.dPcdS[vs.gg]);
	bc.CalKrPcDerW(vs.S[w], vs.kr[w], vs.Pc[w], vs.dKrdS[vs.ww], vs.dPcdS[vs.ww]);
}


/////////////////////////////////////////////////////
// OCPFlowGW
/////////////////////////////////////////////////////


void OCPFlowGW::Setup(const ParamReservoir& rs_param, const USI& i)
{
	if (rs_param.BCparam.size() > 0) {
		pfMethod = new OCPGWFMethod01(rs_param.BCparam[i], vs);
	}
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/30/2023      Create file                          */
/*----------------------------------------------------------------------------*/