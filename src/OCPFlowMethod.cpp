/*! \file    OCPFlowMethod.cpp
 *  \brief   OCPFlowMethod class declaration
 *  \author  Shizhe Li
 *  \date    Oct/04/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPFlowMethod.hpp"

 /////////////////////////////////////////////////////
 // OCPFlowMethod_OGW01
 /////////////////////////////////////////////////////


OCPFlowMethod_OGW01::OCPFlowMethod_OGW01(const vector<vector<OCP_DBL>>& SGOFin,
	const vector<vector<OCP_DBL>>& SWOFin,
	const USI& i, OCPFlowVarSet& vs)
{
	vs.Init(OCPFlowType::OGW, 3, 3);

	SGOF.Setup(SGOFin);
	SWOF.Setup(SWOFin);
	vs.krocw = SWOF.GetKrocw();
	vs.Swco = SWOF.GetSwco();

	Generate_SWPCWG();

	opC.Setup(1);
}


void OCPFlowMethod_OGW01::Generate_SWPCWG()
{
	const std::vector<OCP_DBL> Sw(SWOF.GetSw());
	std::vector<OCP_DBL>       Pcow(SWOF.GetPcow());

	for (USI i = 0; i < Sw.size(); i++) {
		OCP_DBL Pcgo = SGOF.CalPcgo(1 - Sw[i]);
		Pcow[i] += Pcgo; // Pcgw
	}

	SWPCGW.Setup(vector<vector<OCP_DBL>>{Sw, Pcow});
}


void OCPFlowMethod_OGW01::CalKrPc(OCPFlowVarSet& vs)
{
	const INT& g = vs.g;
	const INT& w = vs.w;

	SWOF.CalKrwKrowPcwo(vs.S[w], vs.kr[w], vs.krow, vs.Pc[w]);

	SGOF.CalKrgKrogPcgo(vs.S[g], vs.kr[g], vs.krog, vs.Pc[g]);

	opC.CalOilPer(vs);
}


void OCPFlowMethod_OGW01::CalKrPcDer(OCPFlowVarSet& vs)
{
	const INT& g = vs.g;
	const INT& w = vs.w;

	SWOF.CalKrwKrowPcwoDer(vs.S[w], vs.kr[w], vs.krow, vs.Pc[w], vs.dKrdS[vs.ww], vs.dKrowdSw, vs.dPcdS[vs.ww]);

	SGOF.CalKrgKrogPcgoDer(vs.S[g], vs.kr[g], vs.krog, vs.Pc[g], vs.dKrdS[vs.gg], vs.dKrogdSg, vs.dPcdS[vs.gg]);

	opC.CalOilPerDer(vs);
}

/////////////////////////////////////////////////////
// OCPFlowMethod_OGW02
/////////////////////////////////////////////////////

OCPFlowMethod_OGW02::OCPFlowMethod_OGW02(const vector<vector<OCP_DBL>>& SOF3in,
	const vector<vector<OCP_DBL>>& SGFNin,
	const vector<vector<OCP_DBL>>& SWFNin,
	const USI& i, OCPFlowVarSet& vs)
{
	vs.Init(OCPFlowType::OGW, 3, 3);

	SOF3.Setup(SOF3in);
	SGFN.Setup(SGFNin);
	SWFN.Setup(SWFNin);

	vs.krocw = SOF3.GetKrocw();
	vs.Swco = SWFN.GetSwco();

	Generate_SWPCWG();

	opC.Setup(1);
}


void OCPFlowMethod_OGW02::CalKrPc(OCPFlowVarSet& vs)
{
	const INT& o = vs.o;
	const INT& g = vs.g;
	const INT& w = vs.w;

	SWFN.CalKrwPcwo(vs.S[w], vs.kr[w], vs.Pc[w]);

	SGFN.CalKrgPcgo(vs.S[g], vs.kr[g], vs.Pc[g]);

	SOF3.CalKrowKrog(vs.S[o], vs.krow, vs.krog);

	opC.CalOilPer(vs);
}


void OCPFlowMethod_OGW02::CalKrPcDer(OCPFlowVarSet& vs)
{
	const INT& o = vs.o;
	const INT& g = vs.g;
	const INT& w = vs.w;

	SWFN.CalKrwPcwoDer(vs.S[w], vs.kr[w], vs.Pc[w], vs.dKrdS[vs.ww], vs.dPcdS[vs.ww]);

	SGFN.CalKrgPcgoDer(vs.S[g], vs.kr[g], vs.Pc[g], vs.dKrdS[vs.gg], vs.dPcdS[vs.gg]);

	SOF3.CalKrowKrogDer(vs.S[o], vs.krow, vs.krog, vs.dKrowdSo, vs.dKrogdSo);

	opC.CalOilPerDer(vs);

}


void OCPFlowMethod_OGW02::Generate_SWPCWG()
{

	const std::vector<OCP_DBL> Sw(SWFN.GetSw());
	std::vector<OCP_DBL> Pcow(SWFN.GetPcow());
	USI                  n = Sw.size();
	for (USI i = 0; i < n; i++) {
		// OCP_DBL Pcgo = SGFN.Eval(0, 1 - Sw[i], 2);
		OCP_DBL Pcgo = SGFN.CalPcgo(1 - Sw[i]);
		Pcow[i] += Pcgo; // pcgw
	}

	SWPCGW.Setup(vector<vector<OCP_DBL>>{Sw, Pcow});
}


/////////////////////////////////////////////////////
// OCPFlowMethod_OW01
/////////////////////////////////////////////////////

OCPFlowMethod_OW01::OCPFlowMethod_OW01(const vector<vector<OCP_DBL>>& SWOFin, OCPFlowVarSet& vs)
{
	vs.Init(OCPFlowType::OW, 2, 2);

	SWOF.Setup(SWOFin);

	vs.Swco = SWOF.GetSwco();
}


void OCPFlowMethod_OW01::CalKrPc(OCPFlowVarSet& vs)
{
	const INT& o = vs.o;
	const INT& w = vs.w;
	SWOF.CalKrwKrowPcwo(vs.S[w], vs.kr[w], vs.kr[o], vs.Pc[w]);
}


void OCPFlowMethod_OW01::CalKrPcDer(OCPFlowVarSet& vs)
{
	const INT& o = vs.o;
	const INT& w = vs.w;
	SWOF.CalKrwKrowPcwoDer(vs.S[w], vs.kr[w], vs.kr[o], vs.Pc[w], vs.dKrdS[vs.ww], vs.dKrdS[vs.ow], vs.dPcdS[vs.ww]);
}


/////////////////////////////////////////////////////
// OCPFlowMethod_OG01
/////////////////////////////////////////////////////

OCPFlowMethod_OG01::OCPFlowMethod_OG01(const vector<vector<OCP_DBL>>& SWOFin, OCPFlowVarSet& vs)
{
	vs.Init(OCPFlowType::OG, 2, 2);
	SGOF.Setup(SWOFin);
}


void OCPFlowMethod_OG01::CalKrPc(OCPFlowVarSet& vs)
{
	const INT& o = vs.o;
	const INT& g = vs.g;
	SGOF.CalKrgKrogPcgo(vs.S[g], vs.kr[g], vs.kr[o], vs.Pc[g]);
}


void OCPFlowMethod_OG01::CalKrPcDer(OCPFlowVarSet& vs)
{
	const INT& o = vs.o;
	const INT& g = vs.g;
	SGOF.CalKrgKrogPcgoDer(vs.S[g], vs.kr[g], vs.kr[o], vs.Pc[g], vs.dKrdS[vs.gg], vs.dKrdS[vs.og], vs.dPcdS[vs.gg]);
}


/////////////////////////////////////////////////////
// OCPFlowMethod_GW01
/////////////////////////////////////////////////////


void OCPFlowMethod_GW01::CalKrPc(OCPFlowVarSet& vs)
{
	const INT& g = vs.g;
	const INT& w = vs.w;

	// wetting phase is water
	bc.CalKrPcN(vs.S[g], vs.kr[g], vs.Pc[w]);
	bc.CalKrPcW(vs.S[w], vs.kr[w], vs.Pc[g]);
}


void OCPFlowMethod_GW01::CalKrPcDer(OCPFlowVarSet& vs)
{
	const INT& g = vs.g;
	const INT& w = vs.w;

	// wetting phase is water
	bc.CalKrPcDerN(vs.S[g], vs.kr[g], vs.Pc[w], vs.dKrdS[vs.gg], vs.dPcdS[vs.wg]);
	bc.CalKrPcDerW(vs.S[w], vs.kr[w], vs.Pc[g], vs.dKrdS[vs.ww], vs.dPcdS[vs.gw]);
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/04/2023      Create file                          */
/*----------------------------------------------------------------------------*/