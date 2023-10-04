/*! \file    OCPFlowOGW.cpp
 *  \brief   OCPFlowOGW class declaration
 *  \author  Shizhe Li
 *  \date    Jul/08/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPFlowOGW.hpp"


 /////////////////////////////////////////////////////
 // OCP3POilPerMethod01
 /////////////////////////////////////////////////////


void OCP3POilPerMethod01::CalOilPer(OCPFlowVarSet& vs)
{
	const INT& o = vs.o;
	const INT& g = vs.g;
	const INT& w = vs.w;

	vs.kr[o] = vs.krocw 
		      * ((vs.krow / vs.krocw + vs.kr[w]) * (vs.krog / vs.krocw + vs.kr[g]) 
		      - (vs.kr[w] + vs.kr[g]));
	if (vs.kr[o] < 0) vs.kr[o] = 0;
}


void OCP3POilPerMethod01::CalOilPerDer(OCPFlowVarSet& vs)
{
	const INT& o = vs.o;
	const INT& g = vs.g;
	const INT& w = vs.w;

	vs.kr[o] = vs.krocw
		* ((vs.krow / vs.krocw + vs.kr[w]) * (vs.krog / vs.krocw + vs.kr[g])
			- (vs.kr[w] + vs.kr[g]));

	if (vs.kr[o] < 0) {
		vs.kr[o]     = 0;
		vs.dKrodSo = 0;
		vs.dKrodSg = 0;
		vs.dKrodSw = 0;
	}
	else {
		vs.dKrodSo = vs.dKrowdSo * (vs.krog / vs.krocw + vs.kr[g]) 
			          + vs.dKrogdSo * (vs.krow / vs.krocw + vs.kr[w]);
		vs.dKrodSw = vs.krocw * ((vs.dKrowdSw / vs.krocw + vs.dKrwdSw) 
			          * (vs.krog / vs.krocw + vs.kr[g]) - (vs.dKrwdSw));
		vs.dKrodSg = vs.krocw * ((vs.krow / vs.krocw + vs.kr[w])
			          * (vs.dKrogdSg / vs.krocw + vs.dKrgdSg) - (vs.dKrgdSg));
	}
}

/////////////////////////////////////////////////////
// OCP3POilPerMethod02
/////////////////////////////////////////////////////


void OCP3POILPerMethod02::CalOilPer(OCPFlowVarSet& vs)
{
	const INT& o = vs.o;
	const INT& g = vs.g;
	const INT& w = vs.w;

	const OCP_DBL tmp = vs.S[g] + vs.S[w] - vs.Swco;
	if (tmp <= TINY) {
		vs.kr[o] = vs.krocw;
	}
	else {
		vs.kr[o] = (vs.S[g] * vs.krog + (vs.S[w] - vs.Swco) * vs.krow) / tmp;
	}	
}


void OCP3POILPerMethod02::CalOilPerDer(OCPFlowVarSet& vs)
{
	const INT& o = vs.o;
	const INT& g = vs.g;
	const INT& w = vs.w;

	const OCP_DBL tmp = vs.S[g] + vs.S[w] - vs.Swco;
	if (tmp <= TINY) {
		vs.kr[o]   = vs.krocw;
		vs.dKrodSo = 0;
		vs.dKrodSg = 0;
		vs.dKrodSw = 0;
	}
	else {
		vs.kr[o] = (vs.S[g] * vs.krog + (vs.S[w] - vs.Swco) * vs.krow) / tmp;
		vs.dKrodSo = (vs.S[g] * vs.dKrogdSo + (vs.S[w] - vs.Swco) * vs.dKrowdSo) / tmp;
		vs.dKrodSg = (vs.krog + vs.S[g] * vs.dKrogdSg - vs.kr[o]) / tmp;
		vs.dKrodSw = (vs.krow + (vs.S[w] - vs.Swco) * vs.dKrowdSw - vs.kr[o]) / tmp;
	}
}


/////////////////////////////////////////////////////
// OCPOGWFMethod01
/////////////////////////////////////////////////////


OCPOGWFMethod01::OCPOGWFMethod01(const vector<vector<OCP_DBL>>& SGOFin,
	const vector<vector<OCP_DBL>>& SWOFin,
	const USI& i, OCPFlowVarSet& vs)
{
	vs.Init(OCPFlowType::OGW, 3, 3);

	SGOF.Setup(SGOFin);
	SWOF.Setup(SWOFin);
	vs.krocw = SWOF.GetKrocw();
	vs.Swco  = SWOF.GetSwco();

	Generate_SWPCWG();

	opC.Setup(1);
}


void OCPOGWFMethod01::Generate_SWPCWG()
{
	const std::vector<OCP_DBL> Sw(SWOF.GetSw());
	std::vector<OCP_DBL>       Pcow(SWOF.GetPcow());

	for (USI i = 0; i < Sw.size(); i++) {
		OCP_DBL Pcgo = SGOF.CalPcgo(1 - Sw[i]);
		Pcow[i] += Pcgo; // Pcgw
	}

	SWPCGW.Setup(vector<vector<OCP_DBL>>{Sw, Pcow});
}


void OCPOGWFMethod01::CalKrPc(OCPFlowVarSet& vs)
{
	const INT& g = vs.g;
	const INT& w = vs.w;

	SWOF.CalKrwKrowPcwo(vs.S[w], vs.kr[w], vs.krow, vs.Pc[w]);

	SGOF.CalKrgKrogPcgo(vs.S[g], vs.kr[g], vs.krog, vs.Pc[g]);

	opC.CalOilPer(vs);
}


void OCPOGWFMethod01::CalKrPcDer(OCPFlowVarSet& vs)
{
	const INT& g = vs.g;
	const INT& w = vs.w;

	SWOF.CalKrwKrowPcwoDer(vs.S[w], vs.kr[w], vs.krow, vs.Pc[w], vs.dKrwdSw, vs.dKrowdSw, vs.dPcwdSw);

	SGOF.CalKrgKrogPcgoDer(vs.S[g], vs.kr[g], vs.krog, vs.Pc[g], vs.dKrgdSg, vs.dKrogdSg, vs.dPcgdSg);

	opC.CalOilPerDer(vs);
}

/////////////////////////////////////////////////////
// OCPOGWFMethod01
/////////////////////////////////////////////////////

OCPOGWFMethod02::OCPOGWFMethod02(const vector<vector<OCP_DBL>>& SOF3in,
	const vector<vector<OCP_DBL>>& SGFNin,
	const vector<vector<OCP_DBL>>& SWFNin,
	const USI& i, OCPFlowVarSet& vs)
{
	vs.Init(OCPFlowType::OGW, 3, 3);

	SOF3.Setup(SOF3in);	
	SGFN.Setup(SGFNin);
	SWFN.Setup(SWFNin);

	vs.krocw = SOF3.GetKrocw();
	vs.Swco  = SWFN.GetSwco();

	Generate_SWPCWG();

	opC.Setup(1);
}


void OCPOGWFMethod02::CalKrPc(OCPFlowVarSet& vs)
{
	const INT& o = vs.o;
	const INT& g = vs.g;
	const INT& w = vs.w;

	SWFN.CalKrwPcwo(vs.S[w], vs.kr[w], vs.Pc[w]);

	SGFN.CalKrgPcgo(vs.S[g], vs.kr[g], vs.Pc[g]);

	SOF3.CalKrowKrog(vs.S[o], vs.krow, vs.krog);

	opC.CalOilPer(vs);
}


void OCPOGWFMethod02::CalKrPcDer(OCPFlowVarSet& vs)
{
	const INT& o = vs.o;
	const INT& g = vs.g;
	const INT& w = vs.w;

	SWFN.CalKrwPcwoDer(vs.S[w], vs.kr[w], vs.Pc[w], vs.dKrwdSw, vs.dPcwdSw);

	SGFN.CalKrgPcgoDer(vs.S[g], vs.kr[g], vs.Pc[g], vs.dKrgdSg, vs.dPcgdSg);

	SOF3.CalKrowKrogDer(vs.S[o], vs.krow, vs.krog, vs.dKrowdSo, vs.dKrogdSo);

	opC.CalOilPerDer(vs);

}


void OCPOGWFMethod02::Generate_SWPCWG()
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
// OCPFlowOGW
/////////////////////////////////////////////////////


void OCPFlowOGW::Setup(const ParamReservoir& rs_param, const USI& i)
{
	if (rs_param.SGOF_T.data.size() > 0 && rs_param.SWOF_T.data.size() > 0) {
		pfMethod = new OCPOGWFMethod01(rs_param.SGOF_T.data[i], rs_param.SWOF_T.data[i], 1, vs);
	}
	else if (rs_param.SOF3_T.data.size() > 0 &&
		rs_param.SGFN_T.data.size() > 0 &&
		rs_param.SWFN_T.data.size() > 0) {
		pfMethod = new OCPOGWFMethod02(rs_param.SOF3_T.data[i], rs_param.SGFN_T.data[i], 
			rs_param.SWFN_T.data[i], 1, vs);
	}
	else {
		OCP_ABORT("NO MATCHED METHOD!");
	}
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/08/2023      Create file                          */
/*----------------------------------------------------------------------------*/