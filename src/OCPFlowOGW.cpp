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


void OCP3POilPerMethod01::CalOilPer(OCP3PFVarSet* vs)
{
	vs->kro = vs->krocw 
		      * ((vs->krow / vs->krocw + vs->krw) * (vs->krog / vs->krocw + vs->krg) 
		      - (vs->krw + vs->krg));
	if (vs->kro < 0) vs->kro = 0;
}


void OCP3POilPerMethod01::CalOilPerDer(OCP3PFVarSet* vs)
{

	vs->kro = vs->krocw
		* ((vs->krow / vs->krocw + vs->krw) * (vs->krog / vs->krocw + vs->krg)
			- (vs->krw + vs->krg));

	if (vs->kro < 0) {
		vs->kro     = 0;
		vs->dKrodSo = 0;
		vs->dKrodSg = 0;
		vs->dKrodSw = 0;
	}
	else {
		vs->dKrodSo = vs->dKrowdSo * (vs->krog / vs->krocw + vs->krg) 
			          + vs->dKrogdSo * (vs->krow / vs->krocw + vs->krw);
		vs->dKrodSw = vs->krocw * ((vs->dKrowdSw / vs->krocw + vs->dKrwdSw) 
			          * (vs->krog / vs->krocw + vs->krg) - (vs->dKrwdSw));
		vs->dKrodSg = vs->krocw * ((vs->krow / vs->krocw + vs->krw)
			          * (vs->dKrogdSg / vs->krocw + vs->dKrgdSg) - (vs->dKrgdSg));
	}
}

/////////////////////////////////////////////////////
// OCP3POilPerMethod02
/////////////////////////////////////////////////////


void OCP3POILPerMethod02::CalOilPer(OCP3PFVarSet* vs)
{
	const OCP_DBL tmp = vs->Sg + vs->Sw - vs->Swco;
	if (tmp <= TINY) {
		vs->kro = vs->krocw;
	}
	else {
		vs->kro = (vs->Sg * vs->krog + (vs->Sw - vs->Swco) * vs->krow) / tmp;
	}	
}


void OCP3POILPerMethod02::CalOilPerDer(OCP3PFVarSet* vs)
{
	const OCP_DBL tmp = vs->Sg + vs->Sw - vs->Swco;
	if (tmp <= TINY) {
		vs->kro     = vs->krocw;
		vs->dKrodSo = 0;
		vs->dKrodSg = 0;
		vs->dKrodSw = 0;
	}
	else {
		vs->kro = (vs->Sg * vs->krog + (vs->Sw - vs->Swco) * vs->krow) / tmp;
		vs->dKrodSo = (vs->Sg * vs->dKrogdSo + (vs->Sw - vs->Swco) * vs->dKrowdSo) / tmp;
		vs->dKrodSg = (vs->krog + vs->Sg * vs->dKrogdSg - vs->kro) / tmp;
		vs->dKrodSw = (vs->krow + (vs->Sw - vs->Swco) * vs->dKrowdSw - vs->kro) / tmp;
	}
}


/////////////////////////////////////////////////////
// OCP3PFMethod01
/////////////////////////////////////////////////////


OCP3PFMethod01::OCP3PFMethod01(const vector<vector<OCP_DBL>>& SGOFin,
	const vector<vector<OCP_DBL>>& SWOFin,
	const USI& i, OCP3PFVarSet* vsin) 
{
	SGOF.Setup(SGOFin);
	SWOF.Setup(SWOFin);
	vs = vsin;
	vs->krocw = SWOF.GetKrocw();
	vs->Swco  = SWOF.GetSwco();

	Generate_SWPCWG();

	switch (i)
	{
	case 1:
		opMethod = new OCP3POilPerMethod01();
		break;
	case 2:
	default:
		opMethod = new OCP3POILPerMethod02();
		break;
	}
}


void OCP3PFMethod01::Generate_SWPCWG()
{
	const std::vector<OCP_DBL> Sw(SWOF.GetSw());
	std::vector<OCP_DBL>       Pcow(SWOF.GetPcow());

	for (USI i = 0; i < Sw.size(); i++) {
		OCP_DBL Pcgo = SGOF.CalPcgo(1 - Sw[i]);
		Pcow[i] += Pcgo; // Pcgw
	}

	SWPCGW.Setup(vector<vector<OCP_DBL>>{Sw, Pcow});
}


void OCP3PFMethod01::CalKrPc()
{
	SWOF.CalKrwKrowPcwo(vs->Sw, vs->krw, vs->krow, vs->Pcwo);

	SGOF.CalKrgKrogPcgo(vs->Sg, vs->krg, vs->krog, vs->Pcgo);

	opMethod->CalOilPer(vs);
}


void OCP3PFMethod01::CalKrPcDer()
{
	SWOF.CalKrwKrowPcwoDer(vs->Sw, vs->krw, vs->krow, vs->Pcwo, vs->dKrwdSw, vs->dKrowdSw, vs->dPcwodSw);

	SGOF.CalKrgKrogPcgoDer(vs->Sg, vs->krg, vs->krog, vs->Pcgo, vs->dKrgdSg, vs->dKrogdSg, vs->dPcgodSg);

	opMethod->CalOilPerDer(vs);
}

/////////////////////////////////////////////////////
// OCP3PFMethod01
/////////////////////////////////////////////////////

OCP3PFMethod02::OCP3PFMethod02(const vector<vector<OCP_DBL>>& SOF3in,
	const vector<vector<OCP_DBL>>& SGFNin,
	const vector<vector<OCP_DBL>>& SWFNin,
	const USI& i, OCP3PFVarSet* vsin)
{
	SOF3.Setup(SOF3in);	
	SGFN.Setup(SGFNin);
	SWFN.Setup(SWFNin);

	vs = vsin;
	vs->krocw = SOF3.GetKrocw();
	vs->Swco  = SWFN.GetSwco();

	Generate_SWPCWG();

	switch (i)
	{
	case 1:
		opMethod = new OCP3POilPerMethod01();
		break;
	case 2:
	default:
		opMethod = new OCP3POILPerMethod02();
		break;
	}
}


void OCP3PFMethod02::CalKrPc()
{

	SWFN.CalKrwPcwo(vs->Sw, vs->krw, vs->Pcwo);

	SGFN.CalKrgPcgo(vs->Sg, vs->krg, vs->Pcgo);

	SOF3.CalKrowKrog(vs->So, vs->krow, vs->krog);

	opMethod->CalOilPer(vs);
}


void OCP3PFMethod02::CalKrPcDer()
{

	SWFN.CalKrwPcwoDer(vs->Sw, vs->krw, vs->Pcwo, vs->dKrwdSw, vs->dPcwodSw);

	SGFN.CalKrgPcgoDer(vs->Sg, vs->krg, vs->Pcgo, vs->dKrgdSg, vs->dPcgodSg);

	SOF3.CalKrowKrogDer(vs->So, vs->krow, vs->krog, vs->dKrowdSo, vs->dKrogdSo);

	opMethod->CalOilPerDer(vs);

}


void OCP3PFMethod02::Generate_SWPCWG()
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
	if (rs_param.SGOF_T.data.size() > 0 && rs_param.SWOF_T.data.size()) {
		pfMethod = new OCP3PFMethod01(rs_param.SGOF_T.data[i], rs_param.SWOF_T.data[i], 1, &vs);
	}
	else if (rs_param.SOF3_T.data.size() > 0 &&
		rs_param.SGFN_T.data.size() > 0 &&
		rs_param.SWFN_T.data.size() > 0) {
		pfMethod = new OCP3PFMethod02(rs_param.SOF3_T.data[i], rs_param.SGFN_T.data[i], 
			rs_param.SWFN_T.data[i], 1, &vs);
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