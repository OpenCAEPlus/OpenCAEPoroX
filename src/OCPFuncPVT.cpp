/*! \file    OCPFuncPVT.cpp
 *  \brief   Functions for PVT in OCP
 *  \author  Shizhe Li
 *  \date    Jun/18/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


 // OpenCAEPoroX header files
#include "OCPFuncPVT.hpp"


/////////////////////////////////////////////////////
// PVTW
/////////////////////////////////////////////////////


OCP_DBL OCP_PVTW::CalBw(const OCP_DBL& P) const
{
	table.GetCloseRow(0, P, data);
	const OCP_DBL Pref  = data[0];
	const OCP_DBL bref  = data[1];
	const OCP_DBL Cbref = data[2];
	return bref * (1 - Cbref * (P - Pref));
}


void OCP_PVTW::CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu,
	OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const
{
	OCP_DBL b, bp;
	CalBwMuwDer(P, b, mu, bp, muP);
	xi   = 1 / CONV1 / b;
	xiP  = -1 / CONV1 * bp / (b * b);

	rho  = stdRhoW / b;
	rhoP = -stdRhoW * bp / (b * b);
}


void OCP_PVTW::CalBwMuwDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu, OCP_DBL& bP, OCP_DBL& muP) const
{
	table.GetCloseRow(0, P, data);
	const OCP_DBL Pref   = data[0];
	const OCP_DBL bref   = data[1];
	const OCP_DBL Cbref  = data[2];
	const OCP_DBL muref  = data[3];
	const OCP_DBL Cmuref = data[4];
	b   = bref * (1 - Cbref * (P - Pref));
	bP  = -Cbref * bref;
	mu  = muref * (1 + Cmuref * (P - Pref));
	muP = Cmuref * muref;
}

/////////////////////////////////////////////////////
// PVCO
/////////////////////////////////////////////////////


OCP_DBL OCP_PVCO::CalRhoO(const OCP_DBL& P, const OCP_DBL& Pb) const
{
	table.Eval_All(0, Pb, data);
	const OCP_DBL rssat = data[1];
	const OCP_DBL bref  = data[2];
	const OCP_DBL Cbref = data[4];
	const OCP_DBL b     = bref * (1 - Cbref * (P - Pb));
	const OCP_DBL rhoO  = (stdRhoO + (1000 / CONV1) * rssat * stdRhoG) / b;
	return rhoO;
}


OCP_DBL OCP_PVCO::CalXiO(const OCP_DBL& P, const OCP_DBL& Pb) const
{
	table.Eval_All(0, Pb, data);
	const OCP_DBL rs    = data[1];
	const OCP_DBL bref  = data[2];
	const OCP_DBL Cbref = data[4];
	const OCP_DBL b     = bref * (1 - Cbref * (P - Pb));
	return (1 + rs) / (b * CONV1);
}


void OCP_PVCO::CalRhoXiMuRsDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu, OCP_DBL& rs,
	OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP, OCP_DBL& rsP) const
{
	OCP_DBL b, bP;
	CalRsBoMuoDer(P, b, rs, mu, bP, rsP, muP);

	xi  = (1 + rs) / (b * CONV1);
	xiP = (rsP * b - (1 + rs) * bP) / (b * b * CONV1);

	rho = (stdRhoO + (1000 / CONV1) * rs * stdRhoG) / b;
	rhoP = (1000 / CONV1) * stdRhoG * rsP / b - (stdRhoO + (1000 / CONV1) * rs * stdRhoG) * bP / (b * b);
}


void OCP_PVCO::CalRsBoMuoDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& rs, OCP_DBL& mu,
	OCP_DBL& bP, OCP_DBL& rsP, OCP_DBL& muP) const
{
	table.Eval_All(0, P, data, cdata);
	rs  = data[1];
	b   = data[2];
	mu  = data[3];
	rsP = cdata[1];
	bP  = cdata[2];
	muP = cdata[3];
}


void OCP_PVCO::CalRhoXiMuDer(const OCP_DBL& rs, const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu,
	OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP, OCP_DBL& rhoRs, OCP_DBL& xiRs, OCP_DBL& muRs) const
{
	OCP_DBL b, bP, bRs;
	CalBoMuoDer(rs, P, b, mu, bP, muP, bRs, muRs);

	xi    = (1 + rs) / (CONV1 * b);
	xiP   = -(1 + rs) * bP / (CONV1 * b * b);
	xiRs  = (b - (1 + rs) * bRs) / (CONV1 * b * b);

	rho   = (stdRhoO + (1000 / CONV1) * rs * stdRhoG) / b;	
	rhoP  = -(stdRhoO + (1000 / CONV1) * rs * stdRhoG) * bP / (b * b);	
	rhoRs = ((1000 / CONV1) * stdRhoG * b - (stdRhoO + (1000 / CONV1) * rs * stdRhoG) * bRs) / (b * b);
}


void OCP_PVCO::CalBoMuoDer(const OCP_DBL& rs, const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu,
						     OCP_DBL& bP, OCP_DBL& muP, OCP_DBL& bRs, OCP_DBL& muRs) const
{
	table.Eval_All(1, rs, data, cdata);
	const OCP_DBL Pref     = data[0];
	const OCP_DBL bref     = data[2];
	const OCP_DBL muref    = data[3];
	const OCP_DBL Cbref    = data[4];
	const OCP_DBL Cmuref   = data[5];
	const OCP_DBL PrefRs   = cdata[0];
	const OCP_DBL brefRs   = cdata[2];
	const OCP_DBL murefRs  = cdata[3];
	const OCP_DBL CbrefRs  = cdata[4];
	const OCP_DBL CmurefRs = cdata[5];

	b    = bref * (1 - Cbref * (P - Pref));
	bP   = -bref * Cbref;
	bRs  = brefRs * (1 - Cbref * (P - Pref)) + bref * (-CbrefRs * (P - Pref) + Cbref * PrefRs);

	mu   = muref * (1 + Cmuref * (P - Pref));
	muP  = muref * Cmuref;
	muRs = murefRs * (1 + Cmuref * (P - Pref)) + muref * (CmurefRs * (P - Pref) - Cmuref * PrefRs);

}


/////////////////////////////////////////////////////
// PVDG
/////////////////////////////////////////////////////


void OCP_PVDG::CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu,
	OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const
{
	OCP_DBL b, bp;
	CalBgMugDer(P, b, mu, bp, muP);
	xi   = 1 / CONV1 / b;
	xiP  = -1 / CONV1 * bp / (b * b);

	rho  = 1000 / CONV1 * stdRhoG / b;	
	rhoP = -1000 / CONV1 * stdRhoG * bp / (b * b);
}

OCP_DBL OCP_PVDG::CalBg(const OCP_DBL& P) const
{
	return table.Eval(0, P, 1);
}

void OCP_PVDG::CalBgMugDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu, OCP_DBL& bP, OCP_DBL& muP) const
{
	table.Eval_All(0, P, data, cdata);
	b   = data[1];
	mu  = data[2];
	bP  = cdata[1];
	muP = cdata[2];
}


/////////////////////////////////////////////////////
// PVDO
/////////////////////////////////////////////////////


OCP_DBL OCP_PVDO::CalBo(const OCP_DBL& P) const
{
	return table.Eval(0, P, 1);
}


void OCP_PVDO::CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu,
	OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const
{
	OCP_DBL b, bp;
	CalBoMuoDer(P, b, mu, bp, muP);
	xi   = 1 / CONV1 / b;
	xiP  = -1 / CONV1 * bp / (b * b);

	rho  = stdRhoO / b;	
	rhoP = -stdRhoO * bp / (b * b);
}


void OCP_PVDO::CalBoMuoDer(const OCP_DBL& P, OCP_DBL& bo, OCP_DBL& muo, OCP_DBL& dBodP, OCP_DBL& dMuodP) const
{
	table.Eval_All(0, P, data, cdata);
	bo     = data[1];
	muo    = data[2];
	dBodP  = cdata[1];
	dMuodP = cdata[2];
}


/////////////////////////////////////////////////////
// ViscosityCalculation
/////////////////////////////////////////////////////


ViscosityMethod01::ViscosityMethod01(const TableSet& ts)
{
	table.Setup(ts);
	nc = ts.colNum - 1;
	muc.resize(nc);
	mucP.resize(nc);
	mucT.resize(nc);
}


OCP_DBL ViscosityMethod01::CalViscosity(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi)
{
	OCP_DBL mu = 0;
	table.Eval(P, T - CONV5, muc);
	for (USI i = 0; i < nc; i++)
		mu += zi[i] * log(muc[i]);
	return exp(mu);
}


OCP_DBL ViscosityMethod01::CalViscosity(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* muz)
{
	OCP_DBL mu = 0;
	muP = 0;
	muT = 0;
	table.Eval(P, T - CONV5, muc, mucP, mucT);
	for (USI i = 0; i < nc; i++) {
		muz[i] = zi[i] * log(muc[i]);
		mu     += muz[i];
	}
	mu = exp(mu);

	for (USI i = 0; i < nc; i++) {
		muP    += zi[i] / muc[i] * mucP[i];
		muT    += zi[i] / muc[i] * mucT[i];
		muz[i] *= mu;
	}
	muP *= mu;
	muT *= mu;
		
	return mu;
}


ViscosityMethod02::ViscosityMethod02(const vector<OCP_DBL>& av, const vector<OCP_DBL>& bv)
{
	avisc = av;
	bvisc = bv;
	nc    = avisc.size();
}


OCP_DBL ViscosityMethod02::CalViscosity(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi)
{
	OCP_DBL mu = 0;
	for (USI i = 0; i < nc; i++) {
		muc[i]  = avisc[i] * exp(bvisc[i] / T); 
		mu      += zi[i] * log(muc[i]);
	}
	return exp(mu);
}


OCP_DBL ViscosityMethod02::CalViscosity(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* muz)
{
	OCP_DBL mu = 0;
	muP = 0;
	muT = 0;
	for (USI i = 0; i < nc; i++) {
		muc[i]  = avisc[i] * exp(bvisc[i] / T);
		mucT[i] = -muc[i] * (bvisc[i] / T) / T;
		muz[i]  = zi[i] * log(muc[i]);
		mu      += muz[i];
	}
	mu = exp(mu);

	for (USI i = 0; i < nc; i++) {
		muP    += zi[i] / muc[i] * mucP[i];
		muT    += zi[i] / muc[i] * mucT[i];
		muz[i] *= mu;
	}
	muP *= mu;
	muT *= mu;

	return mu;
}


ViscosityMethod03::ViscosityMethod03(const ComponentParam& param, const USI& tarId)
{
	nc = param.numCom;
	if (param.Tc.activity)    Tc = param.Tc.data[tarId];
	else                      OCP_ABORT("TCRIT hasn't been input!");
	if (param.Pc.activity)    Pc = param.Pc.data[tarId];
	else                      OCP_ABORT("PCRIT hasn't been input!");
	if (param.MW.activity)    MWC = param.MW.data[tarId];
	else                      OCP_ABORT("MW hasn't been input!");

	coef = param.LBCcoef;
	for (auto& lbc : coef)    lbc *= 10; 

	if (param.Vcvis.activity) Vcvis = param.Vcvis.data[tarId];
	else if (param.Zcvis.activity) {
		Vcvis.resize(nc);
		const vector<OCP_DBL> Zcvis = param.Zcvis.data[tarId];
		for (USI i = 0; i < nc; i++) {
			Vcvis[i] = GAS_CONSTANT * Zcvis[i] * Tc[i] / Pc[i];
		}
	}
	else {
		if (param.Vc.activity)
			Vcvis = param.Vc.data[tarId];
		else if (param.Zc.activity) {
			Vcvis.resize(nc);
			const vector<OCP_DBL> Zc = param.Zc.data[tarId];
			for (USI i = 0; i < nc; i++) {
				Vcvis[i] = 10.73159 * Zc[i] * Tc[i] / Pc[i];
			}
		}
		else
			OCP_ABORT("VCRIT or ZCRIT hasn't been input!");
	}		

	muAux.resize(5);
}


OCP_DBL ViscosityMethod03::CalViscosity(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi)
{

	//OCP_DBL MW = 0;
	//for (USI i = 0; i < nc; i++) MW += zi[i] * MWC[i];

	//OCP_DBL tmp;
	//OCP_DBL Tri;
	//OCP_DBL xijT;
	//OCP_DBL xijP;
	//OCP_DBL xijV;

	//fill(muAux.begin(), muAux.end(), 0.0);
	//xijT = 0;
	//xijP = 0;
	//xijV = 0;

	//for (USI i = 0; i < nc; i++) {
	//	tmp = 5.4402 * pow(Tc[i], 1.0 / 6) / sqrt(MW) / pow(Pc[i], 2.0 / 3);
	//	Tri = T / Tc[i];
	//	if (Tri <= 1.5) {
	//		tmp = 34 * 1E-5 * pow(Tri, 0.94) / tmp;
	//	}
	//	else {
	//		tmp = 17.78 * 1E-5 * pow((4.58 * Tri - 1.67), 0.625) / tmp;
	//	}
	//	muAux[0] += zi[i] * sqrt(MWC[i]) * tmp;
	//	muAux[1] += zi[i] * sqrt(MWC[i]);
	//	xijT     += zi[i] * Tc[i];
	//	xijP     += zi[i] * Pc[i];
	//	xijV     += zi[i] * Vcvis[i];
	//}
	//muAux[2] = 5.4402 * pow(xijT, 1.0 / 6) / sqrt(MW) / pow(xijP, 2.0 / 3);
	//muAux[3] = xi * xijV;

	//if (muAux[3] <= 0.18 && OCP_FALSE) {
	//	return muAux[0] / muAux[1] + 2.05 * 1E-4 * muAux[3] / muAux[2];
	//}
	//else {
	//	muAux[4] = muAux[3] * (muAux[3] * (muAux[3] * (coef[4] * muAux[3] + coef[3]) + coef[2]) + coef[1]) + coef[0];
	//	return muAux[0] / muAux[1] + 1E-4 * (pow(muAux[4], 4) - 1) / muAux[2];
	//}
}


OCP_DBL ViscosityMethod03::CalViscosity(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* muz)
{

}


void ViscosityCalculation::Setup(const ComponentParam& param, const USI& tarId)
{
	if (param.avisc.activity) {
		vM = new ViscosityMethod02(param.avisc.data[tarId], param.bvisc.data[tarId]);
	}
	else if (param.viscTab.data.size() > 0) {
		vM = new ViscosityMethod01(param.viscTab);
	}
	else {
		OCP_ABORT("WRONG Viscosity Calculation Params!");
	}
}


/////////////////////////////////////////////////////
// EnthalpyCalculation
/////////////////////////////////////////////////////


EnthalpyMethod01::EnthalpyMethod01(const OCP_DBL& Trefin, const vector<OCP_DBL>& cpl1in, const vector<OCP_DBL>& cpl2in,
	const vector<OCP_DBL>& cpl3in, const vector<OCP_DBL>& cpl4in)
{
	Tref = Trefin + CONV5;
	cpl1 = cpl1in;
	cpl2 = cpl2in;
	cpl3 = cpl3in;
	cpl4 = cpl4in;
	nc   = cpl1.size();
}


OCP_DBL EnthalpyMethod01::CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const
{
	OCP_DBL H = 0;
	OCP_DBL Hz;
	const OCP_DBL dT = T - Tref;		
	for (USI i = 0; i < nc; i++) {
		Hz = (cpl1[i] * dT + 1.0 / 2 * cpl2[i] * pow(dT, 2) 
			+ 1.0 / 3 * cpl3[i] * pow(dT, 3) + 1.0 / 4 * cpl4[i] * pow(dT, 4));
		H   += zi[i] * Hz;
	}
	return H;
}

OCP_DBL EnthalpyMethod01::CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& HT, OCP_DBL* Hz) const
{
	OCP_DBL H = 0;
	HT = 0;
	const OCP_DBL dT = T - Tref;
	for (USI i = 0; i < nc; i++) {

		Hz[i] = (cpl1[i] * dT + 1.0 / 2 * cpl2[i] * pow(dT, 2) 
			  + 1.0 / 3 * cpl3[i] * pow(dT, 3) + 1.0 / 4 * cpl4[i] * pow(dT, 4));

		H     += zi[i] * Hz[i];
		HT    += zi[i] * (cpl1[i] + cpl2[i] * dT + cpl3[i] * pow(dT, 2) + cpl4[i] * pow(dT, 3));
	}
	return H;
}


EnthalpyMethod02::EnthalpyMethod02(const OCP_DBL& Trefin, const vector<OCP_DBL>& Tcritin,
	const vector<OCP_DBL>& cpg1in, const vector<OCP_DBL>& cpg2in,
	const vector<OCP_DBL>& cpg3in, const vector<OCP_DBL>& cpg4in,
	const vector<OCP_DBL>& hvaprin, const vector<OCP_DBL>& hvrin,
	const vector<OCP_DBL>& evin)
{
	Tref  = Trefin + CONV5;
	Tcrit = Tcritin;
	cpg1  = cpg1in;
	cpg2  = cpg2in;
	cpg3  = cpg3in;
	cpg4  = cpg4in;
	hvapr = hvaprin;
	hvr   = hvrin;
	ev    = evin;
	nc    = cpg1.size();
}


OCP_DBL EnthalpyMethod02::CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const
{
	OCP_DBL H = 0;
	OCP_DBL Hz;
	const OCP_DBL dT = T - Tref;
	
	for (USI i = 0; i < nc; i++) {
		Hz = cpg1[i] * dT + 1.0 / 2 * cpg2[i] * pow(dT, 2) 
			+ 1.0 / 3 * cpg3[i] * pow(dT, 3) + 1.0 / 4 * cpg4[i] * pow(dT, 4);
		H   += zi[i] * Hz;
	}
	return H;
}


OCP_DBL EnthalpyMethod02::CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& HT, OCP_DBL* Hz) const
{
	OCP_DBL H = 0;
	HT = 0;
	const OCP_DBL dT = T - Tref;
	for (USI i = 0; i < nc; i++) {
		Hz[i] = cpg1[i] * dT + 1.0 / 2 * cpg2[i] * pow(dT, 2) 
			  + 1.0 / 3 * cpg3[i] * pow(dT, 3) + 1.0 / 4 * cpg4[i] * pow(dT, 4);
		H     += zi[i] * Hz[i];
		HT    += zi[i] * (cpg1[i] + cpg2[i] * dT +
				cpg3[i] * pow(dT, 2) + cpg4[i] * pow(dT, 3));

		if (T < Tcrit[i]) {
			Hz[i] -= hvr[i] * pow((Tcrit[i] - T), ev[i]);
			H     -= zi[i] * hvr[i] * pow((Tcrit[i] - T), ev[i]);
			HT    += zi[i] * hvr[i] * ev[i] * pow((Tcrit[i] - T), ev[i] - 1);
		}
	}
	return H;
}


void EnthalpyCalculation::Setup(const ComponentParam& param, const USI& tarId)
{
	if (param.cpl1.activity) {
		eM = new EnthalpyMethod01(param.Tref[tarId],
			param.cpl1.data[tarId], param.cpl2.data[tarId],
			param.cpl3.data[tarId], param.cpl4.data[tarId]);
	}
	else if (param.cpg1.activity) {
		eM = new EnthalpyMethod02(param.Tref[tarId], param.Tc.data[tarId],
			param.cpg1.data[tarId], param.cpg2.data[tarId],
			param.cpg3.data[tarId], param.cpg4.data[tarId],
			param.hvapr.data[tarId], param.hvr.data[tarId],
			param.ev.data[tarId]);
	}
	else {
		OCP_ABORT("WRONG Enthalpy Calculation Params!");
	}
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jun/18/2023      Create file                          */
/*----------------------------------------------------------------------------*/

