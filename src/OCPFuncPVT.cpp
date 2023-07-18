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


OCP_DBL OCP_PVTW::CalBw(const OCP_DBL& P)
{
	table.GetCloseRow(0, P, data);
	const OCP_DBL Pref  = data[0];
	const OCP_DBL bref  = data[1];
	const OCP_DBL Cbref = data[2];
	return bref * (1 - Cbref * (P - Pref));
}


void OCP_PVTW::CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu,
	OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP)
{
	OCP_DBL b, bp;
	CalBwMuwDer(P, b, mu, bp, muP);
	xi   = 1 / CONV1 / b;
	xiP  = -1 / CONV1 * bp / (b * b);

	rho  = stdRhoW / b;
	rhoP = -stdRhoW * bp / (b * b);
}


void OCP_PVTW::CalBwMuwDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu, OCP_DBL& bP, OCP_DBL& muP)
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


OCP_DBL OCP_PVCO::CalRhoO(const OCP_DBL& P, const OCP_DBL& Pb)
{
	table.Eval_All(0, Pb, data);
	const OCP_DBL rssat = data[1];
	const OCP_DBL bref  = data[2];
	const OCP_DBL Cbref = data[4];
	const OCP_DBL b     = bref * (1 - Cbref * (P - Pb));
	const OCP_DBL rhoO  = (stdRhoO + (1000 / CONV1) * rssat * stdRhoG) / b;
	return rhoO;
}


OCP_DBL OCP_PVCO::CalXiO(const OCP_DBL& P)
{
	table.Eval_All(0, P, data);
	const OCP_DBL rs = data[1];
	const OCP_DBL b  = data[2];
	return (1 + rs) / (b * CONV1);
}


OCP_DBL OCP_PVCO::CalXiO(const OCP_DBL& P, const OCP_DBL& Pb)
{
	table.Eval_All(0, Pb, data);
	const OCP_DBL rs    = data[1];
	const OCP_DBL bref  = data[2];
	const OCP_DBL Cbref = data[4];
	const OCP_DBL b     = bref * (1 - Cbref * (P - Pb));
	return (1 + rs) / (b * CONV1);
}


void OCP_PVCO::CalRhoXiMuRsDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu, OCP_DBL& rs,
	OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP, OCP_DBL& rsP)
{
	OCP_DBL b, bP;
	CalRsBoMuoDer(P, b, rs, mu, bP, rsP, muP);

	xi  = (1 + rs) / (b * CONV1);
	xiP = (rsP * b - (1 + rs) * bP) / (b * b * CONV1);

	rho = (stdRhoO + (1000 / CONV1) * rs * stdRhoG) / b;
	rhoP = (1000 / CONV1) * stdRhoG * rsP / b - (stdRhoO + (1000 / CONV1) * rs * stdRhoG) * bP / (b * b);
}


void OCP_PVCO::CalRsBoMuoDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& rs, OCP_DBL& mu,
	OCP_DBL& bP, OCP_DBL& rsP, OCP_DBL& muP)
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
	OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP, OCP_DBL& rhoRs, OCP_DBL& xiRs, OCP_DBL& muRs)
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
						     OCP_DBL& bP, OCP_DBL& muP, OCP_DBL& bRs, OCP_DBL& muRs)
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
	OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP)
{
	OCP_DBL b, bp;
	CalBgMugDer(P, b, mu, bp, muP);
	xi   = 1 / CONV1 / b;
	xiP  = -1 / CONV1 * bp / (b * b);

	rho  = 1000 / CONV1 * stdRhoG / b;	
	rhoP = -1000 / CONV1 * stdRhoG * bp / (b * b);
}

OCP_DBL OCP_PVDG::CalBg(const OCP_DBL& P)
{
	return table.Eval(0, P, 1);
}

void OCP_PVDG::CalBgMugDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu, OCP_DBL& bP, OCP_DBL& muP)
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


OCP_DBL OCP_PVDO::CalBo(const OCP_DBL& P)
{
	return table.Eval(0, P, 1);
}


void OCP_PVDO::CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu,
	OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP)
{
	OCP_DBL b, bp;
	CalBoMuoDer(P, b, mu, bp, muP);
	xi   = 1 / CONV1 / b;
	xiP  = -1 / CONV1 * bp / (b * b);

	rho  = stdRhoO / b;	
	rhoP = -stdRhoO * bp / (b * b);
}


void OCP_PVDO::CalBoMuoDer(const OCP_DBL& P, OCP_DBL& bo, OCP_DBL& muo, OCP_DBL& dBodP, OCP_DBL& dMuodP)
{
	table.Eval_All(0, P, data, cdata);
	bo     = data[1];
	muo    = data[2];
	dBodP  = cdata[1];
	dMuodP = cdata[2];
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jun/18/2023      Create file                          */
/*----------------------------------------------------------------------------*/

