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


OCP_DBL OCP_PVTW::CalBw(const OCP_DBL& Pw)
{
	table.GetCloseRow(0, Pw, data);
	const OCP_DBL Pref   = data[0];
	const OCP_DBL Bwref  = data[1];
	const OCP_DBL cBwref = data[2];
	return Bwref * (1 - cBwref * (Pw - Pref));
}

void OCP_PVTW::CalBwMuwDer(const OCP_DBL& Pw, OCP_DBL& bw, OCP_DBL& muw, OCP_DBL& dBwdPw, OCP_DBL& dMuwdPw)
{
	table.GetCloseRow(0, Pw, data);
	const OCP_DBL Pref    = data[0];
	const OCP_DBL Bwref   = data[1];
	const OCP_DBL cBwref  = data[2];
	const OCP_DBL Muwref  = data[3];
	const OCP_DBL cMuwref = data[4];
	bw      = Bwref * (1 - cBwref * (Pw - Pref));
	dBwdPw  = -cBwref * Bwref;
	muw     = Muwref * (1 + cMuwref * (Pw - Pref));
	dMuwdPw = cMuwref * cMuwref;
}

/////////////////////////////////////////////////////
// PVCO
/////////////////////////////////////////////////////






/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jun/18/2023      Create file                          */
/*----------------------------------------------------------------------------*/

