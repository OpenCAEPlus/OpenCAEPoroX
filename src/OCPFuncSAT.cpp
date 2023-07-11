/*! \file    OCPFuncSAT.cpp
 *  \brief   Functions for Saturation in OCP
 *  \author  Shizhe Li
 *  \date    Jun/29/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


 // OpenCAEPoroX header files
#include "OCPFuncSAT.hpp"

/////////////////////////////////////////////////////
// SWOF
/////////////////////////////////////////////////////


OCP_DBL OCP_SWOF::GetSwcr() const
{
	 const vector<OCP_DBL>& Sw  = table.GetCol(0);
	 const vector<OCP_DBL>& krw = table.GetCol(1);
	 for (USI i = 0; i < krw.size(); i++) {
		 if (krw[i] >= TINY) {
			 return Sw[i];
		 }
	 }
	 OCP_ABORT("WATER is immobilable in SWOF!");
	 return -1;
}


/////////////////////////////////////////////////////
// SGOF
/////////////////////////////////////////////////////


OCP_DBL OCP_SGOF::GetSgcr() const
{
	const vector<OCP_DBL>& Sg = table.GetCol(0);
	const vector<OCP_DBL>& krg = table.GetCol(1);
	for (USI i = 0; i < krg.size(); i++) {
		if (krg[i] >= TINY) {
			return Sg[i];
		}
	}
	OCP_ABORT("GAS is immobilable in SGOF!");
	return -1;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jun/29/2023      Create file                          */
/*----------------------------------------------------------------------------*/

