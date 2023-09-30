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


/////////////////////////////////////////////////////
// Brooks-Corey
/////////////////////////////////////////////////////


void BrooksCorey::Setup(const BrooksCoreyParam& bcp)
{
	sw_imm = bcp.sw_imm;
	sn_imm = bcp.sn_imm;
	Pentry = bcp.Pentry;
	Pcmax  = bcp.Pcmax;
	Cw_kr  = bcp.Cw_kr;
	Cn_kr  = bcp.Cn_kr;
	C_pc   = bcp.C_pc;
}


OCP_DBL BrooksCorey::CalKrN(const OCP_DBL& sn) const
{
	if (sn < sn_imm) return 0;
	else             return pow((sn - sn_imm) / (1 - sn_imm), Cn_kr);
}


OCP_DBL BrooksCorey::CalKrW(const OCP_DBL& sw) const
{
	if (sw < sw_imm) return 0;
	else             return pow((sw - sw_imm) / (1 - sw_imm), Cw_kr);
}


void BrooksCorey::CalKrPcN(const OCP_DBL& sn, OCP_DBL& kr, OCP_DBL& pc) const
{
	pc = 0;
	kr = CalKrN(sn);
}


void BrooksCorey::CalKrPcW(const OCP_DBL& sw, OCP_DBL& kr, OCP_DBL& pc) const
{

	if (sw < sw_imm) {
		kr = 0;
		pc = Pcmax;
	}
	else {
		const OCP_DBL swn = (sw - sw_imm) / (1 - sw_imm);
		kr = pow(swn, Cw_kr);
		pc = Pentry * pow(swn, (-1.0 / C_pc));
		pc = Pcmax * erf(pc / Pcmax * sqrt(PI) / 2);
	}
}


void BrooksCorey::CalKrPcDerN(const OCP_DBL& sn, OCP_DBL& kr, OCP_DBL& pc, OCP_DBL& dKrdSn, OCP_DBL& dPcdSn) const
{
	if (sn < sn_imm) {
		kr     = 0;
		dKrdSn = 0;
	}
	else {
		const OCP_DBL snn = (sn - sn_imm) / (1 - sn_imm);
		kr     = pow(snn, Cn_kr);
		dKrdSn = Cn_kr * pow(snn, Cn_kr - 1) / (1 - sn_imm);
	}
	pc     = 0;
	dPcdSn = 0;
}


void BrooksCorey::CalKrPcDerW(const OCP_DBL& sw, OCP_DBL& kr, OCP_DBL& pc, OCP_DBL& dKrdSw, OCP_DBL& dPcdSw) const
{
	if (sw < sw_imm) {
		kr     = 0;
		dKrdSw = 0;
		pc     = Pcmax;
		dPcdSw = 0;

	}
	else {
		const OCP_DBL swn = (sw - sw_imm) / (1 - sw_imm);
		kr     = pow(swn, Cw_kr);
		dKrdSw = Cw_kr * pow(swn, Cw_kr - 1) / (1 - sw_imm);
		const OCP_DBL tmp = Pentry * pow(swn, (-1.0 / C_pc)) / Pcmax * sqrt(PI) / 2;
		pc                = Pcmax * erf(tmp);
		dPcdSw            = Pentry * exp(-tmp * tmp) * (-1.0 / C_pc) * pow(swn, (-1.0 / C_pc - 1)) / (1 - sw_imm);
	}
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jun/29/2023      Create file                          */
/*----------------------------------------------------------------------------*/

