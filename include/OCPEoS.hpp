/*! \file    OCPEoS.hpp
 *  \brief   OCPEoS class declaration
 *  \author  Shizhe Li
 *  \date    Jul/23/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPEOS_HEADER__
#define __OCPEOS_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include <vector>

using namespace std;

class OCPEoS
{
public:
	OCPEoS() = default;
	virtual OCP_DBL SolEoS(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi) = 0;
};


class OCPEoS_PR : public OCPEoS
{
public:
	OCPEoS_PR(const ComponentParam& param, const USI& tarId);
	OCP_DBL SolEoS(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi) override;

protected:
	void CalAiBi(const OCP_DBL& P, const OCP_DBL& T);

protected:
	USI             nc;
	/// Critical pressure of components
	vector<OCP_DBL> Pc;
	/// Critical temperature of components
	vector<OCP_DBL> Tc;
	/// OMEGA_A of components
	vector<OCP_DBL> OmegaA;
	/// OMEGA_B of components
	vector<OCP_DBL> OmegaB;
	/// Acentric factor of components
	vector<OCP_DBL> acf;
	/// Binary interaction between components
	vector<OCP_DBL> BIC;
	/// 
	OCP_DBL Aj;
	OCP_DBL Bj;
	OCP_DBL Zj;
	/// Auxliary variable A for components
	vector<OCP_DBL> Ai;
	/// Auxliary variable B for components
	vector<OCP_DBL> Bi;


	const OCP_DBL delta1 = 2.41421356237;
	const OCP_DBL delta2 = -0.41421356237;
};

#endif /* end if __OCPEOS_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/23/2023      Create file                          */
/*----------------------------------------------------------------------------*/