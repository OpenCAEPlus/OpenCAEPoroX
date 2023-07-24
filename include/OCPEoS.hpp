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
#include "UtilMath.hpp"
#include "DenseMat.hpp"
#include "ParamReservoir.hpp"
#include <vector>

using namespace std;

class EoS
{
public:
	EoS() = default;
	/// Calculate fugacity
	virtual void CalFug(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x, 
	vector<OCP_DBL>& fug) = 0;
	/// Calculate fugacity and fugacity coefficient
	virtual void CalFugPhi(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x, 
						   vector<OCP_DBL>& fug, vector<OCP_DBL>& phi) = 0;
	/// Calculate d(lnfug) / dx
	virtual void CalFugX(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x,
		vector<OCP_DBL>& fugx) = 0;
	/// Calculate d(lnfug) / dn
	virtual void CalFugN(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x,
		vector<OCP_DBL>& fugn) = 0;
	/// Calculate d(lnfug) / dP
	virtual void CalFugP(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x,
		vector<OCP_DBL>& fugP) = 0;

public:
	/// Calculate molar volume
	virtual OCP_DBL CalVm(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x) = 0;
	/// Calculate molar volume and derivatives
	virtual OCP_DBL CalVmDer(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x, OCP_DBL& vmP, vector<OCP_DBL>& vmx) = 0;

};


///////////////////////////////////////////////
// EoS_PR
///////////////////////////////////////////////


class EoS_PR : public EoS
{
public:
	EoS_PR(const ComponentParam& param, const USI& tarId);
	/// Calculate fugacity
	void CalFug(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x,
		vector<OCP_DBL>& fug) override;
	/// Calculate fugacity and fugacity coefficient
	void CalFugPhi(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x,
				   vector<OCP_DBL>& fug, vector<OCP_DBL>& phi) override;
	/// Calculate d(lnfug) / dx
	void CalFugX(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x, 
				 vector<OCP_DBL>& fugx) override;
	/// Calculate d(lnfug) / dn
	void CalFugN(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& n,
		vector<OCP_DBL>& fugn) override;
	/// Calculate d(lnfug) / dP
	void CalFugP(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x,
		vector<OCP_DBL>& fugP) override;

public:
	/// Calculate molar volume
	OCP_DBL CalVm(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x) override;
	/// Calculate molar volume and derivatives
	OCP_DBL CalVmDer(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x, OCP_DBL& vmP, vector<OCP_DBL>& vmx) override;


protected:
	/// Calculate Ai, Bi
	void CalAiBi(const OCP_DBL& P, const OCP_DBL& T);
	/// Calculate Aj, Bj
	void CalAjBj(const vector<OCP_DBL>& x);
	/// Calculate compressible factor
	void CalZj(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x);
	/// Calculate Aj Bj Zj
	void CalAjBjZj(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x);
	/// Calculate Ax Bx Zx
	void CalAxBxZx(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x);
	/// Calculate An Bn Zn
	void CalAnBnZn(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x, const OCP_DBL& nt);
	/// Calculate Ap Bp Zp
	void CalApBpZp(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x);
	
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
	/// Volume shift of components
	vector<OCP_DBL> Vshift; 
	/// PR-EoS params
	const OCP_DBL delta1 = 2.41421356237;
	const OCP_DBL delta2 = -0.41421356237;

	/// Auxliary variable A for components
	vector<OCP_DBL> Ai;
	/// Auxliary variable B for components
	vector<OCP_DBL> Bi;
	/// Auxliary variable
	OCP_DBL Aj;
	/// Auxliary variable
	OCP_DBL Bj;
	/// compressible factor
	OCP_DBL Zj;
	/// root container
	vector<OCP_DBL> Z;

	/// Derivatives
	/// dAj / dx
	vector<OCP_DBL> Ax;
	/// dBj / dx
	vector<OCP_DBL> Bx;
	/// dZj / dx
	vector<OCP_DBL> Zx;
	/// dAj / dn
	vector<OCP_DBL> An;
	/// dBj / dn
	vector<OCP_DBL> Bn;
	/// dZj / dn
	vector<OCP_DBL> Zn;
	/// dAj / dP
	OCP_DBL         Ap;
	/// dBj / dP
	OCP_DBL         Bp;
	/// dZj / dP
	OCP_DBL         Zp;
	/// molar fraction
	vector<OCP_DBL> x;
};


class EoSCalculation
{
public:
	EoSCalculation() = default;
	void Setup(const ComponentParam& param, const USI& tarId);
	/// Calculate fugacity
	void CalFug(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x,
		vector<OCP_DBL>& fug) {
		eos->CalFug(P, T, x, fug);
	}
	/// Calculate fugacity and fugacity coefficient
	void CalFugPhi(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x,
		vector<OCP_DBL>& fug, vector<OCP_DBL>& phi) {
		eos->CalFugPhi(P, T, x, fug, phi);
	}
		/// Calculate d(lnfug) / dx
	void CalFugX(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x,
		vector<OCP_DBL>& fugx) {
		eos->CalFugX(P, T, x, fugx);
	}
	/// Calculate d(lnfug) / dn
	void CalFugN(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& n,
		vector<OCP_DBL>& fugn) {
		eos->CalFugN(P, T, n, fugn);
	}
	/// Calculate d(lnfug) / dP
	void CalFugP(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x,
		vector<OCP_DBL>& fugP) {
		eos->CalFugP(P, T, x, fugP);
	}

public:
	/// Calculate molar volume
	OCP_DBL CalVm(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x) { return eos->CalVm(P, T, x); }
	/// Calculate molar volume and derivatives
	OCP_DBL CalVmDer(const OCP_DBL& P, const OCP_DBL& T, const vector<OCP_DBL>& x, OCP_DBL& vmP, vector<OCP_DBL>& vmx) {
		return eos->CalVmDer(P, T, x, vmP, vmx);
	}


protected:
	EoS*  eos;
};


#endif /* end if __OCPEOS_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/23/2023      Create file                          */
/*----------------------------------------------------------------------------*/