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
	virtual void CalFug(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* fug) const = 0;
	/// Calculate fugacity and fugacity coefficient
	virtual void CalFugPhi(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* fug, OCP_DBL* phi) const = 0;
	/// Calculate d(lnfug) / dx
	virtual void CalLnFugX(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* lnfugx) const = 0;
	/// Calculate d(lnfug) / dn
	virtual void CalLnFugN(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, const OCP_DBL& nt, OCP_DBL* lnfugn) const = 0;
	/// Calculate d(lnfug) / dP
	virtual void CalLnFugP(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* lnfugP) const = 0;
	/// Calculate d(lnphi) / dn
	virtual void CalLnPhiN(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, const OCP_DBL& nt, OCP_DBL* lnphin) const = 0;

public:
	/// Calculate molar volume
	virtual OCP_DBL CalVm(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x) const = 0;
	/// Calculate molar volume and derivatives
	virtual OCP_DBL CalVmDer(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL& vmP, OCP_DBL* vmx) const = 0;

};


///////////////////////////////////////////////
// EoS_PR
///////////////////////////////////////////////


class EoS_PR : public EoS
{
public:
	EoS_PR(const ComponentParam& param, const USI& tarId);
	/// Calculate fugacity
	void CalFug(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* fug) const override;
	/// Calculate fugacity and fugacity coefficient
	void CalFugPhi(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* fug, OCP_DBL* phi) const override;
	/// Calculate d(lnfug) / dx
	void CalLnFugX(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* lnfugx) const override;
	/// Calculate d(lnfug) / dn
	void CalLnFugN(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, const OCP_DBL& nt, OCP_DBL* lnfugn) const override;
	/// Calculate d(lnfug) / dP
	void CalLnFugP(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* lnfugP) const override;
	/// Calculate d(lnphi) / dn
	void CalLnPhiN(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, const OCP_DBL& nt, OCP_DBL* lnphin) const override;

public:
	/// Calculate molar volume
	OCP_DBL CalVm(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x) const override;
	/// Calculate molar volume and derivatives
	OCP_DBL CalVmDer(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL& vmP, OCP_DBL* vmx) const override;


protected:
	/// Calculate Ai, Bi
	void CalAiBi(const OCP_DBL& P, const OCP_DBL& T) const;
	/// Calculate Aj, Bj
	void CalAjBj( const OCP_DBL* x) const;
	/// Calculate compressible factor
	void CalZj(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x) const;
	/// Calculate Aj Bj Zj
	void CalAjBjZj(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x) const;
	/// Calculate Ax Bx Zx
	void CalAxBxZx(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x) const;
	/// Calculate An Bn Zn
	void CalAnBnZn(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, const OCP_DBL& nt) const;
	/// Calculate Ap Bp Zp
	void CalApBpZp(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x) const;
	
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
	mutable vector<OCP_DBL> Ai;
	/// Auxliary variable B for components
	mutable vector<OCP_DBL> Bi;
	/// Auxliary variable
	mutable OCP_DBL         Aj;
	/// Auxliary variable
	mutable OCP_DBL         Bj;
	/// compressible factor
	mutable OCP_DBL         Zj;
	/// root container
	mutable vector<OCP_DBL> Z;

	/// Derivatives
	/// dAj / dx
	mutable vector<OCP_DBL> Ax;
	/// dBj / dx
	mutable vector<OCP_DBL> Bx;
	/// dZj / dx
	mutable vector<OCP_DBL> Zx;
	/// dAj / dn
	mutable vector<OCP_DBL> An;
	/// dBj / dn
	mutable vector<OCP_DBL> Bn;
	/// dZj / dn
	mutable vector<OCP_DBL> Zn;
	/// dAj / dP
	mutable OCP_DBL         Ap;
	/// dBj / dP
	mutable OCP_DBL         Bp;
	/// dZj / dP
	mutable OCP_DBL         Zp;
};


class EoSCalculation
{
public:
	EoSCalculation() = default;
	void Setup(const ComponentParam& param, const USI& tarId);
	/// Calculate fugacity
	void CalFug(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* fug) const {
		eos->CalFug(P, T, x, fug);
	}
	/// Calculate fugacity and fugacity coefficient
	void CalFugPhi(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* fug, OCP_DBL* phi) const {
		eos->CalFugPhi(P, T, x, fug, phi);
	}
		/// Calculate d(lnfug) / dx
	void CalLnFugX(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* lnfugx) const {
		eos->CalLnFugX(P, T, x, lnfugx);
	}
	/// Calculate d(lnfug) / dn
	void CalLnFugN(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x,  const OCP_DBL& nt, OCP_DBL* lnfugn) const {
		eos->CalLnFugN(P, T, x, nt, lnfugn);
	}
	/// Calculate d(lnfug) / dP
	void CalLnFugP(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* lnfugP) const {
		eos->CalLnFugP(P, T, x, lnfugP);
	}
	/// Calculate d(lnphi) / dn
	void CalLnPhiN(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, const OCP_DBL& nt, OCP_DBL* lnphin) const {
		eos->CalLnPhiN(P, T, x, nt, lnphin);
	}

public:
	/// Calculate molar volume
	OCP_DBL CalVm(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x) const { return eos->CalVm(P, T, x); }
	/// Calculate molar volume and derivatives
	OCP_DBL CalVmDer(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL& vmP, OCP_DBL* vmx) const {
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