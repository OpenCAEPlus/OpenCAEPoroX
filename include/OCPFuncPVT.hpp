/*! \file    OCPFuncPVT.hpp
 *  \brief   Functions for PVT in OCP
 *  \author  Shizhe Li
 *  \date    Jun/18/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPFUNCPVT_HEADER__
#define __OCPFUNCPVT_HEADER__

// OpenCAEPoroX header files
#include "OCPFuncTable.hpp"
#include "UtilMath.hpp"

using namespace std;


/////////////////////////////////////////////////////
// PVTW
/////////////////////////////////////////////////////


/// Water PVT functions
class OCP_PVTW : public OCPFuncTable
{
	/// 0th column: The reference pressure for items 1 and 3. (Pref), (barsa (METRIC), psia (FIELD))
	///             This reference pressure should be close to field pressures
	/// 1th column: The water formation volume factor at the reference pressure. Bw(Pref), (rm3/sm3(METRIC), rb/stb(FIELD))
	/// 2th column: The water compressibility C = -(dBw / dP) / Bw, (1/bars (METRIC), 1/psi (FIELD))
	/// 3th column: The water viscosity at the reference pressure ¦Ìw(Pref), (cP (METRIC), cP (FIELD))
	/// 4th column: The water ¡°viscosibility¡± Cv = (d¦Ìw / dP) / ¦Ìw, (1/bars (METRIC), 1/psi (FIELD))

public:
	/// default constructor
	OCP_PVTW() = default;
	void Setup(const vector<vector<OCP_DBL>>& src, const OCP_DBL& stdRhoWin, const OCP_DBL& stdVwin) {
		OCPFuncTable::Setup(src);
		stdRhoW = stdRhoWin;
		stdVw   = stdVwin;
	}

	OCP_DBL CalXiW(const OCP_DBL& P) const { return 1 / CONV1 / (CalBw(P) * stdVw);}
	OCP_DBL CalRhoW(const OCP_DBL& P) const { return stdRhoW / CalBw(P); }
	void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu, 
		               OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const;

protected:
	OCP_DBL CalBw(const OCP_DBL& P) const;
	void CalBwMuwDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu, OCP_DBL& bP, OCP_DBL& muP) const;
	
protected:
	/// mass density of water phase in standard condition
	OCP_DBL stdRhoW;
	/// molar volume of water phase in standard conditions (stb/lbmol)
	OCP_DBL stdVw;
};


/////////////////////////////////////////////////////
// PVCO
/////////////////////////////////////////////////////


/// PVT properties of live oil in compressibility form (with dissolved gas)
class OCP_PVCO : public OCPFuncTable
{
	/// 0th column: The bubble point pressure for oil with the dissolved gas-oil ratio given by column 1, (Pbub)
	///             barsa (METRIC), psia (FIELD)
	/// 1th column: The dissolved gas-oil ratio of saturated oil with bubble point given by column 0, (Rs)
	///             (sm3/sm3 (METRIC), Mscf/stb(FIELD))
	/// 2th column: The formation volume factor of saturated oil at Pbub, (rm3/sm3 (METRIC), rb/stb(FIELD))
	/// 3th column: The viscosity of saturated oil at Pbub. (¦Ìo(Pbub)), (cP (METRIC), cP (FIELD))
	/// 4th column: The compressibility of undersaturated oil with a dissolved gas-oil ratio given by column 1.
	///             (C = -(dBo / dP) / Bo), (1/bars (METRIC), 1/psi (FIELD))
	/// 5th column: The ¡°viscosibility¡±, or viscosity compressibility, of undersaturated oil with a dissolved gas-oil ratio given by column 1.
	///             (Cv = (d¦Ìo / dP) / ¦Ìo), (1/bars (METRIC), 1/psi (FIELD))

public:
	/// default constructor
	OCP_PVCO() = default;
	void Setup(const vector<vector<OCP_DBL>>& src, 
		const OCP_DBL& stdRhoOin, const OCP_DBL& stdRhoGin,
		const OCP_DBL& stdVoin,   const OCP_DBL& stdVgin) {
		OCPFuncTable::Setup(src);
		stdRhoO = stdRhoOin;
		stdRhoG = stdRhoGin;	
		stdVo   = stdVoin;
		stdVg   = stdVgin;
	}
	
	OCP_DBL CalRhoO(const OCP_DBL& P, const OCP_DBL& Pb) const;
	OCP_DBL CalXiO(const OCP_DBL& P, const OCP_DBL& Pb) const;
	/// For saturated oil
	void CalRhoXiMuRsDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu, OCP_DBL& rs,
		OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP, OCP_DBL& rsP) const;
	/// For unsaturated oil
	void CalRhoXiMuDer(const OCP_DBL& rs, const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu,
		OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP, OCP_DBL& rhoRs, OCP_DBL& xiRs, OCP_DBL& muRs) const;
	/// For saturated oil
	OCP_DBL CalRs(const OCP_DBL& P) const { return table.Eval(0, P, 1); }	

protected:
	/// For saturated oil
	void CalRsBoMuoDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& rs, OCP_DBL& mu,
		OCP_DBL& bP, OCP_DBL& rsP, OCP_DBL& muP) const;
	/// For unsaturated oil
	void CalBoMuoDer(const OCP_DBL& rs, const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu, 
					   OCP_DBL& bP, OCP_DBL& muP, OCP_DBL& bRs, OCP_DBL& muRs) const;

protected:
	/// mass density of oil phase in standard condition
	OCP_DBL stdRhoO;      
	/// mass density of gas phase in standard condition
	OCP_DBL stdRhoG;   
	/// molar volume of oil in standard condition (stb/lbmol)
	OCP_DBL stdVo;  
	/// molar volume of gas in standard condition (Mscf/lbmol)
	OCP_DBL stdVg;        
};


/////////////////////////////////////////////////////
// PVDG
/////////////////////////////////////////////////////


/// PVT properties of dry gas (no vaporized oil)
class OCP_PVDG : public OCPFuncTable
{
	/// 0th column: The gas phase pressure. (Pg), (barsa (METRIC), psia (FIELD))
	/// 1th column: The corresponding gas formation volume factor. (Bg), (rm3/sm3 (METRIC), rb / Mscf(FIELD))
	/// 2th column: The corresponding gas viscosity. (¦Ìg), (cP (METRIC), cP (FIELD)))

public:
	/// default constructor
	OCP_PVDG() = default;
	void Setup(const vector<vector<OCP_DBL>>& src, const OCP_DBL& stdRhoGin, const OCP_DBL& stdVgin) {
		OCPFuncTable::Setup(src);
		stdRhoG = stdRhoGin;
		stdVg   = stdVgin;
	}
	OCP_DBL CalXiG(const OCP_DBL& P) const { return 1 / CONV1 / (CalBg(P) * stdVg); }
	OCP_DBL CalRhoG(const OCP_DBL& P) const { return (1000 / CONV1) * stdRhoG / CalBg(P);}
	void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu,
					   OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const;

protected:
	OCP_DBL CalBg(const OCP_DBL& P) const;
	void CalBgMugDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu, OCP_DBL& bP, OCP_DBL& muP) const;

protected:
	/// mass density of gas phase in standard condition
	OCP_DBL stdRhoG;
	/// molar volume of gas phase in standard conditions (stb/lbmol)
	OCP_DBL stdVg;
};


/////////////////////////////////////////////////////
// PVDO
/////////////////////////////////////////////////////


/// PVT properties of dead oil (no dissolved gas)
class OCP_PVDO : public OCPFuncTable
{
	/// 0th column: The oil phase pressure. (Po), (barsa (METRIC), psia (FIELD))
	/// 1th column: The corresponding oil formation volume factor. (Bo), (rm3/sm3 (METRIC), rb/stb(FIELD))
	/// 2th column: The corresponding oil viscosity. (¦Ìo), (cP (METRIC), cP (FIELD)))

public:
	/// default constructor
	OCP_PVDO() = default;
	virtual void Setup(const vector<vector<OCP_DBL>>& src, const OCP_DBL& stdRhoOin, const OCP_DBL& stdVoin) {
		OCPFuncTable::Setup(src);
		stdRhoO = stdRhoOin;
		stdVo   = stdVoin;
	}
	virtual OCP_DBL CalXiO(const OCP_DBL& P) const { return 1 / CONV1 / (CalBo(P) * stdVo); }
	virtual OCP_DBL CalRhoO(const OCP_DBL& P) const { return stdRhoO / CalBo(P); }
	virtual void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu,
		OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const;
	
protected:
	OCP_DBL CalBo(const OCP_DBL& P) const;
	void CalBoMuoDer(const OCP_DBL& P, OCP_DBL& bo, OCP_DBL& muo, OCP_DBL& dBodP, OCP_DBL& dMuodP) const;

protected:
	/// mass density of oil phase in standard conditions
	OCP_DBL stdRhoO;   
	/// molar volume of oil phase in standard conditions (stb/lbmol)
	OCP_DBL stdVo;     
};


/////////////////////////////////////////////////////
// PVCDO
/////////////////////////////////////////////////////

/// PVT properties of dead oil (with constant compressibility)
class OCP_PVCDO : public OCP_PVDO
{
	// only 1 row is included
	/// 0th column: The reference pressure (Pref) for items 1 and 3. (Pref), (barsa (METRIC), psia (FIELD))
	/// 1th column: The oil formation volume factor at the reference pressure. (Bo), (rm3/sm3 (METRIC), rb/stb(FIELD))
	/// 2th column: The oil compressibility. (C),  1/bars (METRIC), 1/psi (FIELD)
	/// 3th column: The oil viscosity at the reference pressure. (¦Ìo), (cP (METRIC), cP (FIELD)))
	/// 4th column: The oil ¡°viscosibility¡±. (C¦Ìo), 1/bars (METRIC), 1/psi (FIELD),

public:
	/// default constructor
	OCP_PVCDO() = default;
	void Setup(const vector<vector<OCP_DBL>>& src, const OCP_DBL& stdRhoOin, const OCP_DBL& stdVoin) override {
		Pref  = src[0][0];
		Bref  = src[1][0];
		Cb    = src[2][0];
		muref = src[3][0];
		Cmu   = src[4][0];
		stdRhoO = stdRhoOin;
		stdVo   = stdVoin;
	}
	OCP_DBL CalXiO(const OCP_DBL& P) const override { return 1 / CONV1 / (CalBo(P) * stdVo); }
	OCP_DBL CalRhoO(const OCP_DBL& P) const override { return stdRhoO / CalBo(P); }
	void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu,
		OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const override;

protected:
	OCP_DBL CalBo(const OCP_DBL& P) const;
	void CalBoMuoDer(const OCP_DBL& P, OCP_DBL& bo, OCP_DBL& muo, OCP_DBL& dBodP, OCP_DBL& dMuodP) const;

protected:
	/// reference pressure
	OCP_DBL Pref;
	/// oil formation volume factor at the reference pressure
	OCP_DBL Bref;
	/// The oil compressibility
	OCP_DBL Cb;
	/// The oil viscosity at the reference pressure
	OCP_DBL muref;
	/// The oil ¡°viscosibility¡±
	OCP_DBL Cmu;
};


/////////////////////////////////////////////////////
// Viscosity Calculation
/////////////////////////////////////////////////////

class ViscosityParams
{
public:
	ViscosityParams(const OCP_DBL* Pin, const OCP_DBL* Tin, const OCP_DBL* xin) :
		P(Pin), T(Tin), x(xin) {};
	ViscosityParams(const OCP_DBL* Pin, const OCP_DBL* Tin, const OCP_DBL* xin,
		const OCP_DBL* xiin) :
		P(Pin), T(Tin), x(xin), xi(xiin) {};
	ViscosityParams(const OCP_DBL* Pin, const OCP_DBL* Tin, const OCP_DBL* xin,
		            const OCP_DBL* xiin, const OCP_DBL* xiPin, const OCP_DBL* xiTin,
		            const OCP_DBL* xixin) :
		P(Pin), T(Tin), x(xin), xi(xiin), xiP(xiPin), xiT(xiTin), xix(xixin) {};
public:
	/// Pressure
	const OCP_DBL* P{ nullptr };
	/// Temperature
	const OCP_DBL* T{ nullptr };
	/// molar fraction
	const OCP_DBL* x{ nullptr };
	/// molar density
	const OCP_DBL* xi{ nullptr };
	/// d xi / d P
	const OCP_DBL* xiP{ nullptr };
	/// d xi / d T
	const OCP_DBL* xiT{ nullptr };
	/// d xi / d x
	const OCP_DBL* xix{ nullptr };
};

class ViscosityMethod
{
public:
	ViscosityMethod() = default;
	virtual OCP_DBL CalViscosity(const ViscosityParams& vp) = 0;
	virtual OCP_DBL CalViscosity(const ViscosityParams& vp, OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) = 0;

protected:
	/// num of components
	USI              nc;
	/// viscosity of components
	vector<OCP_DBL>  muc;
	/// d muC / dP 
	vector<OCP_DBL>  mucP;
	/// d muC / dT 
	vector<OCP_DBL>  mucT;
};


/// Use temperature and pressure-dependent viscosity table (for liquid viscosities)
/// Use linear mixing, T is in F(Field)
class ViscosityMethod01 : public ViscosityMethod
{
public:
	ViscosityMethod01(const Table2& tab);
	OCP_DBL CalViscosity(const ViscosityParams& vp) override;
	OCP_DBL CalViscosity(const ViscosityParams& vp, OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) override;

protected:
	/// Viscosity-versus-temperature&pressure dependence
	OCPFuncTable2    table;
};


/// Use tow params: avisc and bvisc (for liquid viscosities)
/// Use linear mixing, T is in R(Field)
class ViscosityMethod02 : public ViscosityMethod
{
public:
	ViscosityMethod02(const vector<OCP_DBL>& av, const vector<OCP_DBL>& bv);
	OCP_DBL CalViscosity(const ViscosityParams& vp) override;
	OCP_DBL CalViscosity(const ViscosityParams& vp, OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) override;

protected:
	/// Coefficients in water and oil viscosity correlation formula
	vector<OCP_DBL>   avisc;
	/// Coefficients in water and oil viscosity correlation formula
	vector<OCP_DBL>   bvisc;
};

/// Lohrenz-Bray-Clark formula 
class ViscosityMethod03 : public ViscosityMethod
{
public:
	ViscosityMethod03(const ComponentParam& param, const USI& tarId);
	OCP_DBL CalViscosity(const ViscosityParams& vp) override;
	OCP_DBL CalViscosity(const ViscosityParams& vp, OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) override;

protected:
	/// num of components
	USI             nc;
	/// LBC coefficients for viscosity calculation
	vector<OCP_DBL> coef;
	/// Critical temperature of components
	vector<OCP_DBL> Tc;     
	/// Critical pressure of components
	vector<OCP_DBL> Pc;
	/// Critical volume used for viscosity calculations only
	vector<OCP_DBL> Vcvis;
	/// Molecular Weight of components
	vector<OCP_DBL> MWC;
	/// Molecular Weight of phase
	OCP_DBL MW;
	/// Auxiliary variables
	vector<OCP_DBL> sqrtMWC;
	/// Auxiliary variables
	OCP_DBL xPc, xTc, xVc;
	/// Auxiliary variables
	vector<OCP_DBL> auxA;
	/// Auxiliary variables
	vector<OCP_DBL> auxB;
};


class ViscosityCalculation
{
public:
	ViscosityCalculation() = default;
	void Setup(const ComponentParam& param, const USI& tarId);
	OCP_DBL CalViscosity(const ViscosityParams& vp) {
		return vM->CalViscosity(vp); 
	}
	OCP_DBL CalViscosity(const ViscosityParams& vp, OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) {
		return vM->CalViscosity(vp, muP, muT, mux);
	}
protected:
	ViscosityMethod* vM;
};


/////////////////////////////////////////////////////
// EnthalpyCalculation
/////////////////////////////////////////////////////


class EnthalpyMethod
{
public:
	EnthalpyMethod() = default;
	virtual OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const = 0;
	virtual OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& HT, OCP_DBL* Hz) const = 0;
};


// liquid_based or simple_hvap, T(R Field)
class EnthalpyMethod01 : public EnthalpyMethod
{
public:
	EnthalpyMethod01(const OCP_DBL& Trefin, const vector<OCP_DBL>& cpl1in, const vector<OCP_DBL>& cpl2in,
				     const vector<OCP_DBL>& cpl3in, const vector<OCP_DBL>& cpl4in);
	OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override;
	OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& HT, OCP_DBL* Hz) const override;

protected:
	/// reference temperature
	OCP_DBL         Tref;
	/// num of components
	USI             nc;
	/// Coefficients in the component liquid enthalpy calculations, Btu/lbmol/F
	vector<OCP_DBL> cpl1;
	/// Coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^2
	vector<OCP_DBL> cpl2;
	/// Coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^3
	vector<OCP_DBL> cpl3;
	/// Coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^4
	vector<OCP_DBL> cpl4;
};


// gas_based, T(R Field)
class EnthalpyMethod02 : public EnthalpyMethod
{
public:
	EnthalpyMethod02(const OCP_DBL& Trefin, const vector<OCP_DBL>& Tcritin,
		const vector<OCP_DBL>& cpg1in, const vector<OCP_DBL>& cpg2in,
		const vector<OCP_DBL>& cpg3in, const vector<OCP_DBL>& cpg4in, 
		const vector<OCP_DBL>& hvaprin, const vector<OCP_DBL>& hvrin, 
		const vector<OCP_DBL>& evin);
	OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override;
	OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& HT, OCP_DBL* Hz) const override;

protected:
	/// reference temperature
	OCP_DBL         Tref;
	/// Critical temperature of hydrocarbon components
	vector<OCP_DBL> Tcrit;
	/// num of components
	USI             nc;
	/// Coefficients in the component liquid enthalpy calculations, Btu/lbmol/F
	vector<OCP_DBL> cpg1;
	/// Coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^2
	vector<OCP_DBL> cpg2;
	/// Coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^3
	vector<OCP_DBL> cpg3;
	/// Coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^4
	vector<OCP_DBL> cpg4;
	/// Coefficients in the component gas enthalpy calculations, Btu/lbmol
	vector<OCP_DBL> hvapr;
	/// Coefficients in the vaporization enthalpy calculations
	vector<OCP_DBL> hvr;
	/// Coefficients in the vaporization enthalpy calculations
	vector<OCP_DBL> ev;
};


class EnthalpyCalculation
{
public:
	void Setup(const ComponentParam& param, const USI& tarId);
	OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const { return eM->CalEnthalpy(T, zi); }
	OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& HT, OCP_DBL* Hz) const { return eM->CalEnthalpy(T, zi, HT, Hz); }
protected:
	EnthalpyMethod* eM;
};





#endif // __OCPFUNCPVT_HEADER__

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jun/18/2023      Create file                          */
/*----------------------------------------------------------------------------*/

