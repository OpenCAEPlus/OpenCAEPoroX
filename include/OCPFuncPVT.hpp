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
	void Setup(const vector<vector<OCP_DBL>>& src, const OCP_DBL& stdRhoWin) {
		OCPFuncTable::Setup(src);
		stdRhoW = stdRhoWin;
	}

	OCP_DBL CalXiW(const OCP_DBL& P) const { return 1 / CONV1 / CalBw(P);}
	OCP_DBL CalRhoW(const OCP_DBL& P) const { return stdRhoW / CalBw(P); }
	void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu, 
		               OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const;

protected:
	OCP_DBL CalBw(const OCP_DBL& P) const;
	void CalBwMuwDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu, OCP_DBL& bP, OCP_DBL& muP) const;
	
protected:
	OCP_DBL stdRhoW;   ///< mass density of water phase in standard condition
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
	void Setup(const vector<vector<OCP_DBL>>& src, const OCP_DBL& stdRhoOin, 
												   const OCP_DBL& stdRhoGin) {
		OCPFuncTable::Setup(src);
		stdRhoO = stdRhoOin;
		stdRhoG = stdRhoGin;	
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
	OCP_DBL stdRhoO;      ///< mass density of oil phase in standard condition
	OCP_DBL stdRhoG;      ///< mass density of gas phase in standard condition
	OCP_DBL stdVo{ 1 };   ///< molar volume of oil in standard condition
	OCP_DBL stdVg{ 1 };   ///< molar volume of gas in standard condition
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
	void Setup(const vector<vector<OCP_DBL>>& src, const OCP_DBL& stdRhoGin) {
		OCPFuncTable::Setup(src);
		stdRhoG = stdRhoGin;
	}
	OCP_DBL CalXiG(const OCP_DBL& P) const { return 1 / CONV1 / CalBg(P); }
	OCP_DBL CalRhoG(const OCP_DBL& P) const { return (1000 / CONV1) * stdRhoG / CalBg(P);}
	void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu,
					   OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const;

protected:
	OCP_DBL CalBg(const OCP_DBL& P) const;
	void CalBgMugDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu, OCP_DBL& bP, OCP_DBL& muP) const;

protected:
	OCP_DBL stdRhoG;   ///< mass density of gas phase in standard condition
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
	void Setup(const vector<vector<OCP_DBL>>& src, const OCP_DBL& stdRhoOin) {
		OCPFuncTable::Setup(src);
		stdRhoO = stdRhoOin;
	}
	OCP_DBL CalXiO(const OCP_DBL& P) const { return 1 / CONV1 / CalBo(P); }
	OCP_DBL CalRhoO(const OCP_DBL& P) const { return stdRhoO / CalBo(P); }
	void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu,
		OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const;
	
protected:
	OCP_DBL CalBo(const OCP_DBL& P) const;
	void CalBoMuoDer(const OCP_DBL& P, OCP_DBL& bo, OCP_DBL& muo, OCP_DBL& dBodP, OCP_DBL& dMuodP) const;

protected:
	OCP_DBL stdRhoO;   ///< mass density of oil phase in standard condition
};


/////////////////////////////////////////////////////
// Viscosity Calculation
/////////////////////////////////////////////////////

class ViscosityMethod
{
public:
	ViscosityMethod() = default;
	virtual OCP_DBL CalViscosity(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi) = 0;
	virtual OCP_DBL CalViscosity(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* muz) = 0;

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
	ViscosityMethod01(const TableSet& ts);
	OCP_DBL CalViscosity(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi) override;
	OCP_DBL CalViscosity(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* muz) override;

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
	OCP_DBL CalViscosity(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi) override;
	OCP_DBL CalViscosity(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* muz) override;

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
	OCP_DBL CalViscosity(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi) override;
	OCP_DBL CalViscosity(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* muz) override;

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
	/// Auxiliary variables
	vector<OCP_DBL> muAux;
};


class ViscosityCalculation
{
public:
	ViscosityCalculation() = default;
	void Setup(const ComponentParam& param, const USI& tarId);
	OCP_DBL CalViscosity(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi) { 
		return vM->CalViscosity(P, T, zi); 
	}
	OCP_DBL CalViscosity(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* zi,
		OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* muz) {
		return vM->CalViscosity(P, T, zi, muP, muT, muz);
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

