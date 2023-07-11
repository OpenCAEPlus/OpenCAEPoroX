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
	OCP_DBL CalBw(const OCP_DBL& Pw);
	void CalBwMuwDer(const OCP_DBL& Pw, OCP_DBL& bw, OCP_DBL& muw, OCP_DBL& dBwdPw, OCP_DBL& dMuwdPw);
	
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
};


/////////////////////////////////////////////////////
// PVDG
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
};


#endif // __OCPFUNCPVT_HEADER__

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jun/18/2023      Create file                          */
/*----------------------------------------------------------------------------*/

