/*! \file    OCPUnits.hpp
 *  \brief   Units used in OpenCAEPoro
 *  \author  Shizhe Li
 *  \date    Oct/22/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OPENCAEPORO_UNITS_HEADER__
#define __OPENCAEPORO_UNITS_HEADER__

// OpenCAEPoroX header files
#include "OCPDataType.hpp"
#include "UtilError.hpp"

using namespace std;


enum class UnitType : USI
{
	/// field units
	FIELD,
	/// metric units
	METRIC,
	/// units for spe11a
	SPE11A,
	/// units for spe11a with mg
	SPE11Amg
};



void SetUnit(const string& ut);


/// OCP unit type
extern UnitType unitType;


////////////////////////////////////////////////
// Field Units
////////////////////////////////////////////////

/// 1 [bbl] = 5.61458 [ft3]
const OCP_DBL FIELD_CONV1           = 5.61458;
/// 1 [Mscf] = 1000 [ft3]
const OCP_DBL FIELD_CONV2           = 1000;
/// 0 [F] = 459.67 [R]
const OCP_DBL FIELD_CONV3           = 459.67;
/// 1 [btu] = 5.40395 [psia]，[ft3]
const OCP_DBL FIELD_CONV4           = 5.40395;
/// 1 [lbm]/[ft3] = 0.0160185 [gm-M]/[cc]
const OCP_DBL FIELD_CONV5           = 0.0160185;
/// Gas constant ([ft3]，[psi]/[R]/[lb-M])
const OCP_DBL FIELD_GAS_CONSTANT    = 10.73159;
/// Gravity constant ([ft2]，[psi]/[lb])
const OCP_DBL FIELD_GRAVITY_FACTOR  = 0.00694444;
/// Darcy constant (([ft3/day]，[cp])/([md]，[ft]，[psi]))
const OCP_DBL FIELD_DARCY           = 0.001127 * FIELD_CONV1;
/// Density of water at standard condition ([lb]/[ft3]) 
const OCP_DBL FIELD_RHOW_STD        = 62.3664;
/// Density of air at standard condition ([lb]/[ft3]) 
const OCP_DBL FIELD_RHOAIR_STD      = 0.076362;
/// Atmospheric pressure ([psia])
const OCP_DBL FIELD_PRESSURE_STD    = 14.7;
/// Standard temperature ([F])
const OCP_DBL FIELD_TEMPERATURE_STD = 60;       
/// Time unit
const string  FIELD_TIME            = "Day";


////////////////////////////////////////////////
// Metric Units
////////////////////////////////////////////////

/// 1 [m3] = 1 [m3]
const OCP_DBL METRIC_CONV1           = 1.0;
/// 0 [C] = 273.15 [K]
const OCP_DBL METRIC_CONV2           = 273.15;
/// 1 [kj] = 0.01 [bar]，[m3]
const OCP_DBL METRIC_CONV3           = 0.01;
/// 1 [kg-M]/[m3] = 0.001 [gm-M]/[cc]
const OCP_DBL METRIC_CONV4           = 0.001;
/// Gas constant ([m3]，[bars]/[K]/[kg-M])
const OCP_DBL METRIC_GAS_CONSTANT    = 0.083143;
/// Gravity constant ([m2]，[bars]/[kg])
const OCP_DBL METRIC_GRAVITY_FACTOR  = 0.0000980665;
/// Darcy constant  (([m3/day]，[cP])/([md]，[m]，[bars]))
// const OCP_DBL METRIC_DARCY           = 0.008527;
const OCP_DBL METRIC_DARCY           = 1.0;
/// Density of water at standard condition ([kg]/[m3]) 
const OCP_DBL METRIC_RHOW_STD        = 999.014;
/// Density of air at standard condition ([kg]/[m3]) 
const OCP_DBL METRIC_RHOAIR_STD      = 1.2232;
/// Atmospheric pressure ([barsa])
const OCP_DBL METRIC_PRESSURE_STD    = 1.01325;
/// Standard temperature ([C])
const OCP_DBL METRIC_TEMPERATURE_STD = 20;
/// Time unit
const string  METRIC_TIME           = "Day";

////////////////////////////////////////////////
// Units for SPE11A (pa,g, 1 bar = 1E5 pa)
////////////////////////////////////////////////

/// 1 [cm3] = 1 [cm3]
const OCP_DBL SPE11A_CONV1           = 1.0;
/// 0 [C] = 273.15 [K]
const OCP_DBL SPE11A_CONV2           = 273.15;
/// 1 [kj] = 0.01 * 1E5 [pa]，[m3]
const OCP_DBL SPE11A_CONV3           = 0.01 * 1E5;
/// 1 [kg-M]/[m3] = 0.001 [gm-M]/[cc]
const OCP_DBL SPE11A_CONV4           = 0.001;
/// Gas constant ([m3]，[pa]/[K]/[kg-M])
const OCP_DBL SPE11A_GAS_CONSTANT    = 0.083143 * 1E5;
/// Gravity constant ([cm2]，[pa]/[g])  (REQUIRED)
const OCP_DBL SPE11A_GRAVITY_FACTOR  = 0.0000980665 * 1E5 *1E1;
/// Darcy constant (s，([cm2]/[cm2])/(pa，s)，pa), thus  (REQUIRED)
const OCP_DBL SPE11A_DARCY           = 1;
/// Density of water at standard condition ([g]/[cm3]) 
const OCP_DBL SPE11A_RHOW_STD        = 999.014 * 1E-3;
/// Density of air at standard condition ([g]/[cm3]) 
const OCP_DBL SPE11A_RHOAIR_STD      = 1.2232 * 1E-3;
/// Atmospheric pressure ([pa])  (REQUIRED)
const OCP_DBL SPE11A_PRESSURE_STD    = 1.01325 * 1E5;
/// Standard temperature ([C])
const OCP_DBL SPE11A_TEMPERATURE_STD = 20;
/// Time unit  (REQUIRED)
const string  SPE11A_TIME            = "Sec";


////////////////////////////////////////////////
// Units for SPE11A with mg (pa,g, 1 bar = 1E5 pa)
////////////////////////////////////////////////

/// 1 [cm3] = 1 [cm3]
const OCP_DBL SPE11Amg_CONV1 = 1.0;
/// 0 [C] = 273.15 [K]
const OCP_DBL SPE11Amg_CONV2 = 273.15;
/// 1 [kj] = 0.01 * 1E5 [pa]，[m3]
const OCP_DBL SPE11Amg_CONV3 = 0.01 * 1E5;
/// 1 [kg-M]/[m3] = 0.001 [gm-M]/[cc]
const OCP_DBL SPE11Amg_CONV4 = 0.001;
/// Gas constant ([m3]，[pa]/[K]/[kg-M])
const OCP_DBL SPE11Amg_GAS_CONSTANT = 0.083143 * 1E5;
/// Gravity constant ([cm2]，[pa]/[mg])  (REQUIRED)
const OCP_DBL SPE11Amg_GRAVITY_FACTOR = 0.0000980665 * 1E5 * 1E1 * 1E-3;
/// Darcy constant (s，([cm2]/[cm2])/(pa，s)，pa), thus  (REQUIRED)
const OCP_DBL SPE11Amg_DARCY = 1;
/// Density of water at standard condition ([mg]/[cm3]) 
const OCP_DBL SPE11Amg_RHOW_STD = 999.014;
/// Density of air at standard condition ([mg]/[cm3]) 
const OCP_DBL SPE11Amg_RHOAIR_STD = 1.2232;
/// Atmospheric pressure ([pa])  (REQUIRED)
const OCP_DBL SPE11Amg_PRESSURE_STD = 1.01325 * 1E5;
/// Standard temperature ([C])
const OCP_DBL SPE11Amg_TEMPERATURE_STD = 20;
/// Time unit  (REQUIRED)
const string  SPE11Amg_TIME = "Sec";



////////////////////////////////////////////////
// Work Units
////////////////////////////////////////////////

// physical variable
/// Darcy constants
extern OCP_DBL CONV_DARCY;
/// Gas constant
extern OCP_DBL GAS_CONSTANT;
/// Gravity constant
extern OCP_DBL GRAVITY_FACTOR;
/// Density of water at standard condition
extern OCP_DBL RHOW_STD;
/// Density of air at standard condition
extern OCP_DBL RHOAIR_STD;
/// Atmospheric pressure
extern OCP_DBL PRESSURE_STD;
/// Standard temperature
extern OCP_DBL TEMPERATURE_STD;

// unit conversion
/// [bbl] -> [ft3], [m3] -> [m3]
extern OCP_DBL CONV1;
/// [Mscf] -> [ft3], [m3] -> [m3]
extern OCP_DBL CONV2;
/// [Mscf] -> [bbl], [m3] -> [m3]
extern OCP_DBL CONV3;
/// [F] -> [R], [C] -> [K]
extern OCP_DBL CONV4;
/// [btu] -> [psia]，[ft3], [kj] -> [bar]，[m3]
extern OCP_DBL CONV5;
/// [lbm]/[ft3] -> [gm-M]/[cc], [kg-M]/[m3] -> [gm-M]/[cc]
extern OCP_DBL CONV6;

extern string TIMEUNIT;

#endif // __OPENCAEPORO_UNITS_HEADER__

  /*----------------------------------------------------------------------------*/
  /*  Brief Change History of This File                                         */
  /*----------------------------------------------------------------------------*/
  /*  Author              Date             Actions                              */
  /*----------------------------------------------------------------------------*/
  /*  Shizhe Li           Oct/22/2023      Create file                          */
  /*----------------------------------------------------------------------------*/