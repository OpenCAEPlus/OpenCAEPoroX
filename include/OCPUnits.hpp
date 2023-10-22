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


////////////////////////////////////////////////
// Field Units
////////////////////////////////////////////////

/// 1 [bbl] = 5.61458 [ft3]
const OCP_DBL FIELD_CONV1           = 5.61458;
/// 1 [btu] = 5.40395 [psia]，[ft3]
const OCP_DBL FIELD_CONV2           = 5.40395;
/// 1 [lbm]/[ft3] = 0.0160185 [gm-M]/[cc]
const OCP_DBL FIELD_CONV3           = 0.0160185;
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


////////////////////////////////////////////////
// Metric Units
////////////////////////////////////////////////

/// 1 [rm3] = 1 [sm3]
const OCP_DBL METRIC_CONV1           = 1.0;
/// 1 [kj] = 0.01 [bar]，[m3]
const OCP_DBL METRIC_CONV2           = 0.01;
/// 1 [kg-M]/[m3] = 0.0160185 [gm-M]/[cc]
const OCP_DBL METRIC_CONV3           = 0.001;
/// Gas constant ([m3]，[bars]/[K]/[kg-M])
const OCP_DBL METRIC_GAS_CONSTANT    = 0.083143;
/// Gravity constant ([m2]，[bars]/[kg])
const OCP_DBL METRIC_GRAVITY_FACTOR  = 0.0000980665;
/// Darcy constant  (([m3/day]，[cP])/([md]，[m]，[bars]))
const OCP_DBL METRIC_DARCY           = 0.008527 * METRIC_CONV1;
/// Density of water at standard condition ([kg]/[m3]) 
const OCP_DBL METRIC_RHOW_STD        = 999.014;
/// Density of air at standard condition ([kg]/[m3]) 
const OCP_DBL METRIC_RHOAIR_STD      = 1.2232;
/// Atmospheric pressure ([barsa])
const OCP_DBL METRIC_PRESSURE_STD    = 1.013;
/// Standard temperature ([C])
const OCP_DBL METRIC_TEMPERATURE_STD = 20;




// Physical consts
const OCP_DBL GAS_CONSTANT    = 10.73159;    ///< Gas Constant
const OCP_DBL GRAVITY_FACTOR  = 0.00694444; ///< 0.00694444 ft2 psi / lb
const OCP_DBL RHOW_STD        = 62.3664;    ///< Water density at surface cond: lb/ft3
const OCP_DBL RHOAIR_STD      = 0.076362;   ///< Air density at surface cond : lb/ft3
const OCP_DBL PRESSURE_STD    = 14.7;       ///< 14.6959 psia = 1 atm
const OCP_DBL TEMPERATURE_STD = 60;         ///< Standard temperature: F


// Unit conversion consts
const OCP_DBL CONV1 = 5.61458;               ///< 1 bbl = CONV1 ft3
const OCP_DBL CONV2 = 1.12712E-3;            ///< Darcy constant in field unit
const OCP_DBL CONV3 = 0.45359237;            ///< 1 lb = CONV3 kg
const OCP_DBL CONV4 = 0.02831685;            ///< 1 ft3 = CONV4 m3
const OCP_DBL CONV5 = 459.67;                ///< 0 F = CONV5 R
const OCP_DBL CONV6 = 778.172448;            ///< 1[Btu] = 778.172448 [ft]，[lbf]
const OCP_DBL CONV7 = CONV3 / (CONV4 * 1E3); ///< lbm/ft3 -> gm-M/cc

/// Darcy constants
const OCP_DBL CONV_DARCY = CONV1 * CONV2;


#endif // __OPENCAEPORO_UNITS_HEADER__

  /*----------------------------------------------------------------------------*/
  /*  Brief Change History of This File                                         */
  /*----------------------------------------------------------------------------*/
  /*  Author              Date             Actions                              */
  /*----------------------------------------------------------------------------*/
  /*  Shizhe Li           Oct/22/2023      Create file                          */
  /*----------------------------------------------------------------------------*/