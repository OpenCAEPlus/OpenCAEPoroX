/*! \file    OCPUnits.cpp
 *  \brief   Units used in OpenCAEPoro
 *  \author  Shizhe Li
 *  \date    Oct/22/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

 // OpenCAEPoroX header files
#include "OCPUnits.hpp"


void SetUnit(const string& ut)
{
	if (ut == "FIELD") {

		unitType        = UnitType::FIELD;

		CONV_DARCY      = FIELD_DARCY;
		GAS_CONSTANT    = FIELD_GAS_CONSTANT;
		GRAVITY_FACTOR  = FIELD_GRAVITY_FACTOR;
		RHOW_STD        = FIELD_RHOW_STD;
		RHOAIR_STD      = FIELD_RHOAIR_STD;
		PRESSURE_STD    = FIELD_PRESSURE_STD;
		TEMPERATURE_STD = FIELD_TEMPERATURE_STD;
		CONV1           = FIELD_CONV1;
		CONV2           = FIELD_CONV2;
		CONV3           = FIELD_CONV2 / FIELD_CONV1;
		CONV4           = FIELD_CONV3;
		CONV5           = FIELD_CONV4;
		CONV6           = FIELD_CONV5;
		TIMEUNIT        = FIELD_TIME;
	}
	else if (ut == "METRIC") {

		unitType        = UnitType::METRIC;

		CONV_DARCY      = METRIC_DARCY;
		GAS_CONSTANT    = METRIC_GAS_CONSTANT;
		GRAVITY_FACTOR  = METRIC_GRAVITY_FACTOR;
		RHOW_STD        = METRIC_RHOW_STD;
		RHOAIR_STD      = METRIC_RHOAIR_STD;
		PRESSURE_STD    = METRIC_PRESSURE_STD;
		TEMPERATURE_STD = METRIC_TEMPERATURE_STD;
		CONV1           = METRIC_CONV1;
		CONV2           = METRIC_CONV1;
		CONV3           = METRIC_CONV1;
		CONV4           = METRIC_CONV2;
		CONV5           = METRIC_CONV3;
		CONV6           = METRIC_CONV4;
		TIMEUNIT        = METRIC_TIME;
	}
	else if (ut == "SPE11A") {
		unitType        = UnitType::SPE11A;

		CONV_DARCY      = SPE11A_DARCY;
		GAS_CONSTANT    = SPE11A_GAS_CONSTANT;
		GRAVITY_FACTOR  = SPE11A_GRAVITY_FACTOR;
		RHOW_STD        = SPE11A_RHOW_STD;
		RHOAIR_STD      = SPE11A_RHOAIR_STD;
		PRESSURE_STD    = SPE11A_PRESSURE_STD;
		TEMPERATURE_STD = SPE11A_TEMPERATURE_STD;
		CONV1           = SPE11A_CONV1;
		CONV2           = SPE11A_CONV1;
		CONV3           = SPE11A_CONV1;
		CONV4           = SPE11A_CONV2;
		CONV5           = SPE11A_CONV3;
		CONV6           = SPE11A_CONV4;
		TIMEUNIT        = SPE11A_TIME;
	}
	else if (ut == "SPE11Amg") {
		unitType        = UnitType::SPE11Amg;

		CONV_DARCY      = SPE11Amg_DARCY;
		GAS_CONSTANT    = SPE11Amg_GAS_CONSTANT;
		GRAVITY_FACTOR  = SPE11Amg_GRAVITY_FACTOR;
		RHOW_STD        = SPE11Amg_RHOW_STD;
		RHOAIR_STD      = SPE11Amg_RHOAIR_STD;
		PRESSURE_STD    = SPE11Amg_PRESSURE_STD;
		TEMPERATURE_STD = SPE11Amg_TEMPERATURE_STD;
		CONV1           = SPE11Amg_CONV1;
		CONV2           = SPE11Amg_CONV1;
		CONV3           = SPE11Amg_CONV1;
		CONV4           = SPE11Amg_CONV2;
		CONV5           = SPE11Amg_CONV3;
		CONV6           = SPE11Amg_CONV4;
		TIMEUNIT        = SPE11Amg_TIME;
	}
	else {
		OCP_ABORT("Unit Type is not available!");
	}
}


/// OCP unit type
UnitType unitType = UnitType::FIELD;


////////////////////////////////////////////////
// Work Units
////////////////////////////////////////////////

/// Darcy constants
OCP_DBL CONV_DARCY       = FIELD_DARCY;
/// Gas constant
OCP_DBL GAS_CONSTANT     = FIELD_GAS_CONSTANT;
/// Gravity constant
OCP_DBL GRAVITY_FACTOR   = FIELD_GRAVITY_FACTOR;
/// Density of water at standard condition
OCP_DBL RHOW_STD         = FIELD_RHOW_STD;
/// Density of air at standard condition
OCP_DBL RHOAIR_STD       = FIELD_RHOAIR_STD;
/// Atmospheric pressure
OCP_DBL PRESSURE_STD     = FIELD_PRESSURE_STD;
/// Standard temperature
OCP_DBL TEMPERATURE_STD  = FIELD_TEMPERATURE_STD;
/// [bbl] -> [ft3], [m3] -> [m3]
OCP_DBL CONV1            = FIELD_CONV1;
/// [Mscf] -> [ft3], [m3] -> [m3]
OCP_DBL CONV2            = FIELD_CONV2;
/// [Mscf] -> [bbl], [m3] -> [m3]
OCP_DBL CONV3            = FIELD_CONV2 / FIELD_CONV1;
/// [F] -> [R], [C] -> [K]
OCP_DBL CONV4            = FIELD_CONV3;
/// [btu] -> [psia]¡¤[ft3], [kj] -> [bar]¡¤[m3]
OCP_DBL CONV5            = FIELD_CONV4;
/// [lbm]/[ft3] -> [gm-M]/[cc], [kg-M]/[m3] -> [gm-M]/[cc]
OCP_DBL CONV6            = FIELD_CONV5;

string TIMEUNIT          = FIELD_TIME;


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/22/2023      Create file                          */
/*----------------------------------------------------------------------------*/