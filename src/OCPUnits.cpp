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


void SetUnit(const string& unitType)
{
	if (unitType == "FIELD") {

	}
	else if (unitType == "METRIC") {

	}
	else {
		OCP_ABORT("Unit Type is not available!");
	}
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/22/2023      Create file                          */
/*----------------------------------------------------------------------------*/