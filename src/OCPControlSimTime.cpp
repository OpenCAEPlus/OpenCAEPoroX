/*! \file    OCPControlSimTime.cpp
 *  \brief   OCPControlSimTime class definition
 *  \author  Shizhe Li
 *  \date    Mar/22/2024
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPControlSimTime.hpp"

OCP_BOOL ControlSimTime::SetCurSimTime(const OCP_DBL& t)
{ 
	return t > maxSimTime;
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Mar/22/2024      Create file                          */
/*----------------------------------------------------------------------------*/