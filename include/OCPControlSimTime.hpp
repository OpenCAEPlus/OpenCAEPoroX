/*! \file    OCPControlSimTime.hpp
 *  \brief   OCPControlSimTime class declaration
 *  \author  Shizhe Li
 *  \date    Mar/22/2024
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPCONTROLSIMTIME_HEADER__
#define __OCPCONTROLSIMTIME_HEADER__

 // Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "OCPTimeRecord.hpp"
#include "ParamControl.hpp"
#include "OCPControlFast.hpp"

using namespace std;


/// OCP simulation time controler, unit:second
class ControlSimTime
{
public:
	void SetMaxSimTime(const OCP_DBL& t) { maxSimTime = t; }
	OCP_BOOL SetCurSimTime(const OCP_DBL& t);

protected:
	/// current simulation time
	OCP_DBL  curSimTime{ 0 };
	/// max simulation time
	OCP_DBL  maxSimTime{ 1E20 };
};



#endif /* end if __OCPControlSimTime_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Mar/22/2024      Create file                          */
/*----------------------------------------------------------------------------*/