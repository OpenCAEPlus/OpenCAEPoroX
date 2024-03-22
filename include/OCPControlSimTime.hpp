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
#include <mpi.h>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "ParamControl.hpp"
#include "OCPControlFast.hpp"
#include "UtilTiming.hpp"

using namespace std;


/// OCP simulation time controler, start from the first time step, end with the last time step
/// unit:second, 
class ControlSimTime
{
public:
	void SetNextSimTime(const OCP_DBL& t) { nextSimTime = t; }
	void Initialize();
	OCP_BOOL IfStop();
	OCP_DBL GetTotalSimTime() const { return totalSimTime; }
	OCP_DBL GetWaitingTime() const {return waitTime;}

protected:
	/// total simulation time
	OCP_DBL     totalSimTime;
	/// next simulation time
	OCP_DBL     nextSimTime{ 1E20 };
	/// current simulation time from last setup
	OCP_DBL     curSimTime{ 0 };
	/// total waiting time
	OCP_DBL     waitTime;
	/// timer
	GetWallTime timer;
};



#endif /* end if __OCPControlSimTime_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Mar/22/2024      Create file                          */
/*----------------------------------------------------------------------------*/