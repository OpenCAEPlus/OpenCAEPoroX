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


void ControlSimTime::Initialize()
{ 
	timer.Start(); 
	totalSimTime = 0;
	curSimTime   = 0;
	waitTime     = 0;
}


OCP_BOOL ControlSimTime::IfStop()
{ 
	curSimTime += timer.Stop();

	if (curSimTime > fabs(nextSimTime)) {

		totalSimTime += curSimTime;
		curSimTime   = 0;

		if (nextSimTime < 0) {
			timer.Start();
			if (CURRENT_RANK == MASTER_PROCESS) {
				cout << "Input the time(sec) to continue running(check again if negative): " << endl;
				cin >> nextSimTime;
			}
			MPI_Bcast(&nextSimTime, 1, OCPMPI_DBL, MASTER_PROCESS, MPI_COMM_WORLD);

			waitTime += timer.Stop();
			timer.Start();
			return OCP_FALSE;
		}
		else {
			return OCP_TRUE;
		}
	}
	else {
		timer.Start();
		return OCP_FALSE;
	}
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Mar/22/2024      Create file                          */
/*----------------------------------------------------------------------------*/