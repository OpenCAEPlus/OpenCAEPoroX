/*! \file    BoundaryCondition.hpp
 *  \brief   BoundaryCondition class declaration
 *  \author  Shizhe Li
 *  \date    Aug/25/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BOUNDARYCONDITION_HEADER__
#define __BOUNDARYCONDITION_HEADER__


// OpenCAEPoroX header files
#include "HeatLoss.hpp"


#include <vector>

using namespace std;


class BoundaryCondition
{
public:
	BoundaryCondition() = default;
	void Setup(const ParamReservoir& rs_param, const OCP_USI& nb) {
		heatLoss.Setup(rs_param, nb);
	}
	void CalHeatLoss(const BulkVarSet& bvs, const OCP_DBL& t, const OCP_DBL& dt) {
		heatLoss.CalHeatLoss(bvs, t, dt);
	}
	void ResetToLastTimeStep() { heatLoss.ResetToLastTimeStep(); }
	void UpdateLastTimeStep() { heatLoss.UpdateLastTimeStep(); }

public: // temp
	/// heat loss term
	HeatLoss heatLoss;
};


#endif /* end if __BOUNDARYCONDITION_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/25/2023      Create file                          */
/*----------------------------------------------------------------------------*/
