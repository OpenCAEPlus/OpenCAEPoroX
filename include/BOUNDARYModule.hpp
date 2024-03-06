/*! \file    BOUNDARYModule.hpp
 *  \brief   BOUNDARYModule class declaration
 *  \author  Shizhe Li
 *  \date    Sep/22/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BOUNDARYMODULE_HEADER__
#define __BOUNDARYMODULE_HEADER__


// OpenCAEPoroX header files
#include "ParamReservoir.hpp"
#include "BulkVarSet.hpp"
#include "HeatLoss.hpp"
#include "BoundaryFlow.hpp"


#include <vector>

using namespace std;


class BoundaryVarSet
{
public:
	vector<string>  boundName;
	vector<USI>     boundIndex;
	vector<OCP_DBL> boundArea;
};


class BOUNDARYModule
{
	friend class BulkAccumuTerm01;

public:
	void Setup(const ParamReservoir& rs_param, const BulkVarSet& bvs) {
		heatLoss.Setup(rs_param, bvs, vs.boundIndex);
		boundaryFlow.Setup(rs_param, bvs, vs.boundName, vs.boundIndex);
	}
	void ResetToLastTimeStep() { heatLoss.ResetToLastTimeStep(); }
	void UpdateLastTimeStep() { heatLoss.UpdateLastTimeStep(); }


public:
	/// Heat loss term
	HeatLoss     heatLoss;
	/// Boundary flow term
	BoundaryFlow boundaryFlow;

public:
	auto& GetBoundaryIndex() { return vs.boundIndex; }
	auto& GetBoundaryArea() { return vs.boundArea; }
	auto GetBoundaryArea(const OCP_USI& n) const { return vs.boundArea[n]; }
	auto& GetBoundName() { return vs.boundName; }


protected:
	BoundaryVarSet  vs;
};



#endif /* end if __BOUNDARYModuleHEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/22/2023      Create file                          */
/*----------------------------------------------------------------------------*/
