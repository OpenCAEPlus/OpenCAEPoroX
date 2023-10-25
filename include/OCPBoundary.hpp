/*! \file    OCPBoundary.hpp
 *  \brief   OCPBoundary class declaration
 *  \author  Shizhe Li
 *  \date    Sep/22/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPBOUNDARY_HEADER__
#define __OCPBOUNDARY_HEADER__


// OpenCAEPoroX header files
#include "ParamReservoir.hpp"
#include "BulkVarSet.hpp"
#include "HeatLoss.hpp"
#include "BoundaryFlow.hpp"


#include <vector>

using namespace std;


class OCPBoundary
{
	friend class BulkAccumuTerm01;

public:
	void Setup(const ParamReservoir& rs_param, const BulkVarSet& bvs) {
		heatLoss.Setup(rs_param, bvs, boundIndex);
		boundaryFlow.Setup(rs_param, bvs, boundName, boundIndex);
	}
	void ResetToLastTimeStep() { heatLoss.ResetToLastTimeStep(); }
	void UpdateLastTimeStep() { heatLoss.UpdateLastTimeStep(); }


public:
	/// Heat loss term
	HeatLoss     heatLoss;
	/// Boundary flow term
	BoundaryFlow boundaryFlow;

public:
	auto& GetBoundaryIndex() { return boundIndex; }
	auto& GetBoundaryArea() { return boundArea; }
	auto& GetBoundName() { return boundName; }

protected:


protected:
	vector<string>  boundName;
	vector<USI>     boundIndex;
	vector<OCP_DBL> boundArea;
};



#endif /* end if __OCPBoundary_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/22/2023      Create file                          */
/*----------------------------------------------------------------------------*/
