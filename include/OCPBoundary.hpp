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


#include <vector>

using namespace std;


class OCPBoundary
{

public:
	void Setup(const ParamReservoir& rs_param, const BulkVarSet& bvs) {
		heatLoss.Setup(rs_param, bvs, boundIndex);
	}
	void ResetToLastTimeStep() { heatLoss.ResetToLastTimeStep(); }
	void UpdateLastTimeStep() { heatLoss.UpdateLastTimeStep(); }


public:

	HeatLoss heatLoss;

public:
	auto& GetBoundaryIndex() { return boundIndex; }

protected:
	vector<string> boundName;
	vector<USI>    boundIndex;
};



#endif /* end if __OCPBoundary_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/22/2023      Create file                          */
/*----------------------------------------------------------------------------*/
