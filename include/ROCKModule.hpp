/*! \file    ROCKModule.hpp
 *  \brief   ROCKModule class declaration
 *  \author  Shizhe Li
 *  \date    Aug/21/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __ROCKMODULE_HEADER__
#define __ROCKMODULE_HEADER__

 // Standard header files
#include <cassert>

// OpenCAEPoroX header files
#include "OCPRock.hpp"
#include "BulkVarSet.hpp"
#include "OptionalFeatures.hpp"

using namespace std;


class ROCKModule
{

public:
    void Setup(const ParamReservoir& rs_param, const OCP_USI& nb, OptionalFeatures& opts)
	{
	
		NTROCK = rs_param.NTROOC;

		if (rs_param.thermal) {
			for (USI i = 0; i < NTROCK; i++) {
				if (rs_param.rockSet[i].type == "LINEAR") {
					ROCKs.push_back(new OCPRockT_Linear(rs_param.rockSet[i]));
				}
				else {
					ROCKs.push_back(new OCPRockT_Exp(rs_param.rockSet[i]));
				}
			}
		}
		else {
			for (USI i = 0; i < NTROCK; i++) {
				ROCKs.push_back(new OCPRockIsoT_Linear(rs_param.rockSet[i]));
			}
		}

		if (ROCKNUM.empty()) ROCKNUM.resize(nb, 0);
	
    }
    auto GetROCK(const OCP_USI& n) const { return ROCKs[ROCKNUM[n]]; }
    auto& GetROCKNUM() { return ROCKNUM; }
    auto GetNTROCK() { return NTROCK; }

protected:
    /// num of Rock regions
    USI              NTROCK;
    /// Index of ROCK region for each bulk
    vector<USI>      ROCKNUM;
    /// Rock modules
    vector<OCPRock*> ROCKs;

};



#endif /* end if __ROCKMODULE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/21/2023      Create file                          */
/*----------------------------------------------------------------------------*/