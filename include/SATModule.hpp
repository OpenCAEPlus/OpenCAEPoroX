/*! \file    SATModule.hpp
 *  \brief   SATModule class declaration
 *  \author  Shizhe Li
 *  \date    Aug/21/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __SATMODULE_HEADER__
#define __SATMODULE_HEADER__

 // Standard header files
#include <cassert>

// OpenCAEPoroX header files
#include "FlowUnit.hpp"
#include "BulkVarSet.hpp"

using namespace std;


class SATModule
{

public:
    void Setup(const ParamReservoir& rs_param, const BulkVarSet& bvs, OptionalModules& opts)
    {
        NTSFUN = rs_param.NTSFUN;

        // Setup Saturation function
        for (USI i = 0; i < NTSFUN; i++)
            SATs.push_back(new FlowUnit(rs_param, i, opts));


        if (SATNUM.empty() || NTSFUN == 1) {
            SATNUM.clear();
            SATNUM.resize(bvs.nb, 0);
        }
    }
    auto GetSAT(const OCP_USI& n) const { return SATs[SATNUM[n]]; }
    auto& GetSATNUM() { return SATNUM; }
    auto GetNTSFUN() { return NTSFUN; }

protected:
    /// number of SAT regions
    USI               NTSFUN;
    /// Index of SAT region for each bulk
    vector<USI>       SATNUM;
    /// SAT modules
    vector<FlowUnit*> SATs;
};



#endif /* end if __SATMODULE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/21/2023      Create file                          */
/*----------------------------------------------------------------------------*/