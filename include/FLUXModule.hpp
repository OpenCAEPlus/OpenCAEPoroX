/*! \file    FLUXModule.hpp
 *  \brief   FLUXModule class declaration
 *  \author  Shizhe Li
 *  \date    Feb/25/2024
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __FLUXMODULE_HEADER__
#define __FLUXMODULE_HEADER__

 // Standard header files
#include <cassert>

// OpenCAEPoroX header files
#include "FluxUnit.hpp"

using namespace std;


class FLUXModule
{

public:
    void Setup(const ParamReservoir& rs_param, const vector<BulkConnPair>& iterConn, const Bulk& bk, BulkConnOptionalModules& optMs)
    {
        const USI np = bk.GetPhaseNum();
        const USI nc = bk.GetComNum();
        if (rs_param.thermal) {
            FLUXs.push_back(FluxUnit(2, np, nc, optMs));
            FLUXNUM.resize(iterConn.size(), 0);
        }
        else {
            if (rs_param.GRAVDR) {
                FLUXs.push_back(FluxUnit(0, np, nc, optMs));
                FLUXs.push_back(FluxUnit(1, np, nc, optMs));
                FLUXNUM.resize(iterConn.size(), 0);
                for (OCP_USI c = 0; c < iterConn.size(); c++) {
					if (iterConn[c].Direction() == ConnDirect::mf ||
                        iterConn[c].Direction() == ConnDirect::fm) {
                        FLUXNUM[c] = 1;
					}
					else {
                        FLUXNUM[c] = 0;
					}
                }
            }
            else {
                FLUXs.push_back(FluxUnit(0, np, nc, optMs));
                FLUXNUM.resize(iterConn.size(), 0);
            }
        }
    }
    auto GetFlux(const OCP_USI& c) const { return &FLUXs[FLUXNUM[c]]; }

protected:
    /// number of PVT regions
    USI                  NTFLUX;
    /// Index of PVT region for each bulk
    vector<USI>          FLUXNUM;
    /// PVT modules
    vector<FluxUnit>     FLUXs;
};



#endif /* end if __FLUXMODULE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Feb/25/2024      Create file                          */
/*----------------------------------------------------------------------------*/