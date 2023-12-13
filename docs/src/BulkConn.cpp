/*! \file    BulkConn.cpp
 *  \brief   BulkConn class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// Standard header files
#include <cassert>
#include <cmath>
#include <ctime>

// OpenCAEPoroX header files
#include "BulkConn.hpp"

/////////////////////////////////////////////////////////////////////
// General
/////////////////////////////////////////////////////////////////////


void BulkConn::InputParam(const ParamReservoir& rs_param, const Bulk& bk)
{
    if (rs_param.thermal) {
        flux.push_back(new OCPFluxT01(bk.vs.np, bk.vs.nc));
        for (auto& iter : iteratorConn) {
            iter.SetFluxNum(0);
            flux[iter.FluxNum()]->CalBulkConnArea(iter, bk);
        }
    }
    else {
        if (rs_param.GRAVDR) {
            flux.push_back(new OCPFlux01(bk.vs.np, bk.vs.nc));
            flux.push_back(new OCPFlux02(bk.vs.np, bk.vs.nc));
            for (auto& iter : iteratorConn) {
                if (iter.Direction() == ConnDirect::mf ||
                    iter.Direction() == ConnDirect::fm) {
                    iter.SetFluxNum(1);
                }
                else {
                    iter.SetFluxNum(0);
                }             
                flux[iter.FluxNum()]->CalBulkConnArea(iter, bk);
            }
        }
        else {
            flux.push_back(new OCPFlux01(bk.vs.np, bk.vs.nc));
            for (auto& iter : iteratorConn) {
                iter.SetFluxNum(0);
                flux[iter.FluxNum()]->CalBulkConnArea(iter, bk);
            }
        }
    }

    
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/