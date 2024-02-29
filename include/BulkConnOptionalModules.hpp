/*! \file    BulkConnOptionalModules.hpp
 *  \brief   BulkConnOptionalModules class declaration
 *  \author  Shizhe Li
 *  \date    Feb/28/2024
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BULKCONNOPTIONALMODULE_HEADER__
#define __BULKCONNOPTIONALMODULE_HEADER__

#include "HeatConduct.hpp"
#include "OCPDiffusion.hpp"

class BulkConnOptionalModules
{

public:
    void ResetToLastTimeStep()
    {
        heatConduct.ResetToLastTimeStep();
    }
    void UpdateLastTimeStep()
    {
        heatConduct.UpdateLastTimeStep();
    }

public:

    void Setup(const ParamReservoir& rs_param, const Bulk& bk) {
        heatConduct.Setup(rs_param, bk.GetVarSet());
        diffusion.Setup(rs_param, bk.GetVarSet());
    }

public:

    /// Heat Conduct
    HeatConduct    heatConduct;
    /// Mass diffusion
    OCPDiffusion   diffusion;

};

#endif /* end if __BulkConnOptionalFeatures_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Feb/28/2024      Create file                          */
/*----------------------------------------------------------------------------*/