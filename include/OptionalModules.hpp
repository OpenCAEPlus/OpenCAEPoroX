/*! \file    OptionalModules.hpp
 *  \brief   OptionalModules class declaration
 *  \author  Shizhe Li
 *  \date    Dec/25/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OPTIONALMODULE_HEADER__
#define __OPTIONALMODULE_HEADER__

#include "AcceleratePEC.hpp"
#include "OCPMiscible.hpp"
#include "OCPScalePcow.hpp"
#include "HeatLoss.hpp"
#include "HeatConduct.hpp"

class OptionalModules
{
    friend class MixtureUnitComp;
    friend class FlowUnit_OGW;
    friend class FlowUnit_OW;
    friend class Reservoir;

    // For Output
    friend class Out4RPT;

public:
    void ResetToLastTimeStep()
    {
        skipPSA.ResetToLastTimeStep();
        surTen.ResetTolastTimeStep();
        misFac.ResetTolastTimeStep();
        misCur.ResetTolastTimeStep();
        scalePcow.ResetTolastTimeStep();
        heatLoss.ResetToLastTimeStep();
        heatConduct.ResetToLastTimeStep();
    }
    void UpdateLastTimeStep()
    {
        skipPSA.UpdateLastTimeStep();
        surTen.UpdateLastTimeStep();
        misFac.UpdateLastTimeStep();
        misCur.UpdateLastTimeStep();
        scalePcow.UpdateLastTimeStep();
        heatLoss.UpdateLastTimeStep();
        heatConduct.UpdateLastTimeStep();
    }

    /////////////////////////////////////////////////////////////////////
    // Common vars
    /////////////////////////////////////////////////////////////////////
protected:
    /// num of bulks
    OCP_USI nb;

/////////////////////////////////////////////////////////////////////
// Attachment Module
/////////////////////////////////////////////////////////////////////


    /////////////////////////////////////////////////////////////////////
    // Accelerate PVT
    /////////////////////////////////////////////////////////////////////

protected:
    /// Skip Stability Analysis
    SkipPSA skipPSA;

    /////////////////////////////////////////////////////////////////////
    // Phase Permeability Curve
    /////////////////////////////////////////////////////////////////////

protected:
    /// Surface Tension Calculation
    SurfaceTension surTen;
    /// Miscible Factore Calculation
    MiscibleFactor misFac;
    /// Miscible Curve Correct
    MiscibleCurve  misCur;
    /// Scale water-oil capillary pressure
    ScalePcow      scalePcow;



/////////////////////////////////////////////////////////////////////
// Independent Module
/////////////////////////////////////////////////////////////////////

public:

    void SetupIndependentModule(const ParamReservoir& rs_param, const BulkVarSet& bvs) {
        heatLoss.Setup(rs_param, bvs.nb);
        heatConduct.Setup(rs_param, bvs);
    }

    /////////////////////////////////////////////////////////////////////
    // Thermal
    /////////////////////////////////////////////////////////////////////

public:

    /// Heat Loss
    HeatLoss       heatLoss;
    /// Heat Conduct
    HeatConduct    heatConduct;

};

#endif /* end if __OptionalFeatures_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/25/2022      Create file                          */
/*----------------------------------------------------------------------------*/