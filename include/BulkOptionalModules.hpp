/*! \file    BulkOptionalModules.hpp
 *  \brief   BulkOptionalModules class declaration
 *  \author  Shizhe Li
 *  \date    Dec/25/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BULKOPTIONALMODULE_HEADER__
#define __BULKOPTIONALMODULE_HEADER__

#include "AcceleratePEC.hpp"
#include "OCPMiscible.hpp"
#include "OCPScalePcow.hpp"
#include "OCPBoundary.hpp"
#include "HeatConduct.hpp"

class BulkOptionalModules
{
    friend class MixtureUnit;
    friend class OCPMixtureComp;
    friend class FlowUnit;
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
        boundary.ResetToLastTimeStep();
        heatConduct.ResetToLastTimeStep();
    }
    void UpdateLastTimeStep()
    {
        skipPSA.UpdateLastTimeStep();
        surTen.UpdateLastTimeStep();
        misFac.UpdateLastTimeStep();
        misCur.UpdateLastTimeStep();
        scalePcow.UpdateLastTimeStep();
        boundary.UpdateLastTimeStep();
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
    SkipPSA        skipPSA;

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
        boundary.Setup(rs_param, bvs);
        heatConduct.Setup(rs_param, bvs);
    }

public:

    /// Boundary condition handler
    OCPBoundary    boundary;
    /// Heat Conduct
    HeatConduct    heatConduct;

};

#endif /* end if __BulkOptionalFeatures_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/25/2022      Create file                          */
/*----------------------------------------------------------------------------*/