/*! \file    OptionalFeatures.hpp
 *  \brief   OptionalFeatures class declaration
 *  \author  Shizhe Li
 *  \date    Dec/25/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OPTIONALFEATURES_HEADER__
#define __OPTIONALFEATURES_HEADER__

#include "AcceleratePVT.hpp"
#include "OCPMiscible.hpp"
#include "OCPScalePcow.hpp"

class OptionalFeatures
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
    }
    void UpdateLastTimeStep()
    {
        skipPSA.UpdateLastTimeStep();
        surTen.UpdateLastTimeStep();
        misFac.UpdateLastTimeStep();
        misCur.UpdateLastTimeStep();
    }

    /////////////////////////////////////////////////////////////////////
    // Common vars
    /////////////////////////////////////////////////////////////////////
protected:
    /// num of bulks
    OCP_USI numBulk;


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
};

#endif /* end if __OptionalFeatures_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/25/2022      Create file                          */
/*----------------------------------------------------------------------------*/