/*! \file    OCPScalePcow.hpp
 *  \brief   OCPScalePcow class declaration
 *  \author  Shizhe Li
 *  \date    Jul/01/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPSCALEPCOW__HEADER__
#define __OCPSCALEPCOW__HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPFlowOGW.hpp"

#include <vector>

using namespace std;


/////////////////////////////////////////////////////////////////////
// Scale The Water-Oil Capillary Pressure Curves (From SWATINIT)
/////////////////////////////////////////////////////////////////////


class ScalePcowVarSet
{
public:
    void SetNb(const OCP_USI& nbin) { nb = nbin; }
public:
    /// number of bulks
    OCP_USI         nb;
    /// initial water distribution
    vector<OCP_DBL> swatInit;
    /// Scale values for Pcow, it will be calculated from swatInit      
    vector<OCP_DBL> scaleVal;      
};


/// Method used to scale the water permeability curve.
class ScalePcowMethod
{
public:
    /// Default constructor
    ScalePcowMethod() = default;
    /// Setup Scale coefficient
    virtual void SetScaleVal(const OCP_USI& bId, ScalePcowVarSet& spvs, OCP_DBL& Swinout, const OCP_DBL& Pcowin) const = 0;
    /// Scale Pcow and dPcowdSw
    virtual void ScaleDer(const OCP_USI& bId, const ScalePcowVarSet& spvs) const = 0;
    /// Scale Pcow
    virtual void Scale(const OCP_USI& bId, const ScalePcowVarSet& spvs) const = 0;
};


/// Method used to scale the water permeability curve.
class ScalePcowMethod01 : public ScalePcowMethod
{
public:
    /// Construct with Pmax, Pmin
    ScalePcowMethod01(OCPFlow* flowin, ScalePcowVarSet& spvs);
    /// Setup Scale coefficient
    void SetScaleVal(const OCP_USI& bId, ScalePcowVarSet& spvs, OCP_DBL& Swinout, const OCP_DBL& Pcowin) const override;
    /// Scale Pcow and dPcowdSw
    void ScaleDer(const OCP_USI& bId, const ScalePcowVarSet& spvs) const override;
    /// Scale Pcow
    void Scale(const OCP_USI& bId, const ScalePcowVarSet& spvs) const override;

protected:
    /// saturation of connate water
    OCP_DBL        Swco;
    /// maximum capillary pressure: Po - Pw
    OCP_DBL        maxPcow;
    /// minimum capillary pressure: Po - Pw
    OCP_DBL        minPcow;

    // Dependent modules
    /// flow model
    OCPFlow*       flow;
};


/// The Class scales the Pcow
class ScalePcow
{
    // For Output
    friend class Out4RPT;

public:
    /// Setup ScalePcow term
    USI Setup(OCPFlow* flowin);
    /// Setup Scale coefficient
    void SetScaleVal(const OCP_USI& bId, const USI& mIndex, OCP_DBL& Swinout, const OCP_DBL& Pcowin);
    /// Scale Pcow and dPcowdSw
    void ScaleDer(const OCP_USI& bId, const USI& mIndex) const;
    /// Scale Pcow
    void Scale(const OCP_USI& bId, const USI& mIndex) const;
    /// Get swatinit
    auto& GetSwatInit() { return vs.swatInit; }
    /// Reset to last time step
    void ResetTolastTimeStep() { }
    /// Update last time step
    void UpdateLastTimeStep() { }

protected:
    /// If scale
    OCP_BOOL ifScale{ OCP_FALSE };
    /// ScalePcow variable set
    ScalePcowVarSet vs;
    /// Method for scaling Pcow
    vector<ScalePcowMethod*> scalePcowMethod;
};

#endif /* end if __OCPSCALEPCOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/01/2023      Create file                          */
/*----------------------------------------------------------------------------*/