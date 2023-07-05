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

#include <vector>
#include <functional>

using namespace std;


/////////////////////////////////////////////////////////////////////
// Scale The Water-Oil Capillary Pressure Curves (From SWATINIT)
/////////////////////////////////////////////////////////////////////

/// Params used to scale the water permeability curve.
class ScalePcowMethodParams
{
    friend class ScalePcowMethod01;
public:
    ScalePcowMethodParams() = default;
    ScalePcowMethodParams(const OCP_DBL& Swco,
        const OCP_DBL& Pcmax,
        const OCP_DBL& Pcmin,
        const function<OCP_DBL(const OCP_DBL&)>& CalPcowf,
        OCP_DBL* Pcow,
        OCP_DBL* dPcowdSw)
    {
        swco         = Swco;
        maxPcow      = Pcmax;
        minPcow      = Pcmin;
        CalPcow      = CalPcowf;
        Pcow_out     = Pcow;
        dPcowdSw_out = dPcowdSw;
    }

protected:
    OCP_DBL                           swco;         ///< connate water saturation
    OCP_DBL                           maxPcow;      ///< maximum Pcow
    OCP_DBL                           minPcow;      ///< minimum Pcow
    function<OCP_DBL(const OCP_DBL&)> CalPcow;      ///< Functions to calculate Pcow by Sw
    OCP_DBL*                          Pcow_out;     ///< Pcow(from FlowUnit, will be scaled)
    OCP_DBL*                          dPcowdSw_out; ///< dPcowdSw(from FlowUnit, will be scaled)
};


/// Method used to scale the water permeability curve.
class ScalePcowMethod
{
public:
    /// Default constructor
    ScalePcowMethod() = default;
    /// Setup Scale coefficient
    virtual OCP_DBL SetScaleVal(const OCP_DBL& Swinit, OCP_DBL& Swinout, const OCP_DBL& Pcowin) const = 0;
    /// Scale Pcow and dPcowdSw
    virtual void ScaleDer(const OCP_DBL& sv) const = 0;
    /// Scale Pcow
    virtual void Scale(const OCP_DBL& sv) const = 0;

protected:
    OCP_DBL* Pcow_out;                   ///< Pcow(from FlowUnit, will be scaled)
    OCP_DBL* dPcowdSw_out;               ///< dPcowdSw(from FlowUnit, will be scaled)
};


/// Method used to scale the water permeability curve.
class ScalePcowMethod01 : public ScalePcowMethod
{
public:
    /// Default constructor
    ScalePcowMethod01() = default;
    /// Construct with Pmax, Pmin
    ScalePcowMethod01(const ScalePcowMethodParams& params) {
        swco         = params.swco;
        maxPcow      = params.maxPcow;
        minPcow      = params.minPcow;
        CalPcow      = params.CalPcow;
        Pcow_out     = params.Pcow_out;
        dPcowdSw_out = params.dPcowdSw_out;
    }      
    /// Setup Scale coefficient
    OCP_DBL SetScaleVal(const OCP_DBL& Swinit, OCP_DBL& Swinout, const OCP_DBL& Pcowin) const override;
    /// Scale Pcow and dPcowdSw
    void ScaleDer(const OCP_DBL& sv) const override;
    /// Scale Pcow
    void Scale(const OCP_DBL& sv) const override;

protected:
    OCP_DBL swco;                               ///< connate water saturation
    OCP_DBL maxPcow;                            ///< maximum Pcow
    OCP_DBL minPcow;                            ///< minimum Pcow
    function<OCP_DBL(const OCP_DBL&)> CalPcow;  ///< Functions to calculate Pcow by Sw
};


/// The Class scales the Pcow
class ScalePcow
{
    friend class Reservoir;
    // For Output
    friend class Out4RPT;

public:
    /// Input option param(data was distributed if have) from input file
    void InputParam(const OCP_BOOL& ifscale) { ifScale = ifscale; }
    /// Setup ScalePcow term
    USI Setup(const ScalePcowMethodParams& params);
    /// Setup Scale coefficient
    void SetScaleVal(const OCP_USI& bId, const USI& mIndex, OCP_DBL& Swinout, const OCP_DBL& Pcowin);
    /// Scale Pcow and dPcowdSw
    void ScaleDer(const OCP_USI& bId, const USI& mIndex) const;
    /// Scale Pcow
    void Scale(const OCP_USI& bId, const USI& mIndex) const;

protected:
    OCP_BOOL ifScale{ OCP_FALSE }; ///< If true, then Scale will be used
    vector<OCP_DBL> swatInit;      ///< Initial water distribution
    vector<OCP_DBL> scaleVal;      ///< Scale values for Pcow, it will be calculated from swatInit
    
    /// Method
protected:
    vector<ScalePcowMethod*> scalePcowMethod; ///< Method for scaling Pcow
};

#endif /* end if __OCPSCALEPCOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/01/2023      Create file                          */
/*----------------------------------------------------------------------------*/