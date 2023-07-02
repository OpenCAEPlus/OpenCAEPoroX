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


class ScalePcowMethodParams
{
    friend class ScalePcowMethod01;
public:
    ScalePcowMethodParams() = default;
    ScalePcowMethodParams(const OCP_DBL& Swco,
        const OCP_DBL& Pcmax,
        const OCP_DBL& Pcmin,
        const function<OCP_DBL(const OCP_DBL&)>& CalPcowf)
    {
        swco    = Swco;
        maxPcow = Pcmax;
        minPcow = Pcmin;
        CalPcow = CalPcowf;
    }

protected:
    OCP_DBL swco;                        ///< connate water saturation
    OCP_DBL maxPcow;                     ///< maximum Pcow
    OCP_DBL minPcow;                     ///< minimum Pcow
    function<OCP_DBL(OCP_DBL)> CalPcow;  ///< Functions to calculate Pcow by Sw
};


/// Method for ScalePcow
class ScalePcowMethod
{
public:
    /// Default constructor
    ScalePcowMethod() = default;
    /// Setup Scale coefficient
    virtual OCP_DBL SetScaleVal(const OCP_DBL& Swinit, OCP_DBL& Swin, const OCP_DBL& Pcowin) const = 0;
    /// Scale Pcow and dPcowdSw
    virtual void Scale(const OCP_DBL& sv, OCP_DBL& Pcow, OCP_DBL& dPcowdSw) const = 0;
    /// Scale Pcow
    virtual void Scale(const OCP_DBL& sv, OCP_DBL& Pcow) const = 0;
};

class ScalePcowMethod01 : public ScalePcowMethod
{
public:
    /// Default constructor
    ScalePcowMethod01() = default;
    /// Construct with Pmax, Pmin
    ScalePcowMethod01(const ScalePcowMethodParams& params) {
        swco    = params.swco;
        maxPcow = params.maxPcow;
        minPcow = params.minPcow;
        CalPcow = params.CalPcow;
    }      
    /// Setup Scale coefficient
    OCP_DBL SetScaleVal(const OCP_DBL& Swinit, OCP_DBL& Swinout, const OCP_DBL& Pcowin) const override;
    /// Scale Pcow and dPcowdSw
    void Scale(const OCP_DBL& sv, OCP_DBL& Pcow, OCP_DBL& dPcowdSw) const override;
    /// Scale Pcow
    void Scale(const OCP_DBL& sv, OCP_DBL& Pcow) const override;

protected:
    OCP_DBL swco;                        ///< connate water saturation
    OCP_DBL maxPcow;                     ///< maximum Pcow
    OCP_DBL minPcow;                     ///< minimum Pcow
    function<OCP_DBL(OCP_DBL)> CalPcow;  ///< Functions to calculate Pcow by Sw
};


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
    void SetScaleVal(const USI& mIndex, const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin);
    /// Scale Pcow and dPcowdSw
    void Scale(const USI& mIndex, const OCP_USI& bId, OCP_DBL& Pcow, OCP_DBL& dPcowdSw) const;
    /// Scale Pcow
    void Scale(const USI& mIndex, const OCP_USI& bId, OCP_DBL& Pcow) const;
protected:
    OCP_BOOL ifScale{ OCP_FALSE }; ///< If true, then Scale will be used
    vector<OCP_DBL> swatInit;      ///< Initial water distribution
    vector<OCP_DBL> scaleVal;      ///< Scale values for Pcow, it will be calculated from swatInit
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