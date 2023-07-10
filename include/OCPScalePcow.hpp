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
};


/// Method used to scale the water permeability curve.
class ScalePcowMethod01 : public ScalePcowMethod
{
public:
    /// Default constructor
    ScalePcowMethod01() = default;
    /// Construct with Pmax, Pmin
    ScalePcowMethod01(OCPFlow* flowin);
    /// Setup Scale coefficient
    OCP_DBL SetScaleVal(const OCP_DBL& Swinit, OCP_DBL& Swinout, const OCP_DBL& Pcowin) const override;
    /// Scale Pcow and dPcowdSw
    void ScaleDer(const OCP_DBL& sv) const override;
    /// Scale Pcow
    void Scale(const OCP_DBL& sv) const override;

protected:
    OCP_DBL        Swco;
    OCP_DBL        maxPcow;
    OCP_DBL        minPcow;
    OCPFlow*       flow;
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
    USI Setup(OCPFlow* flowin);
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