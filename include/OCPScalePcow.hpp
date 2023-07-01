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

using namespace std;


/////////////////////////////////////////////////////////////////////
// Scale The Water-Oil Capillary Pressure Curves (From SWATINIT)
/////////////////////////////////////////////////////////////////////

class ScalePcow
{
    friend class Reservoir;
    // For Output
    friend class Out4RPT;

public:
    /// Input option param(data was distributed if have) from input file
    void InputParam(const OCP_BOOL& ifscale) { ifScale = ifscale; }
    /// Setup ScalePcow term
    void Setup();
    /// Return ifScale
    OCP_BOOL IfScale() const { return ifScale; }
    /// Assign value to scaleVal
    void    AssignScaleValue(const OCP_USI& n, const OCP_DBL& v) { scaleVal[n] = v; }
    OCP_DBL GetSwInit(const OCP_USI& n) const { return swatInit[n]; }
    /// Return scaleVal
    OCP_DBL GetScaleVal(const OCP_USI& n) const { return scaleVal[n]; }

protected:
    OCP_BOOL ifSetup{ OCP_FALSE }; ///< Only one setup is needed.
    OCP_BOOL ifScale{ OCP_FALSE }; ///< If true, then Scale will be used
    vector<OCP_DBL>
        scaleVal; ///< Scale values for Pcow, it will be calculated from swatInit
    vector<OCP_DBL> swatInit; ///< Initial water distribution
};

#endif /* end if __OCPSCALEPCOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/01/2023      Create file                          */
/*----------------------------------------------------------------------------*/