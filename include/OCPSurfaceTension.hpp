/*! \file    OCPSurfacTension.hpp
 *  \brief   OCPSurfacTension class declaration
 *  \author  Shizhe Li
 *  \date    Jul/03/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPSURFACETENSION_HEADER__
#define __OCPSURFACETENSION_HEADER__

#include "OCPConst.hpp"

#include <vector>

using namespace std;


/// Method uesd to calculate surface tension
class SurTenMethod
{
public:
    /// Default constructor
    SurTenMethod() = default;
    /// Calculate surface tensions
    virtual OCP_DBL CalSurfaceTension() const = 0;

protected:
    const USI* NP;     ///< num of present phases
};


/// Params uesd to calculate surface tension
class SurTenMethod01Params
{
    friend class SurTenMethod01;
public:
    /// Default constructor
    SurTenMethod01Params() = default;
    SurTenMethod01Params(const vector<OCP_DBL>& parachorin,
                   const OCP_DBL* xvin,
                   const OCP_DBL* xlin,
                   const OCP_DBL* bvin,
                   const OCP_DBL* blin,
                   const USI* NPin
    )
    {
        ifUse    = OCP_TRUE;
        parachor = parachorin;
        xv       = xvin;
        xl       = xlin;
        bv       = bvin;
        bl       = blin;
        NP       = NPin;
    }
    OCP_BOOL IfUse() const { return ifUse; }
protected:
    OCP_BOOL        ifUse{ OCP_FALSE };
    vector<OCP_DBL> parachor;  ///< used to calculate oil-gas surface tension by Macleod-Sugden correlation
    const OCP_DBL*  xv;        ///< molar fractions of vapour phase
    const OCP_DBL*  xl;        ///< molar fractions of liquid phase
    const OCP_DBL*  bv;        ///< molar density of vapour phase
    const OCP_DBL*  bl;        ///< molar density of liquid phase
    const USI*      NP;        ///< num of present phases
};


/// Macleod - Sugden correlation
class SurTenMethod01 : public SurTenMethod
{
public:
    /// Default constructor
    SurTenMethod01() = default;
    /// Construct with parachor
    SurTenMethod01(const SurTenMethod01Params& params)
    {
        parachor = params.parachor;
        xv = params.xv;
        xl = params.xl;
        bv = params.bv;
        bl = params.bl;
        NP = params.NP;
        NC = parachor.size();
    }
    OCP_DBL CalSurfaceTension() const;
    vector<OCP_DBL> parachor; ///< used to calculate oil-gas surface tension by Macleod-Sugden correlation
    const OCP_DBL*  xv;       ///< molar fractions of vapour phase
    const OCP_DBL*  xl;       ///< molar fractions of liquid phase
    const OCP_DBL*  bv;       ///< molar density of vapour phase
    const OCP_DBL*  bl;       ///< molar density of liquid phase
    USI             NC;       ///< num of components
};


class SurTenMethodParams
{
public:
    SurTenMethodParams(const SurTenMethod01Params& params01in) {
        params01 = params01in;
    }

public:
    SurTenMethod01Params params01;
};


class SurfaceTension
{
public:
    /// Default constructor
    SurfaceTension() = default;
    /// Setuo surface tension method
    USI Setup(const SurTenMethodParams& params);
    /// Calculate surface tension with specified method
    OCP_DBL CalSurfaceTension(const USI& mIndex) const {
        return stMethod[mIndex]->CalSurfaceTension();
    }
protected:
    vector<SurTenMethod*> stMethod;
};



#endif /* end if __PHASEPERMEABILITY_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/03/2022      Create file                          */
/*----------------------------------------------------------------------------*/