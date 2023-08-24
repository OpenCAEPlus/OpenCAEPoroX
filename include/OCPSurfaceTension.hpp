/*! \file    OCPSurfacTension.hpp
 *  \brief   OCPSurfacTension class declaration
 *  \author  Shizhe Li
 *  \date    Jul/03/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPSURFACETENSION_HEADER__
#define __OCPSURFACETENSION_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPMixtureComp.hpp"

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
};



/// Macleod - Sugden correlation
class SurTenMethod01 : public SurTenMethod
{
public:
    /// Default constructor
    SurTenMethod01(const vector<OCP_DBL>& parachorin, OCPMixtureComp* mix) {
        parachor = parachorin;
        NC       = parachor.size();
        mixture  = mix;
    }
    OCP_DBL CalSurfaceTension() const;

protected:
    /// used to calculate oil-gas surface tension by Macleod-Sugden correlation
    vector<OCP_DBL> parachor; 
    /// num of components
    USI             NC;
    /// mixture model
    OCPMixtureComp* mixture;
};


class SurfaceTension
{
public:
    /// Default constructor
    SurfaceTension() = default;
    /// Setup surface tension method(different mixture use different setup)
    USI Setup(const ParamReservoir& rs_param, const USI& i, const OCP_USI& nb, OCPMixtureComp* mix);
    /// Calculate surface tension with specified method
    void CalSurfaceTension(const OCP_USI& bId, const USI& mIndex) {
        if (ifUse) {
            surTen[bId] = stMethod[mIndex]->CalSurfaceTension();
        }        
    }
    const auto IfUse() const { return ifUse; }
    const auto GetSurfaceTension(const OCP_USI& bId) const {
        OCP_ASSERT(ifUse, "Surface Tension is not available!");
        return surTen[bId];
    }
    void ResetTolastTimeStep() { }
    void UpdateLastTimeStep()  { }

protected:
    /// If calculate surface tension
    OCP_BOOL              ifUse{ OCP_FALSE };
    /// surface tension calculation methods
    vector<SurTenMethod*> stMethod;
    /// surface tension
    vector<OCP_DBL>       surTen;
};



#endif /* end if __PHASEPERMEABILITY_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/03/2023      Create file                          */
/*----------------------------------------------------------------------------*/