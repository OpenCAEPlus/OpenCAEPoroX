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


class SurTenVarSet
{
public:
    void SetNb(const OCP_USI& nbin) { nb = nbin; }
public:
    /// num of bulks
    OCP_USI         nb;
    /// Surface Tension
    vector<OCP_DBL> surTen;
};


/// Method uesd to calculate surface tension
class SurTenMethod
{
public:
    /// Default constructor
    SurTenMethod() = default;
    /// Calculate surface tensions
    virtual void CalSurfaceTension(const OCP_USI& bId, SurTenVarSet& stvs) const = 0;
};



/// Macleod - Sugden correlation
class SurTenMethod01 : public SurTenMethod
{
public:
    /// Default constructor
    SurTenMethod01(const vector<OCP_DBL>& parachorin, OCPMixtureComp* mix, SurTenVarSet& stvs) {
        parachor = parachorin;
        NC       = parachor.size();
        mixture  = mix;

        stvs.surTen.resize(stvs.nb);
    }
    void CalSurfaceTension(const OCP_USI& bId, SurTenVarSet& stvs) const;

protected:
    /// used to calculate oil-gas surface tension by Macleod-Sugden correlation
    vector<OCP_DBL> parachor; 
    /// num of components
    USI             NC;

    // Dependent modules
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
            stMethod[mIndex]->CalSurfaceTension(bId, vs);
        }        
    }
    const auto& GetVS()const { return vs; }
    const auto IfUse() const { return ifUse; }
    const auto GetSurfaceTension(const OCP_USI& bId) const {
        OCP_ASSERT(ifUse, "Surface Tension is not available!");
        return vs.surTen[bId];
    }
    void ResetTolastTimeStep() { }
    void UpdateLastTimeStep()  { }

protected:
    /// If calculate surface tension
    OCP_BOOL              ifUse{ OCP_FALSE };
    /// surface tension variable set
    SurTenVarSet          vs;
    /// surface tension calculation methods
    vector<SurTenMethod*> stMethod;
};



#endif /* end if __PHASEPERMEABILITY_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/03/2023      Create file                          */
/*----------------------------------------------------------------------------*/