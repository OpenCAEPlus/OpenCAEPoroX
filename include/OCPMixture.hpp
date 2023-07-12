/*! \file    OCPMixture.hpp
 *  \brief   OCPMixture class declaration
 *  \author  Shizhe Li
 *  \date    Jul/12/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPMIXTURE_HEADER__
#define __OCPMIXTURE_HEADER__

#include "OCPConst.hpp"

#include <vector>

using namespace std;

const USI OCPMIXTURE_OG  = 21;
const USI OCPMIXTURE_OW  = 22;
const USI OCPMIXTURE_GW  = 23;
const USI OCPMIXTURE_OGW = 3;

/// mixture varset
class OCPMixtureVarSet
{
public:
    OCPMixtureVarSet() {  }
    void Init() {
    }

public:
    /// num of phase, components
    USI np, nc;
    /// pressure, temperature
    OCP_DBL P, T;
    /// mole number of components
    vector<OCP_DBL> Ni;
};


/////////////////////////////////////////////////////
// OCPMixture
/////////////////////////////////////////////////////

class OCPMixture
{
public:
    OCPMixture() = default;
    USI FlowType() const { return mixtureType; }
    OCPMixtureVarSet& GetVarSet() { return vs; }

protected:
    USI              mixtureType;
    OCPMixtureVarSet vs;
};


#endif /* end if __OCPMIXTURE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/12/2023      Create file                          */
/*----------------------------------------------------------------------------*/