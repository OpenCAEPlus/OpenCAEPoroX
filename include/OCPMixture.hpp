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

#include "OCPMixtureVarSet.hpp"

using namespace std;

/////////////////////////////////////////////////////
// OCPMixture
/////////////////////////////////////////////////////

class OCPMixture
{
public:
    OCPMixture() = default;
    auto MixtureType() const { return vs.mixtureType; }
    auto IfBlkModel() const { return (vs.mixtureType >= OCPMixtureType::SP) && (vs.mixtureType < OCPMixtureType::COMP); }
    const OCPMixtureVarSet& GetVarSet() const { return vs; }
    virtual void CalVStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) = 0;
    virtual OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;

public:
    auto OilIndex() const { return vs.o; }
    auto GasIndex() const { return vs.g; }
    auto WatIndex() const { return vs.w; }
    auto LiquidIndex() const { return vs.l; }

protected:
    /// mixture variables set
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