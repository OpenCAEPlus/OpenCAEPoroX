/*! \file    MixtureUnitThermal.hpp
 *  \brief   MixtureUnitThermal class declaration
 *  \author  Shizhe Li
 *  \date    Nov/10/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __MIXTURETHERMAL_HEADER__
#define __MIXTURETHERMAL_HEADER__

#include <cmath>

// OpenCAEPoroX header files
#include "MixtureUnit.hpp"
#include "OCPFuncPVT.hpp"

/// MixtureUnitThermal is inherited class of Mixture, it's used for ifThermal model.
/// K-value Model
class MixtureUnitThermal : public MixtureUnit
{
public:
    MixtureUnitThermal() = default;
    void OutMixtureIters() const override{};
};

class MixtureUnitThermal_OW : public MixtureUnitThermal
{
public:
    MixtureUnitThermal_OW() = default;
    MixtureUnitThermal_OW(const ParamReservoir& param, const USI& tarId, OptionalModules& opts);
    OCPMixture* GetMixture() override { return &OWTM; }
    void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) override;
    /// flash calculation with saturation of phases.
    void InitFlashIMPEC(const OCP_USI& bId, const BulkVarSet& bvs) override;
    void InitFlashFIM(const OCP_USI& bId, const BulkVarSet& bvs) override;
    /// Flash calculation with moles of components.
    void FlashIMPEC(const OCP_USI& bId, const BulkVarSet& bvs) override;
    /// Flash calculation with moles of components and Calculate the derivative
    void FlashFIM(const OCP_USI& bId, const BulkVarSet& bvs) override;
    /// Return molar density of phase, it's used to calculate the molar density of
    /// injection fluids in injection wells.
    OCP_DBL XiPhase(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const vector<OCP_DBL>& Ziin,
                    const PhaseType& pt) override;

    /// return mass density of phase.
    OCP_DBL RhoPhase(const OCP_DBL& Pin,
                     const OCP_DBL& Pbb,
                     const OCP_DBL& Tin,
                     const vector<OCP_DBL>& Ziin,
                     const PhaseType& pt) override;

    // for Well
    OCP_DBL CalInjWellEnthalpy(const OCP_DBL& Tin, const OCP_DBL* Ziin) override;

protected:
    OCPMixtureUnitThermalOW  OWTM;
};

#endif /* end if __MIXTURETHERMAL_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           NOV/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/
