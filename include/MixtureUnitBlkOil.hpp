/*! \file    MixtureUnitBlkOil.hpp
 *  \brief   MixtureUnitBlkOil class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __MIXTUREBO_HEADER__
#define __MIXTUREBO_HEADER__

#include <cmath>

// OpenCAEPoroX header files
#include "MixtureUnit.hpp"
#include "OCPFuncPVT.hpp"
#include "OCPMixtureBlkOilOW.hpp"
#include "OCPMixtureBlkOilGW.hpp"
#include "OCPMixtureBlkOilOGW.hpp"

/// MixtureUnitBlkOil is inherited class of Mixture, it's used for black oil model.
class MixtureUnitBlkOil : public MixtureUnit
{
public:
    MixtureUnitBlkOil() = default;

    // For Well
    OCP_DBL CalInjWellEnthalpy(const OCP_DBL& Tin, const OCP_DBL* Ziin) override
    {
        OCP_ABORT("Can not be used in Black Oil Model!");
    }
    void OutMixtureIters() const override{};

};

///////////////////////////////////////////////
// MixtureUnitBlkOil_SP
///////////////////////////////////////////////

// Single Phase
class MixtureUnitBlkOil_SP : public MixtureUnitBlkOil
{
public:
    MixtureUnitBlkOil_SP() = default;
    MixtureUnitBlkOil_SP(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts)
    {
        OCP_ABORT("Not Completed!");
    };
    void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) override
    {
        OCP_ABORT("Not Completed!");
    }
    void InitFlashIMPEC(const OCP_DBL& Pin,
                        const OCP_DBL& Pbbin,
                        const OCP_DBL& Tin,
                        const OCP_DBL* Sjin,
                        const OCP_DBL& Vpore,
                        const OCP_DBL* Ziin,
                        const OCP_USI& bId) override
    {
        OCP_ABORT("Not Completed!");
    };
    void InitFlashFIM(const OCP_DBL& Pin,
                      const OCP_DBL& Pbbin,
                      const OCP_DBL& Tin,
                      const OCP_DBL* Sjin,
                      const OCP_DBL& Vpore,
                      const OCP_DBL* Ziin,
                      const OCP_USI& bId) override
    {
        OCP_ABORT("Not Completed!");
    };
    void FlashIMPEC(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const OCP_DBL* Niin,
                    const USI&     lastNP,
                    const OCP_DBL* xijin,
                    const OCP_USI& bId) override
    {
        OCP_ABORT("Not Completed!");
    };
    void FlashFIM(const OCP_DBL& Pin,
                  const OCP_DBL& Tin,
                  const OCP_DBL* Niin,
                  const OCP_DBL* Sjin,
                  const USI&     lastNP,
                  const OCP_DBL* xijin,
                  const OCP_USI& bId) override
    {
        OCP_ABORT("Not Completed!");
    };
    OCP_DBL XiPhase(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const vector<OCP_DBL>& Ziin,
                    const PhaseType& pt) override
    {
        OCP_ABORT("Not Completed!");
        return 0;
    };
    OCP_DBL RhoPhase(const OCP_DBL& Pin,
                     const OCP_DBL& Pbb,
                     const OCP_DBL& Tin,
                     const vector<OCP_DBL>& Ziin,
                     const PhaseType& pt) override
    {
        OCP_ABORT("Not Completed!");
        return 0;
    };
};

///////////////////////////////////////////////
// MixtureUnitBlkOil_OW
///////////////////////////////////////////////

class MixtureUnitBlkOil_OW : public MixtureUnitBlkOil
{
public:
    MixtureUnitBlkOil_OW() = default;
    MixtureUnitBlkOil_OW(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts);
    OCPMixture* GetMixture() override { return &OWM; }
    void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) override;
    void InitFlashIMPEC(const OCP_DBL& Pin,
                        const OCP_DBL& Pbbin,
                        const OCP_DBL& Tin,
                        const OCP_DBL* Sjin,
                        const OCP_DBL& Vpore,
                        const OCP_DBL* Ziin,
                        const OCP_USI& bId) override;
    void InitFlashFIM(const OCP_DBL& Pin,
                      const OCP_DBL& Pbbin,
                      const OCP_DBL& Tin,
                      const OCP_DBL* Sjin,
                      const OCP_DBL& Vpore,
                      const OCP_DBL* Ziin,
                      const OCP_USI& bId) override;
    void FlashIMPEC(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const OCP_DBL* Niin,
                    const USI&     lastNP,
                    const OCP_DBL* xijin,
                    const OCP_USI& bId) override;
    void FlashFIM(const OCP_DBL& Pin,
                  const OCP_DBL& Tin,
                  const OCP_DBL* Niin,
                  const OCP_DBL* Sjin,
                  const USI&     lastNP,
                  const OCP_DBL* xijin,
                  const OCP_USI& bId) override;
    OCP_DBL XiPhase(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const vector<OCP_DBL>& Ziin,
                    const PhaseType& pt) override;
    OCP_DBL RhoPhase(const OCP_DBL& Pin,
                     const OCP_DBL& Pbb,
                     const OCP_DBL& Tin,
                     const vector<OCP_DBL>& Ziin,
                     const PhaseType& pt) override;

protected:
    OCPMixtureBlkOilOW OWM;
};


///////////////////////////////////////////////
// MixtureUnitBlkOil_GW
///////////////////////////////////////////////

class MixtureUnitBlkOil_GW : public MixtureUnitBlkOil
{
public:
    MixtureUnitBlkOil_GW() = default;
    MixtureUnitBlkOil_GW(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts);
    OCPMixture* GetMixture() override { return &GWM; }
    void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) override;
    void InitFlashIMPEC(const OCP_DBL& Pin,
        const OCP_DBL& Pbbin,
        const OCP_DBL& Tin,
        const OCP_DBL* Sjin,
        const OCP_DBL& Vpore,
        const OCP_DBL* Ziin,
        const OCP_USI& bId) override;
    void InitFlashFIM(const OCP_DBL& Pin,
        const OCP_DBL& Pbbin,
        const OCP_DBL& Tin,
        const OCP_DBL* Sjin,
        const OCP_DBL& Vpore,
        const OCP_DBL* Ziin,
        const OCP_USI& bId) override;
    void FlashIMPEC(const OCP_DBL& Pin,
        const OCP_DBL& Tin,
        const OCP_DBL* Niin,
        const USI& lastNP,
        const OCP_DBL* xijin,
        const OCP_USI& bId) override;
    void FlashFIM(const OCP_DBL& Pin,
        const OCP_DBL& Tin,
        const OCP_DBL* Niin,
        const OCP_DBL* Sjin,
        const USI& lastNP,
        const OCP_DBL* xijin,
        const OCP_USI& bId) override;
    OCP_DBL XiPhase(const OCP_DBL& Pin,
        const OCP_DBL& Tin,
        const vector<OCP_DBL>& Ziin,
        const PhaseType& pt) override;
    OCP_DBL RhoPhase(const OCP_DBL& Pin,
        const OCP_DBL& Pbb,
        const OCP_DBL& Tin,
        const vector<OCP_DBL>& Ziin,
        const PhaseType& pt) override;

private:
    OCPMixtureBlkOilGW GWM;
};



///////////////////////////////////////////////
// MixtureUnitBlkOil_OGW
///////////////////////////////////////////////

class MixtureUnitBlkOil_OGW : public MixtureUnitBlkOil
{
public:
    MixtureUnitBlkOil_OGW() = default;
    MixtureUnitBlkOil_OGW(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts);
    OCPMixture* GetMixture() override { return &OGWM; }
    void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) override;
    void InitFlashIMPEC(const OCP_DBL& Pin,
                        const OCP_DBL& Pbbin,
                        const OCP_DBL& Tin,
                        const OCP_DBL* Sjin,
                        const OCP_DBL& Vpore,
                        const OCP_DBL* Ziin,
                        const OCP_USI& bId) override;
    void InitFlashFIM(const OCP_DBL& Pin,
                      const OCP_DBL& Pbbin,
                      const OCP_DBL& Tin,
                      const OCP_DBL* Sjin,
                      const OCP_DBL& Vpore,
                      const OCP_DBL* Ziin,
                      const OCP_USI& bId) override;
    void FlashIMPEC(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const OCP_DBL* Niin,
                    const USI&     lastNP,
                    const OCP_DBL* xijin,
                    const OCP_USI& bId) override;
    void FlashFIM(const OCP_DBL& Pin,
                  const OCP_DBL& Tin,
                  const OCP_DBL* Niin,
                  const OCP_DBL* Sjin,
                  const USI&     lastNP,
                  const OCP_DBL* xijin,
                  const OCP_USI& bId) override;
    OCP_DBL XiPhase(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const vector<OCP_DBL>& Ziin,
                    const PhaseType& pt) override;
    OCP_DBL RhoPhase(const OCP_DBL& Pin,
                     const OCP_DBL& Pbb,
                     const OCP_DBL& Tin,
                     const vector<OCP_DBL>& Ziin,
                     const PhaseType& pt) override;

protected:
    OCPMixtureBlkOilOGW OGWM;
};

#endif /* end if __MIXTUREBO_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/