/*! \file    MixtureUnitComp.hpp
 *  \brief   MixtureUnitComp class declaration
 *  \author  Shizhe Li
 *  \date    Nov/30/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __MIXTURECOMP_HEADER__
#define __MIXTURECOMP_HEADER__

// Standard header files
#include <algorithm>
#include <math.h>
#include <vector>

// OpenCAEPoroX header files
#include "DenseMat.hpp"
#include "UtilMath.hpp"
#include "MixtureUnit.hpp"
#include "OCPFuncPVT.hpp"
#include "OCPPhaseEquilibrium.hpp"
#include "OCPMixtureComp.hpp"

using namespace std;


class MixtureUnitComp : public MixtureUnit
{

public:
    MixtureUnitComp() = default;

    MixtureUnitComp(const ParamReservoir& rs_param, const USI& i, OptionalFeatures& opts);
    OCPMixture* GetMixture() override { return &compM; }
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

    // ftype = 0, flash from single phase
    // ftype = 1, skip phase stability analysis and num of phase = 1
    // ftype = 1, skip phase stability analysis and num of phase = 2
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


    OCP_DBL
    XiPhase(const OCP_DBL& Pin,
            const OCP_DBL& Tin,
            const vector<OCP_DBL>& Ziin,
            const PhaseType& pt) override;

    OCP_DBL
    RhoPhase(const OCP_DBL& Pin,
             const OCP_DBL& Pbb,
             const OCP_DBL& Tin,
             const vector<OCP_DBL>& Ziin,
             const PhaseType& pt) override;

    // For Well

    OCP_DBL CalInjWellEnthalpy(const OCP_DBL& Tin, const OCP_DBL* Ziin) override
    {
        OCP_ABORT("Can not be used in Compositional Model!");
    }

    void OutMixtureIters() const override { compM.OutputIters(); }


protected:
    /// mixture of components
    OCPMixtureComp  compM;
    /// Skip stability analysis
    SkipPSA*        skipPSA;
    USI             skipMethodIndex;
    /// Surface tension
    SurfaceTension* surTen;
    USI             stMethodIndex;
    /// Miscible Factor
    MiscibleFactor* misFac;
    USI             mfMethodIndex;

};


#endif //__MIXTURECOMP_HEADER__