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
#include "OCPMixtureThermalOW.hpp"

/// MixtureUnitThermal is inherited class of Mixture, it's used for ifThermal model.
/// K-value Model
class MixtureUnitThermal : public MixtureUnit
{
public:
    MixtureUnitThermal() = default;
    void Allocate()
    {
        Ni.resize(numCom);
        phaseExist.resize(numPhase);
        S.resize(numPhase);
        vj.resize(numPhase);
        nj.resize(numPhase);
        x.resize(numPhase * numCom);
        rho.resize(numPhase);
        xi.resize(numPhase);
        mu.resize(numPhase);
        vfi.resize(numCom);
        rhoP.resize(numPhase);
        rhoT.resize(numPhase);
        rhox.resize(numPhase * numCom);
        xiP.resize(numPhase);
        xiT.resize(numPhase);
        xix.resize(numPhase * numCom);
        muP.resize(numPhase);
        muT.resize(numPhase);
        mux.resize(numPhase * numCom);
        dXsdXp.resize((numCom + 1) * numPhase * (numCom + 2));
        Ufi.resize(numCom);
        H.resize(numPhase);
        HT.resize(numPhase);
        Hx.resize(numPhase * numCom);
    }
    void SetupOptionalFeatures(OptionalFeatures& optFeatures) override{};
    void OutMixtureIters() const override{};
};

class MixtureUnitThermal_OW : public MixtureUnitThermal
{
public:
    MixtureUnitThermal_OW() = default;
    MixtureUnitThermal_OW(const ParamReservoir& param, const USI& tarId);
    void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) override;
    /// flash calculation with saturation of phases.
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
    /// Flash calculation with moles of components.
    void FlashIMPEC(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const OCP_DBL* Niin,
                    const USI&     lastNP,
                    const OCP_DBL* xijin,
                    const OCP_USI& bId) override;
    /// Flash calculation with moles of components and Calculate the derivative
    void FlashFIM(const OCP_DBL& Pin,
                  const OCP_DBL& Tin,
                  const OCP_DBL* Niin,
                  const OCP_DBL* Sjin,
                  const USI&     lastNP,
                  const OCP_DBL* xijin,
                  const OCP_USI& bId) override;
    /// Return molar density of phase, it's used to calculate the molar density of
    /// injection fluids in injection wells.
    OCP_DBL XiPhase(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const vector<OCP_DBL>& Ziin,
                    const USI&     tarPhase) override;

    /// return mass density of phase.
    OCP_DBL RhoPhase(const OCP_DBL& Pin,
                     const OCP_DBL& Pbb,
                     const OCP_DBL& Tin,
                     const vector<OCP_DBL>& Ziin,
                     const USI&     tarPhase) override;

    // for Well
    void CalProdWeight(const OCP_DBL&         Pin,
                       const OCP_DBL&         Tin,
                       const OCP_DBL*         Niin,
                       const vector<OCP_DBL>& prodPhase,
                       vector<OCP_DBL>&       prodWeight) override;

    void CalProdRate(const OCP_DBL&   Pin,
                     const OCP_DBL&   Tin,
                     const OCP_DBL*   Niin,
                     vector<OCP_DBL>& prodRate) override;

    void    SetupWellOpt(WellOpt&                  wellopt,
                         const vector<SolventINJ>& sols,
                         const OCP_DBL&            Psurf,
                         const OCP_DBL&            Tsurf) override;
    OCP_DBL CalInjWellEnthalpy(const OCP_DBL& Tin, const OCP_DBL* Ziin) override;

    const OCP_DBL& GetNt() const override { return OWTM.GetVarSet().Nt; }
    const OCP_DBL& GetNi(const USI& i) const override {
        return OWTM.GetVarSet().Ni[i];
    }
    const OCP_DBL& GetVf() const override { return OWTM.GetVarSet().Vf; }
    const OCP_BOOL& GetPhaseExist(const USI& j) const override { return OWTM.GetVarSet().phaseExist[j]; }
    const OCP_DBL& GetS(const USI& j) const override { return OWTM.GetVarSet().S[j]; }
    const OCP_DBL& GetVj(const USI& j) const override { return OWTM.GetVarSet().vj[j]; }
    const OCP_DBL& GetNj(const USI& j) const override { return OWTM.GetVarSet().nj[j]; }
    const OCP_DBL& GetXij(const USI& j, const USI& i) const override
    {
        return OWTM.GetVarSet().xij[j * numCom + i];
    }
    const OCP_DBL& GetRho(const USI& j) const override { return OWTM.GetVarSet().rho[j]; }
    const OCP_DBL& GetXi(const USI& j) const override { return OWTM.GetVarSet().xi[j]; }
    const OCP_DBL& GetMu(const USI& j) const override { return OWTM.GetVarSet().mu[j]; }
    const OCP_DBL& GetVfP() const override { return OWTM.GetVarSet().vfP; }
    const OCP_DBL& GetVfT() const override { return OWTM.GetVarSet().vfT; }
    const OCP_DBL& GetVfi(const USI& i) const override { return OWTM.GetVarSet().vfi[i]; }
    const OCP_DBL& GetRhoP(const USI& j) const override { return OWTM.GetVarSet().rhoP[j]; }
    const OCP_DBL& GetRhoT(const USI& j) const override { return OWTM.GetVarSet().rhoT[j]; }
    const OCP_DBL& GetXiP(const USI& j) const override { return OWTM.GetVarSet().xiP[j]; }
    const OCP_DBL& GetXiT(const USI& j) const override { return OWTM.GetVarSet().xiT[j]; }
    const OCP_DBL& GetMuP(const USI& j) const override { return OWTM.GetVarSet().muP[j]; }
    const OCP_DBL& GetMuT(const USI& j) const override { return OWTM.GetVarSet().muT[j]; }
    const OCP_DBL& GetRhoX(const USI& j, const USI& i) const override
    {
        return OWTM.GetVarSet().rhox[j * numCom + i];
    }
    const OCP_DBL& GetXiX(const USI& j, const USI& i) const override
    {
        return OWTM.GetVarSet().xix[j * numCom + i];
    }
    const OCP_DBL& GetMuX(const USI& j, const USI& i) const override
    {
        return OWTM.GetVarSet().mux[j * numCom + i];
    }
    const vector<OCP_DBL>& GetDXsDXp() const override { return OWTM.GetVarSet().dXsdXp; }
    const OCP_DBL          GetUf() const override { return OWTM.GetVarSet().Uf; }
    const OCP_DBL          GetUfP() const override { return OWTM.GetVarSet().UfP; }
    const OCP_DBL          GetUfT() const override { return OWTM.GetVarSet().UfT; }
    const OCP_DBL          GetUfi(const USI& i) const override { return OWTM.GetVarSet().Ufi[i]; }
    const OCP_DBL          GetH(const USI& j) const override { return OWTM.GetVarSet().H[j]; }
    const OCP_DBL          GetHT(const USI& j) const override { return OWTM.GetVarSet().HT[j]; }
    const OCP_DBL& GetHx(const USI& j, const USI& i) const override
    {
        return OWTM.GetVarSet().Hx[j * numCom + i];
    }

protected:
    OCPMixtureUnitThermalOW  OWTM;
    EnthalpyCalculation  eC;
};

#endif /* end if __MIXTURETHERMAL_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           NOV/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/
