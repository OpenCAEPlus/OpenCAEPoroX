/*! \file    OCPMixtureBlkOilGW.hpp
 *  \brief   OCPMixtureBlkOilGW class declaration
 *  \author  Shizhe Li
 *  \date    Sep/22/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPMIXTUREBLKOILGW_HEADER__
#define __OCPMIXTUREBLKOILGW_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPFuncPVT.hpp"
#include "OCPMixture.hpp"

#include <vector>

using namespace std;


/////////////////////////////////////////////////////
// OCPMixtureBlkOilGWMethod
/////////////////////////////////////////////////////


class OCPMixtureBlkOilGWMethod
{
public:
    OCPMixtureBlkOilGWMethod() = default;
    virtual void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    virtual void Flash(OCPMixtureVarSet& vs) = 0;
    virtual void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    virtual void FlashDer(OCPMixtureVarSet& vs) = 0;
    virtual OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt) = 0;
    virtual OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt) = 0;
    virtual void CalVStd(OCPMixtureVarSet& vs) = 0;
    virtual OCP_DBL CalVmStd(const PhaseType& pt) = 0;
    virtual OCP_BOOL IfWellFriend() const = 0;
};


/////////////////////////////////////////////////////
// OCPMixtureBlkOilGWMethod01
/////////////////////////////////////////////////////


/// Use PVTCO2 and PVTH2O
/// Here, molar density(xi) is the mass density(rho)
class OCPMixtureBlkOilGWMethod01 : public OCPMixtureBlkOilGWMethod
{
public:
    OCPMixtureBlkOilGWMethod01(const ParamReservoir& rs_param, const USI& i, OCPMixtureVarSet& vs);
    void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void Flash(OCPMixtureVarSet& vs) override;
    void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void FlashDer(OCPMixtureVarSet& vs) override;
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt) override;
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt) override;
    OCP_DBL CalVmStd(const PhaseType& pt) override;
    void CalVStd(OCPMixtureVarSet& vs) override;
    OCP_BOOL IfWellFriend() const override { return OCP_FALSE; }

protected:
    void CalNi(const OCP_DBL& Vp, OCPMixtureVarSet& vs);
    OCP_DBL CalXiG(const OCP_DBL& P, const OCP_DBL& T) const { return CalRhoG(P, T); }
    OCP_DBL CalXiW(const OCP_DBL& P, const OCP_DBL& T) const { return CalRhoW(P, T); }
    OCP_DBL CalRhoG(const OCP_DBL& P, const OCP_DBL& T) const { return PVTCO2.CalRho(P, T); }
    OCP_DBL CalRhoW(const OCP_DBL& P, const OCP_DBL& T) const;


protected:
    /// PVT table for CO2
    OCP_PVTCO2    PVTCO2;
    /// PVT table for H2O
    OCP_PVTH2O    PVTH2O;
    /// volume per unit mass of gas in standard condition
    OCP_DBL       stdVg;
    /// volume per unit mass of water in standard condition
    OCP_DBL       stdVw;
    /// if use water density correction
    Garciaw       garciaw;
};


/////////////////////////////////////////////////////
// OCPMixtureBlkOilGW 
/////////////////////////////////////////////////////

class OCPMixtureBlkOilGW : public OCPMixture
{
public:
    void Setup(const ParamReservoir& rs_param, const USI& i);
    void InitFlash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& Sw, const OCP_DBL& Vp) {
        SetPS(P, T, Sw);
        pmMethod->InitFlash(Vp, vs);
    }
    void Flash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) {
        SetPN(P, T, Ni);
        pmMethod->Flash(vs);
    }
    void InitFlashDer(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& Sw, const OCP_DBL& Vp) {
        SetPS(P, T, Sw);
        pmMethod->InitFlashDer(Vp, vs);
    }
    void FlashDer(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) {
        SetPN(P, T, Ni);
        pmMethod->FlashDer(vs);
    }
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt) {
        return pmMethod->CalXi(P, T, pt);
    }
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt) {
        return pmMethod->CalRho(P, T, pt);
    }
    void CalVStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPN(P, T, Ni);
        pmMethod->CalVStd(vs);
    }
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalVmStd(pt);
    }
    OCP_BOOL IfWellFriend() const override { return pmMethod->IfWellFriend(); }

protected:
    void SetPN(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) {
        vs.P = P;
        vs.T = T;
        vs.Ni[0] = Ni[0];
        vs.Ni[1] = Ni[1];
    }
    void SetPS(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& Sw) {
        vs.P = P;
        vs.T = T;
        vs.S[0] = 1 - Sw;
        vs.S[1] = Sw;
    }

protected:
    OCPMixtureBlkOilGWMethod* pmMethod;
};


#endif /* end if __OCPMixtureBlkOilGW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/22/2023      Create file                          */
/*----------------------------------------------------------------------------*/