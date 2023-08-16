/*! \file    OCPMixtureBlkOilOGW.hpp
 *  \brief   OCPMixtureBlkOilOGW class declaration
 *  \author  Shizhe Li
 *  \date    Jul/19/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPMIXTUREBLKOILOGW_HEADER__
#define __OCPMIXTUREBLKOILOGW_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPFuncPVT.hpp"
#include "OCPMixture.hpp"

#include <vector>

using namespace std;


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOGWMethod
/////////////////////////////////////////////////////


class OCPMixtureBlkOilOGWMethod
{
public:
    OCPMixtureBlkOilOGWMethod() = default;
    virtual OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const PhaseType& pt) = 0;
    virtual OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const PhaseType& pt) = 0;
    virtual void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    virtual void Flash(OCPMixtureVarSet& vs) = 0;
    virtual void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    virtual void FlashDer(OCPMixtureVarSet& vs) = 0;
    virtual OCP_DBL CalVmStd(const PhaseType& pt) = 0;
    virtual void CalVStd(OCPMixtureVarSet& vs) = 0;
};


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOGWMethod01
/////////////////////////////////////////////////////


/// Use PVDO and PVTW
// Note that Vo,std, Vg,std, Vw,std are assumed to be 1
class OCPMixtureBlkOilOGWMethod01 : public OCPMixtureBlkOilOGWMethod
{
public:
    OCPMixtureBlkOilOGWMethod01(const ParamReservoir& rs_param, const USI& i,
        OCPMixtureVarSet& vs);
    void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void Flash(OCPMixtureVarSet& vs) override;
    void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void FlashDer(OCPMixtureVarSet& vs) override;
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const PhaseType& pt) override;
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const PhaseType& pt) override;
    OCP_DBL CalVmStd(const PhaseType& pt) override;
    void CalVStd(OCPMixtureVarSet& vs) override;

protected:
    OCP_DBL CalXiO(const OCP_DBL& P, const OCP_DBL& Pb) { return PVCO.CalXiO(P, Pb); }
    OCP_DBL CalXiG(const OCP_DBL& P) { return PVDG.CalXiG(P); }
    OCP_DBL CalXiW(const OCP_DBL& P) { return PVTW.CalXiW(P); }
    OCP_DBL CalRhoO(const OCP_DBL& P, const OCP_DBL& Pb) { return PVCO.CalRhoO(P, Pb); }
    OCP_DBL CalRhoG(const OCP_DBL& P) { return PVDG.CalRhoG(P); }
    OCP_DBL CalRhoW(const OCP_DBL& P) { return PVTW.CalRhoW(P); }

    void CalNi(const OCP_DBL& Vp, OCPMixtureVarSet& vs);

protected:
    /// PVCO table
    OCP_PVCO        PVCO;
    /// PVDG table
    OCP_PVDG        PVDG;
    /// PVTW table
    OCP_PVTW        PVTW;
    /// 1 lbmol oil components could absorb x lbmol gas components
    OCP_DBL         x;
    /// molar volume of oil phase in standard conditions
    const OCP_DBL   stdVo{ 1 };
    /// molar volume of gas phase in standard conditions
    const OCP_DBL   stdVg{ 1 };
    /// molar volume of water phase in standard conditions
    const OCP_DBL   stdVw{ 1 };
};


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOGW 
/////////////////////////////////////////////////////

class OCPMixtureBlkOilOGW : public OCPMixture
{
public:
    OCPMixtureBlkOilOGW() { mixtureType = OCPMixtureType::BO_OGW; }
    void Setup(const ParamReservoir& rs_param, const USI& i);
    void InitFlash(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& Sg, const OCP_DBL& Sw, const OCP_DBL& Vp) {
        SetPS(P, Pb, Sg, Sw);
        pmMethod->InitFlash(Vp, vs);
    }
    void Flash(const OCP_DBL& P, const OCP_DBL* Ni) {
        SetPN(P, Ni);
        pmMethod->Flash(vs);
    }
    void InitFlashDer(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& Sg, const OCP_DBL& Sw, const OCP_DBL& Vp) {
        SetPS(P, Pb, Sg, Sw);
        pmMethod->InitFlashDer(Vp, vs);
    }
    void FlashDer(const OCP_DBL& P, const OCP_DBL* Ni) {
        SetPN(P, Ni);
        pmMethod->FlashDer(vs);
    }
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const PhaseType& pt) {
        return pmMethod->CalXi(P, Pb, pt);
    }
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const PhaseType& pt) {
        return pmMethod->CalRho(P, Pb, pt);
    }
    void CalVStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPN(P, Ni);
        pmMethod->CalVStd(vs);
    }
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalVmStd(pt);
    }

protected:
    void SetPN(const OCP_DBL& P, const OCP_DBL* Ni) {
        vs.P     = P;
        vs.Ni[0] = Ni[0];
        vs.Ni[1] = Ni[1];
        vs.Ni[2] = Ni[2];
    }
    void SetPS(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& Sg, const OCP_DBL& Sw) {
        vs.P    = P;
        vs.Pb   = Pb;
        vs.S[0] = 1 - Sg - Sw;
        vs.S[1] = Sg;
        vs.S[2] = Sw;
    }

protected:
    OCPMixtureBlkOilOGWMethod* pmMethod;
};


#endif /* end if __OCPMIXTUREBLKOILOGW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/19/2023      Create file                          */
/*----------------------------------------------------------------------------*/