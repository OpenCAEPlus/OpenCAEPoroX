/*! \file    OCPMixtureBlkoilOW.hpp
 *  \brief   OCPMixtureBlkoilOW class declaration
 *  \author  Shizhe Li
 *  \date    Jul/13/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPMIXTUREBLKOILOW_HEADER__
#define __OCPMIXTUREBLKOILOW_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPFuncPVT.hpp"
#include "OCPMixture.hpp"

#include <vector>

using namespace std;


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOWMethod
/////////////////////////////////////////////////////


class OCPMixtureBlkOilOWMethod
{
public:
    OCPMixtureBlkOilOWMethod() = default;
    virtual void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    virtual void Flash(OCPMixtureVarSet& vs) = 0;
    virtual void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    virtual void FlashDer(OCPMixtureVarSet& vs) = 0;
    virtual OCP_DBL CalXi(const OCP_DBL& P, const PhaseType& pt) = 0;
    virtual OCP_DBL CalRho(const OCP_DBL& P, const PhaseType& pt) = 0;
    virtual void CalVStd(OCPMixtureVarSet& vs) = 0;
    virtual OCP_DBL CalVmStd(const PhaseType& pt) = 0;
};


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOWMethod01
/////////////////////////////////////////////////////


/// Use PVDO and PVTW
// Note that Vo,std, Vw,std are assumed to be 1
class OCPMixtureBlkOilOWMethod01 : public OCPMixtureBlkOilOWMethod
{
public:
    OCPMixtureBlkOilOWMethod01(const ParamReservoir& rs_param, const USI& i, OCPMixtureVarSet& vs);
    void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void Flash(OCPMixtureVarSet& vs) override;
    void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void FlashDer(OCPMixtureVarSet& vs) override;
    OCP_DBL CalXi(const OCP_DBL& P, const PhaseType& pt) override;
    OCP_DBL CalRho(const OCP_DBL& P, const PhaseType& pt) override;
    OCP_DBL CalVmStd(const PhaseType& pt) override;
    void CalVStd(OCPMixtureVarSet& vs) override;
    
protected:
    OCP_DBL CalXiO(const OCP_DBL& P) { return PVDO.CalXiO(P); }
    OCP_DBL CalXiW(const OCP_DBL& P) { return PVTW.CalXiW(P); }
    OCP_DBL CalRhoO(const OCP_DBL& P) { return PVDO.CalRhoO(P); }
    OCP_DBL CalRhoW(const OCP_DBL& P) { return PVTW.CalRhoW(P); }


protected:
    /// PVDO table
    OCP_PVDO        PVDO;
    /// PVTW table
    OCP_PVTW        PVTW;
    /// molar volume of oil phase in standard conditions (stb/lbmol)
    const OCP_DBL   stdVo{ 1 };
    /// molar volume of water phase in standard conditions (stb/lbmol)
    const OCP_DBL   stdVw{ 1 };
};


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOW 
/////////////////////////////////////////////////////

class OCPMixtureBlkOilOW : public OCPMixture
{
public:
    OCPMixtureBlkOilOW() { mixtureType = OCPMixtureType::BO_OW; }
    void Setup(const ParamReservoir& rs_param, const USI& i);
    void InitFlash(const OCP_DBL& P, const OCP_DBL& Sw, const OCP_DBL& Vp) {
        SetPS(P, Sw);
        pmMethod->InitFlash(Vp, vs);
    }
    void Flash(const OCP_DBL& P, const OCP_DBL* Ni) {
        SetPN(P, Ni);
        pmMethod->Flash(vs);
    }
    void InitFlashDer(const OCP_DBL& P, const OCP_DBL& Sw, const OCP_DBL& Vp) {
        SetPS(P, Sw);
        pmMethod->InitFlashDer(Vp, vs);
    }
    void FlashDer(const OCP_DBL& P, const OCP_DBL* Ni) {
        SetPN(P, Ni);
        pmMethod->FlashDer(vs);
    }
    OCP_DBL CalXi(const OCP_DBL& P, const PhaseType& pt) {
        return pmMethod->CalXi(P, pt);
    }
    OCP_DBL CalRho(const OCP_DBL& P, const PhaseType& pt) {
        return pmMethod->CalRho(P, pt);
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
    }
    void SetPS(const OCP_DBL& P, const OCP_DBL& Sw) {
        vs.P     = P;
        vs.S[0]  = 1 -Sw;
        vs.S[1]  = Sw;      
    }

protected:
    OCPMixtureBlkOilOWMethod* pmMethod;
};


#endif /* end if __OCPMIXTUREBLKOILOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/13/2023      Create file                          */
/*----------------------------------------------------------------------------*/