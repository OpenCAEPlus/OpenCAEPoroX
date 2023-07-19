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


/// Calculate oil, gas, water relative permeability and capillary pressure
class OCPMixtureBlkOilOWMethod
{
public:
    OCPMixtureBlkOilOWMethod() = default;
    virtual OCP_DBL CalRho(const OCP_DBL& P, const USI& tarPhase) = 0;
    virtual OCP_DBL CalXi(const OCP_DBL& P, const USI& tarPhase) = 0;
    virtual void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    virtual void Flash(OCPMixtureVarSet& vs) = 0;
    virtual void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    virtual void FlashDer(OCPMixtureVarSet& vs) = 0;
};


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOWMethod01
/////////////////////////////////////////////////////


/// Use PVDO and PVTW
class OCPMixtureBlkOilOWMethod01 : public OCPMixtureBlkOilOWMethod
{
public:
    OCPMixtureBlkOilOWMethod01(const vector<vector<OCP_DBL>>& PVDOin,
        const vector<vector<OCP_DBL>>& PVTWin,
        const OCP_DBL& stdRhoO,
        const OCP_DBL& stdRhoW,
        OCPMixtureVarSet& vs);
    OCP_DBL CalRho(const OCP_DBL& P, const USI& tarPhase) override;
    OCP_DBL CalXi(const OCP_DBL& P, const USI& tarPhase) override;
    void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void Flash(OCPMixtureVarSet& vs) override;
    void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void FlashDer(OCPMixtureVarSet& vs) override;
    
protected:
    OCP_PVDO        PVDO;
    OCP_PVTW        PVTW;
};


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOW 
/////////////////////////////////////////////////////

class OCPMixtureBlkOilOW : public OCPMixture
{
public:
    OCPMixtureBlkOilOW() { mixtureType = OCPMIXTURE_BO_OW; }
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
    OCP_DBL CalXi(const OCP_DBL& P, const USI& tarPhase) { return pmMethod->CalXi(P, tarPhase); }
    OCP_DBL CalRho(const OCP_DBL& P, const USI& tarPhase) { return pmMethod->CalRho(P, tarPhase); }

protected:
    void GetStdRhoOW(const ParamReservoir& rs_param);
    void SetPN(const OCP_DBL& P, const OCP_DBL* Ni) { 
        vs.P     = P;
        vs.Ni[0] = Ni[0];
        vs.Ni[1] = Ni[1];
    }
    void SetPS(const OCP_DBL& P, const OCP_DBL& Sw) {
        vs.P     = P;
        vs.S[1]  = Sw;
    }

protected:
    ///< Phase and Component index
    OCPMixtureBlkOilOWMethod* pmMethod;
    OCP_DBL                   stdRhoO;
    OCP_DBL                   stdRhoW;
};


#endif /* end if __OCPMIXTUREBLKOILOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/13/2023      Create file                          */
/*----------------------------------------------------------------------------*/