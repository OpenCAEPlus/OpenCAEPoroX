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
    virtual OCP_DBL CalRhoO(const OCP_DBL& P, const OCP_DBL& Pb) = 0;
    virtual OCP_DBL CalXiO(const OCP_DBL& P, const OCP_DBL& Pb) = 0;
    virtual OCP_DBL CalRhoG(const OCP_DBL& P) = 0;
    virtual OCP_DBL CalXiG(const OCP_DBL& P) = 0;
    virtual OCP_DBL CalRhoW(const OCP_DBL& P) = 0;
    virtual OCP_DBL CalXiW(const OCP_DBL& P) = 0;
    virtual void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    virtual void Flash(OCPMixtureVarSet& vs) = 0;
    virtual void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    virtual void FlashDer(OCPMixtureVarSet& vs) = 0;
};


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOGWMethod01
/////////////////////////////////////////////////////


/// Use PVDO and PVTW
class OCPMixtureBlkOilOGWMethod01 : public OCPMixtureBlkOilOGWMethod
{
public:
    OCPMixtureBlkOilOGWMethod01(const vector<vector<OCP_DBL>>& PVCOin,
        const vector<vector<OCP_DBL>>& PVDGin,
        const vector<vector<OCP_DBL>>& PVTWin,
        const OCP_DBL& stdRhoO,
        const OCP_DBL& stdRhoG,
        const OCP_DBL& stdRhoW,
        OCPMixtureVarSet& vs);
    OCP_DBL CalRhoO(const OCP_DBL& P, const OCP_DBL& Pb) override { return PVCO.CalRhoO(P, Pb); }
    OCP_DBL CalXiO(const OCP_DBL& P, const OCP_DBL& Pb) override { return PVCO.CalXiO(P, Pb); }
    OCP_DBL CalRhoG(const OCP_DBL& P) override { return PVDG.CalRhoG(P); }
    OCP_DBL CalXiG(const OCP_DBL& P) override { return PVDG.CalXiG(P); }
    OCP_DBL CalRhoW(const OCP_DBL& P) override { return PVTW.CalRhoW(P); }
    OCP_DBL CalXiW(const OCP_DBL& P) override { return PVTW.CalXiW(P); }
    void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void Flash(OCPMixtureVarSet& vs) override;
    void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void FlashDer(OCPMixtureVarSet& vs) override;

protected:
    void CalNi(const OCP_DBL& Vp, OCPMixtureVarSet& vs);

protected:
    OCP_PVCO        PVCO;
    OCP_PVDG        PVDG;
    OCP_PVTW        PVTW;
};


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOGW 
/////////////////////////////////////////////////////

class OCPMixtureBlkOilOGW : public OCPMixture
{
public:
    OCPMixtureBlkOilOGW() { mixtureType = OCPMIXTURE_BO_OGW; }
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
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const USI& tarPhase) { 
        if (tarPhase == OIL)         return pmMethod->CalXiO(P, Pb);
        else if (tarPhase == GAS)    return pmMethod->CalXiG(P);
        else if (tarPhase == WATER)  return pmMethod->CalXiW(P);
        else                         OCP_ABORT("WRONG TarPhase");
    }
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const USI& tarPhase) {
        if (tarPhase == OIL)         return pmMethod->CalRhoO(P, Pb);
        else if (tarPhase == GAS)    return pmMethod->CalRhoG(P);
        else if (tarPhase == WATER)  return pmMethod->CalRhoW(P);
        else                         OCP_ABORT("WRONG TarPhase");
    }

protected:
    void GetStdRhoOGW(const ParamReservoir& rs_param);
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
    OCP_DBL                    stdRhoO;
    OCP_DBL                    stdRhoG;
    OCP_DBL                    stdRhoW;
};


#endif /* end if __OCPMIXTUREBLKOILOGW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/19/2023      Create file                          */
/*----------------------------------------------------------------------------*/