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


#include "OCPMixtureCompMethod.hpp"
#include "OptionalModules.hpp"

using namespace std;

/////////////////////////////////////////////////////
// OCPMixture
/////////////////////////////////////////////////////

class OCPMixture
{
public:
    OCPMixture() = default;
    auto MixtureType() const { return vs.mixtureType; }
    virtual OCP_BOOL IfWellFriend() const = 0;
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


////////////////////////////////////////////////////////////////
// OCPMixtureComp 
////////////////////////////////////////////////////////////////

class OCPMixtureComp : public OCPMixture
{
public:
    void Setup(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts);
    void Flash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) {
        SetPTN(P, T, Ni);
        pmMethod->Flash(vs);
    }
    void InitFlash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* S, const OCP_DBL* Ni, const OCP_DBL& Vp, const OCP_USI& bId) {
        SetPTSN(P, T, S, Ni);
        pmMethod->InitFlash(Vp, vs);
        skipPSA->CalSkipForNextStep(bId, skipMethodIndex);
    }
    void Flash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, const USI& lNP, const OCP_DBL* lx, const OCP_USI& bId) {    
        SetPTN(P, T, Ni);
        const USI ftype = skipPSA->CalFtype01(bId, skipMethodIndex, vs);
        pmMethod->Flash(vs, ftype, lNP, lx);
        skipPSA->CalSkipForNextStep(bId, skipMethodIndex);
    }
    void InitFlashDer(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* S, const OCP_DBL* Ni, const OCP_DBL& Vp, const OCP_USI& bId) {
        SetPTSN(P, T, S, Ni);
        pmMethod->InitFlashDer(Vp, vs);
        skipPSA->CalSkipForNextStep(bId, skipMethodIndex);
    }
    void FlashDer(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, const OCP_DBL* S, const USI& lNP, const OCP_DBL* lx, const OCP_USI& bId) {      
        SetPTSN(P, T, S, Ni);
        const USI ftype = skipPSA->CalFtype02(bId, skipMethodIndex, vs, lNP);
        pmMethod->FlashDer(vs, ftype, lNP, lx);
        skipPSA->CalSkipForNextStep(bId, skipMethodIndex);
    }
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) {
        return pmMethod->CalXi(P, T + CONV5, z, pt);
    }
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) {
        return pmMethod->CalRho(P, T + CONV5, z, pt);
    }
    void CalVStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPTN(P, T, Ni);
        return pmMethod->CalVStd(vs);
    }
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalVmStd(P, T + CONV5, z, pt);
    }
    void OutputIters() const { pmMethod->OutIters(); }
    OCP_BOOL IfWellFriend() const override { return pmMethod->IfWellFriend(); }

protected:
    void SetPTSN(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* S, const OCP_DBL* Ni) {
        vs.P = P;
        vs.T = T + CONV5;
        copy(S, S + vs.np, vs.S.begin());
        copy(Ni, Ni + vs.nc, vs.Ni.begin());
    }
    void SetPTN(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) {
        vs.P = P;
        vs.T = T + CONV5;
        copy(Ni, Ni + vs.nc, vs.Ni.begin());
    }

protected:
    /// method
    OCPMixtureCompMethod* pmMethod;
    // dependent module
    /// Skip stability analysis
    SkipPSA*              skipPSA;
    USI                   skipMethodIndex;
};


#endif /* end if __OCPMIXTURE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/12/2023      Create file                          */
/*----------------------------------------------------------------------------*/