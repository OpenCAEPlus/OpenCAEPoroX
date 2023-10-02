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
#include "OCPMixtureKMethod.hpp"
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
    OCPMixtureComp(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts);
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


////////////////////////////////////////////////////////////////
// OCPMixtureK 
////////////////////////////////////////////////////////////////


class OCPMixtureK : public OCPMixture
{
public:
    OCPMixtureK() = default;
    OCPMixtureK(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts);
    void Flash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) {
    }
    void InitFlash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* S, const OCP_DBL* Ni, const OCP_DBL& Vp, const OCP_USI& bId) {
    }
    void Flash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, const USI& lNP, const OCP_DBL* lx, const OCP_USI& bId) {
    }
    void InitFlashDer(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* S, const OCP_DBL* Ni, const OCP_DBL& Vp, const OCP_USI& bId) {
    }
    void FlashDer(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, const OCP_DBL* S, const USI& lNP, const OCP_DBL* lx, const OCP_USI& bId) {
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
};


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOGW 
/////////////////////////////////////////////////////

class OCPMixtureBlkOilOGW : public OCPMixture
{
public:
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
        return pmMethod->CalXi(P, Pb, 0, 0, pt);
    }
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const PhaseType& pt) {
        return pmMethod->CalRho(P, Pb, 0, 0, pt);
    }
    void CalVStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPN(P, Ni);
        pmMethod->CalVStd(vs);
    }
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalVmStd(0, 0, 0, 0, pt);
    }
    OCP_BOOL IfWellFriend() const override { return pmMethod->IfWellFriend(); }

protected:
    void SetPN(const OCP_DBL& P, const OCP_DBL* Ni) {
        vs.P = P;
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
    OCPMixtureKMethod* pmMethod;
};


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOW 
/////////////////////////////////////////////////////

class OCPMixtureBlkOilOW : public OCPMixture
{
public:
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
        return pmMethod->CalXi(P, 0, 0, 0, pt);
    }
    OCP_DBL CalRho(const OCP_DBL& P, const PhaseType& pt) {
        return pmMethod->CalRho(P, 0, 0, 0, pt);
    }
    void CalVStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPN(P, Ni);
        pmMethod->CalVStd(vs);
    }
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalVmStd(0, 0, 0, 0, pt);
    }
    OCP_BOOL IfWellFriend() const override { return pmMethod->IfWellFriend(); }

protected:
    void SetPN(const OCP_DBL& P, const OCP_DBL* Ni) {
        vs.P = P;
        vs.Ni[0] = Ni[0];
        vs.Ni[1] = Ni[1];
    }
    void SetPS(const OCP_DBL& P, const OCP_DBL& Sw) {
        vs.P = P;
        vs.S[0] = 1 - Sw;
        vs.S[1] = Sw;
    }

protected:
    OCPMixtureKMethod* pmMethod;
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
        return pmMethod->CalXi(P, 0, T, 0, pt);
    }
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt) {
        return pmMethod->CalRho(P, 0, T, 0, pt);
    }
    void CalVStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPN(P, T, Ni);
        pmMethod->CalVStd(vs);
    }
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalVmStd(0, 0, 0, 0, pt);
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
    OCPMixtureKMethod* pmMethod;
};


/////////////////////////////////////////////////////
// OCPMixtureUnitThermalOW 
/////////////////////////////////////////////////////

class OCPMixtureUnitThermalOW : public OCPMixture
{
public:
    void Setup(const ParamReservoir& rs_param, const USI& i);
    void InitFlash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& Sw, const OCP_DBL& Vp) {
        SetPTS(P, T, Sw);
        pmMethod->InitFlash(Vp, vs);
    }
    void Flash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) {
        SetPTN(P, T, Ni);
        pmMethod->Flash(vs);
    }
    void InitFlashDer(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& Sw, const OCP_DBL& Vp) {
        SetPTS(P, T, Sw);
        pmMethod->InitFlashDer(Vp, vs);
    }
    void FlashDer(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) {
        SetPTN(P, T, Ni);
        pmMethod->FlashDer(vs);
    }
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt) {
        return pmMethod->CalXi(P, 0, T + CONV5, 0, pt);
    }
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt) {
        return pmMethod->CalRho(P, 0, T + CONV5, 0, pt);
    }
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalVmStd(P, 0, T + CONV5, 0, pt);
    }
    void CalVStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPTN(P, T, Ni);
        return pmMethod->CalVStd(vs);
    }

    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) { return pmMethod->CalEnthalpy(T + CONV5, zi); }
    OCP_BOOL IfWellFriend() const override { return pmMethod->IfWellFriend(); }
protected:
    void SetPTN(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) {
        vs.P = P;
        vs.T = T + CONV5;
        vs.Ni[0] = Ni[0];
        vs.Ni[1] = Ni[1];
    }
    void SetPTS(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& Sw) {
        vs.P = P;
        vs.T = T + CONV5;
        vs.S[0] = 1 - Sw;
        vs.S[1] = Sw;
    }

protected:
    OCPMixtureKMethod* pmMethod;
};


#endif /* end if __OCPMIXTURE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/12/2023      Create file                          */
/*----------------------------------------------------------------------------*/