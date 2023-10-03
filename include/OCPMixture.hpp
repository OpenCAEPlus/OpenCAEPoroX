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
    const OCPMixtureVarSet& GetVarSet() const { return vs; }
    virtual void Flash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) = 0;
    virtual void InitFlash(const OCP_USI& bId, const BulkVarSet& bvs) = 0;
    virtual void Flash(const OCP_USI& bId, const BulkVarSet& bvs) = 0;
    virtual void InitFlashDer(const OCP_USI& bId, const BulkVarSet& bvs) = 0;
    virtual void FlashDer(const OCP_USI& bId, const BulkVarSet& bvs) = 0;
    virtual void CalVStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) = 0;
    virtual OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;
    virtual OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;
    virtual OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;
    virtual void OutputIters() const = 0;
    virtual OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const = 0;
    virtual OCP_BOOL IfWellFriend() const = 0;

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
    void Flash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPTN(P, T, Ni);
        pmMethod->Flash(vs);
    }
    void InitFlash(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->InitFlash(bvs.rockVp[bId], vs);
        skipPSA->CalSkipForNextStep(bId, skipMethodIndex, vs);
    }
    void Flash(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        const USI ftype = skipPSA->CalFtype01(bId, skipMethodIndex, vs);
        pmMethod->Flash(vs, ftype);
        skipPSA->CalSkipForNextStep(bId, skipMethodIndex, vs);
    }
    void InitFlashDer(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->InitFlashDer(bvs.rockVp[bId], vs);
        skipPSA->CalSkipForNextStep(bId, skipMethodIndex, vs);
    }
    void FlashDer(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        const USI ftype = skipPSA->CalFtype02(bId, skipMethodIndex, vs);
        pmMethod->FlashDer(vs, ftype);
        skipPSA->CalSkipForNextStep(bId, skipMethodIndex, vs);
    }
    void CalVStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPTN(P, T, Ni);
        return pmMethod->CalVStd(vs);
    }
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalVmStd(P, T + CONV5, z, pt);
    }
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalXi(P, T + CONV5, z, pt);
    }
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalRho(P, T + CONV5, z, pt);
    }
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override { OCP_ABORT("Not Used!"); }
    void OutputIters() const override { pmMethod->OutIters(); }
    OCP_BOOL IfWellFriend() const override { return pmMethod->IfWellFriend(); }

protected:
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


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOGW 
/////////////////////////////////////////////////////

class OCPMixtureBlkOilOGW : public OCPMixture
{
public:
    void Setup(const ParamReservoir& rs_param, const USI& i);
    OCPMixtureBlkOilOGW() = default;
    OCPMixtureBlkOilOGW(const ParamReservoir& rs_param, const USI& i);
    void Flash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPN(P, Ni);
        pmMethod->Flash(vs);
    }
    void InitFlash(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->InitFlash(bvs.rockVp[bId], vs);
    }
    void Flash(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->Flash(vs);
    }
    void InitFlashDer(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->InitFlashDer(bvs.rockVp[bId], vs);
    }
    void FlashDer(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->FlashDer(vs);
    }
    void CalVStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPN(P, Ni);
        pmMethod->CalVStd(vs);
    }
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalVmStd(0, 0, 0, 0, pt);
    }
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalXi(P, Pb, 0, 0, pt);
    }
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalRho(P, Pb, 0, 0, pt);
    }
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override { OCP_ABORT("Not Used!"); }
    void OutputIters() const override { OCP_ABORT("Not Used!"); }
    OCP_BOOL IfWellFriend() const override { return pmMethod->IfWellFriend(); }

protected:
    void SetPN(const OCP_DBL& P, const OCP_DBL* Ni) {
        vs.P = P;
        vs.Ni[0] = Ni[0];
        vs.Ni[1] = Ni[1];
        vs.Ni[2] = Ni[2];
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
    OCPMixtureBlkOilOW() = default;
    OCPMixtureBlkOilOW(const ParamReservoir& rs_param, const USI& i);
    void Flash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPN(P, Ni);
        pmMethod->Flash(vs);
    }
    void InitFlash(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->InitFlash(bvs.rockVp[bId], vs);
    }
    void Flash(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->Flash(vs);
    }
    void InitFlashDer(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->InitFlashDer(bvs.rockVp[bId], vs);
    }
    void FlashDer(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->FlashDer(vs);
    }
    void CalVStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPN(P, Ni);
        pmMethod->CalVStd(vs);
    }
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalVmStd(0, 0, 0, 0, pt);
    }
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalXi(P, 0, 0, 0, pt);
    }
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalRho(P, 0, 0, 0, pt);
    }
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override { OCP_ABORT("Not Used!"); }
    void OutputIters() const override { OCP_ABORT("Not Used!"); }
    OCP_BOOL IfWellFriend() const override { return pmMethod->IfWellFriend(); }

protected:
    void SetPN(const OCP_DBL& P, const OCP_DBL* Ni) {
        vs.P = P;
        vs.Ni[0] = Ni[0];
        vs.Ni[1] = Ni[1];
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
    OCPMixtureBlkOilGW() = default;
    OCPMixtureBlkOilGW(const ParamReservoir& rs_param, const USI& i);
    void Flash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPN(P, T, Ni);
        pmMethod->Flash(vs);
    }
    void InitFlash(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->InitFlash(bvs.rockVp[bId], vs);
    }
    void Flash(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->Flash(vs);
    }
    void InitFlashDer(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->InitFlashDer(bvs.rockVp[bId], vs);
    }
    void FlashDer(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->FlashDer(vs);
    }
    void CalVStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPN(P, T, Ni);
        pmMethod->CalVStd(vs);
    }
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalVmStd(0, 0, 0, 0, pt);
    }
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalXi(P, 0, T, 0, pt);
    }
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalRho(P, 0, T, 0, pt);
    }
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override { OCP_ABORT("Not Used!"); }
    void OutputIters() const override { OCP_ABORT("Not Used!"); }
    OCP_BOOL IfWellFriend() const override { return pmMethod->IfWellFriend(); }

protected:
    void SetPN(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) {
        vs.P = P;
        vs.T = T;
        vs.Ni[0] = Ni[0];
        vs.Ni[1] = Ni[1];
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
    OCPMixtureUnitThermalOW() = default;
    OCPMixtureUnitThermalOW(const ParamReservoir& rs_param, const USI& i);
    void Flash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPTN(P, T, Ni);
        pmMethod->Flash(vs);
    }
    void Flash(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->Flash(vs);
    }
    void InitFlash(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->InitFlash(bvs.rockVp[bId], vs);
    }
    void InitFlashDer(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->InitFlashDer(bvs.rockVp[bId], vs);
    }
    void FlashDer(const OCP_USI& bId, const BulkVarSet& bvs) override {
        pmMethod->SetVarSet(bId, bvs, vs);
        pmMethod->FlashDer(vs);
    }
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalVmStd(P, 0, T + CONV5, 0, pt);
    }
    void CalVStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPTN(P, T, Ni);
        return pmMethod->CalVStd(vs);
    }
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalXi(P, 0, T + CONV5, 0, pt);
    }
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalRho(P, 0, T + CONV5, 0, pt);
    }
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override { return pmMethod->CalEnthalpy(T + CONV5, zi); }
    void OutputIters() const override { OCP_ABORT("Not Used!"); }
    OCP_BOOL IfWellFriend() const override { return pmMethod->IfWellFriend(); }
protected:
    void SetPTN(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) {
        vs.P = P;
        vs.T = T + CONV5;
        vs.Ni[0] = Ni[0];
        vs.Ni[1] = Ni[1];
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