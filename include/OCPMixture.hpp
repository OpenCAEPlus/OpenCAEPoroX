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


#include "OCPMixtureMethodComp.hpp"
#include "OCPMixtureMethodK.hpp"
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
        pmMethod->SetVarSet(P, T, Ni, vs);
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
        pmMethod->SetVarSet(P, T, Ni, vs);
        return pmMethod->CalVStd(vs);
    }
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalVmStd(P, T, z, pt);
    }
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalXi(P, T, z, pt);
    }
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalRho(P, T, z, pt);
    }
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override { OCP_ABORT("Not Used!"); }
    void OutputIters() const override { pmMethod->OutIters(); }
    OCP_BOOL IfWellFriend() const override { return pmMethod->IfWellFriend(); }

protected:
    /// method
    OCPMixtureMethodComp* pmMethod;
    // dependent module
    /// Skip stability analysis
    SkipPSA*              skipPSA;
    USI                   skipMethodIndex;
};


/////////////////////////////////////////////////////
// OCPMixtureK
/////////////////////////////////////////////////////

class OCPMixtureK : public OCPMixture
{
public:
    OCPMixtureK() = default;
    OCPMixtureK(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts);
    void Flash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        pmMethod->SetVarSet(P, T, Ni, vs);
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
        pmMethod->SetVarSet(P, T, Ni, vs);
        pmMethod->CalVStd(vs);
    }
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalVmStd(P, P, T, z, pt);
    }
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalXi(P, Pb, T, z, pt);
    }
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalRho(P, Pb, T, z, pt);
    }
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override { return pmMethod->CalEnthalpy(T, zi); }
    void OutputIters() const override { return pmMethod->OutIters(); }
    OCP_BOOL IfWellFriend() const override { return pmMethod->IfWellFriend(); }


protected:
    OCPMixtureMethodK* pmMethod;
};



#endif /* end if __OCPMIXTURE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/12/2023      Create file                          */
/*----------------------------------------------------------------------------*/