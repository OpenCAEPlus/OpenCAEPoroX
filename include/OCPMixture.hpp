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

using namespace std;

/////////////////////////////////////////////////////
// OCPMixture
/////////////////////////////////////////////////////

class OCPMixture
{
public:
    OCPMixture() = default;
    auto MixtureType() const { return vs.mixtureType; }
    auto IfBlkModel() const { return (vs.mixtureType >= OCPMixtureType::SP) && (vs.mixtureType < OCPMixtureType::COMP); }
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
    void Setup(const ParamReservoir& rs_param, const USI& i);
    void Flash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) {
        SetPTN(P, T, Ni);
        pmMethod->Flash(vs);
    }
    void InitFlash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* S, const OCP_DBL* Ni, const OCP_DBL& Vp) {
        SetPTSN(P, T, S, Ni);
        pmMethod->InitFlash(Vp, vs);
    }
    void Flash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, const USI& ftype, const USI& lNP, const OCP_DBL* lx) {
        SetPTN(P, T, Ni);
        pmMethod->Flash(vs, ftype, lNP, lx);
    }
    void InitFlashDer(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* S, const OCP_DBL* Ni, const OCP_DBL& Vp) {
        SetPTSN(P, T, S, Ni);
        pmMethod->InitFlashDer(Vp, vs);
    }
    void FlashDer(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, const USI& ftype, const USI& lNP, const OCP_DBL* lx) {
        SetPTN(P, T, Ni);
        pmMethod->FlashDer(vs, ftype, lNP, lx);
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
    const auto& GetNCPE() const { return pmMethod->GetNC(); }
    const auto& GetNPPE() const { return pmMethod->GetNP(); }
    const auto& GetNPmaxPE() const { return pmMethod->GetNPmax(); }
    const auto GetEoSPE() const { return pmMethod->GetEoS(); }
    const auto GetFtypePE() const { return pmMethod->GetFtype(); }
    const auto GetNumPhasePE(const USI& np) const { return pmMethod->GetNumPhasePE(np); }
    const auto& GetZiPE() const { return pmMethod->GetZi(); }
    const auto& GetNtPE() const { return pmMethod->GetNt(); }
    const auto& GetP() const { return vs.P; }
    const auto& GetT() const { return vs.T; }

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
    OCPMixtureCompMethod* pmMethod;
};


#endif /* end if __OCPMIXTURE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/12/2023      Create file                          */
/*----------------------------------------------------------------------------*/