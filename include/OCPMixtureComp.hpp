/*! \file    OCPMixtureComp.hpp
 *  \brief   OCPMixtureComp class declaration
 *  \author  Shizhe Li
 *  \date    Jul/31/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPMIXTURECOMP_HEADER__
#define __OCPMIXTURECOMP_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPFuncPVT.hpp"
#include "OCPMixture.hpp"
#include "OCPPhaseEquilibrium.hpp"


#include <vector>

using namespace std;


/////////////////////////////////////////////////////
// OCPMixtureCompMethod
/////////////////////////////////////////////////////


/// OCPMixtureCompMethod is a bsaic class used in compositional model
/// all variables are about components and phases participating in 
/// Phase Equilibrium Calculations
/// For Isothermal Model
class OCPMixtureCompMethod
{
public:
    OCPMixtureCompMethod() = default;
    virtual OCP_DBL CalRhoO(const OCP_DBL& P) = 0;
    virtual OCP_DBL CalXiO(const OCP_DBL& P) = 0;
    virtual OCP_DBL CalRhoW(const OCP_DBL& P) = 0;
    virtual OCP_DBL CalXiW(const OCP_DBL& P) = 0;
    virtual void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    virtual void Flash(OCPMixtureVarSet& vs) = 0;
    virtual void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    virtual void FlashDer(OCPMixtureVarSet& vs) = 0;

////////////////////////////////////////////////////////////////
// Basic variables
////////////////////////////////////////////////////////////////

    void Setup(const ComponentParam& param, const USI& tarId);

protected:
    /// maximum num of allowable phases
    USI             NPmax;
    /// num of components
    USI             NC;
    /// num of existing phases
    USI             NP;
    /// total mole number of components
    OCP_DBL         Nh;
    /// molar fraction of components
    vector<OCP_DBL> zi;


////////////////////////////////////////////////////////////////
// Baisc Components Property
////////////////////////////////////////////////////////////////
    

protected:
    /// Molecular Weight of components
    vector<OCP_DBL> MWC;
    /// critical temperature of components
    vector<OCP_DBL> Tc;
    /// critical pressure of components
    vector<OCP_DBL> Pc;
    /// critical volume of components
    vector<OCP_DBL> Vc;

////////////////////////////////////////////////////////////////
// EoS and Phase Equilibrium Calculation
////////////////////////////////////////////////////////////////

protected:
    /// Copy phase properties from PhaseEquilibrium
    void CopyPhaseFromPE(OCPMixtureVarSet& vs);

protected:
    /// Phase equilibrium calculations
    OCPPhaseEquilibrium PE;
    // EoS Calculations
    EoSCalculation      eos;


////////////////////////////////////////////////////////////////
// Baisc Phases Property
////////////////////////////////////////////////////////////////


protected:
    /// Calculate Molecular Weight
    void CalMW(OCPMixtureVarSet& vs);
    /// Calculate moalr volume and volume of phase
    void CalVmVj(OCPMixtureVarSet& vs);
    /// Calculate phase property
    void CalProperty(OCPMixtureVarSet& vs);
    /// Calculate phases' molar density, mass density and viscosity
    void CalXiRhoMu(OCPMixtureVarSet& vs);
    /// Calculate phase property and derivatives
    void CalPropertyDer(OCPMixtureVarSet& vs);
    /// Calculate phases' molar density, mass density and viscosity and derivatives
    void CalXiRhoMuDer(OCPMixtureVarSet& vs);

protected:

    /// Molecular Weight of phase
    vector<OCP_DBL>         MW;
    /// moalr volume of phase
    vector<OCP_DBL>         vm;
    /// d vm / dP
    vector<OCP_DBL>         vmP;
    /// d vm / dx (x is molar fraction of components)
    vector<vector<OCP_DBL>> vmx;
    // viscosity calculations
    ViscosityCalculation    visCal;


////////////////////////////////////////////////////////////////
// Phase Identification
////////////////////////////////////////////////////////////////

protected:
    /// Identify phases
    void IdentifyPhase(OCPMixtureVarSet& vs);
    /// ReOrder phases
    void ReOrderPhase(OCPMixtureVarSet& vs);
    /// ReOrder phases and ders
    void ReOrderPhaseDer(OCPMixtureVarSet& vs);

protected:
    /// phase labels, used to identify phase
    vector<USI>     phaseLabel;
    /// Index of all existing phases 
    vector<USI>     epIndex;
    /// work space for reorder
    vector<OCP_DBL> rowork;


////////////////////////////////////////////////////////////////
// Special Derivatives Calculations
////////////////////////////////////////////////////////////////

protected:
    // For Phase num : any
    /// Calculate dVf / dP, dVf / dNi (full derivatives)
    void CalVfiVfp_full01(OCPMixtureVarSet& vs);
    /// Assemble Matrix for CalVfiVfp_full01
    void AssembleMatVfiVfp_full01();
    /// Assemble Rhs for CalVfiVfp_full01
    void AssembleRhsVfiVfp_full01();
    /// Calculate d (Sj xij) / d (P , Ni)
    void CaldXsdXp01(OCPMixtureVarSet& vs);

    // For Phase num : <=2
    /// Calculate dVf / dP, dVf / dNi (full derivatives)
    void CalVfiVfp_full02(OCPMixtureVarSet& vs);
    /// Assemble Matrix for CalVfiVfp_full02
    void AssembleMatVfiVfp_full02();
    /// Assemble Rhs for CalVfiVfp_full02
    void AssembleRhsVfiVfp_full02();
    /// Calculate d (Sj xij) / d (P , Ni)
    void CaldXsdXp02(OCPMixtureVarSet& vs);

protected:
    /// d ln fij / d P
    vector<vector<OCP_DBL>> lnfugP;
    /// d ln fij / d nkj
    vector<vector<OCP_DBL>> lnfugN;
    /// Jacobian Mat for calculating derivates
    vector<OCP_DBL>         JmatDer;
    /// rhs or d nij / d Nk, d nij / dP in calVtiVtp
    /// rhs or dXs / dXp in Cal dXsdXp
    vector<OCP_DBL>         rhsDer;
    // for linearsolve with lapack
    /// used in dgesv_ in lapack
    vector<OCP_INT>         pivot;

////////////////////////////////////////////////////////////////
// Skip Stability Analysis
////////////////////////////////////////////////////////////////

protected:

    /// start point of current flash and next flash
    // INPUT:
    // ftype == 0: flash from single phase
    // ftype == 1: single phase(skip phase stability analysis)
    // ftype == 2: two phase(skip phase stability analysis)
    // OUTPUT:
    // ftype == 0: try to skip phase stability analysis in next flash for current bulk,
    //             and recalculate the range for judgement
    // ftype == 1: try to skip phase stability analysis in next flash for current bulk,
    //             and don't recalculate the range for judgement
    // ftype == 1: don't skip phase stability analysis in next flash for current bulk.
    USI             ftype{ 0 };
};


////////////////////////////////////////////////////////////////
// OCPMixtureCompMethod01
////////////////////////////////////////////////////////////////


/// water and hydrocarbon phases are immiscible
class OCPMixtureCompMethod01 : public OCPMixtureCompMethod
{
public:
    OCPMixtureCompMethod01(const ParamReservoir& rs_param, const USI& i);
    OCP_DBL CalRhoW(const OCP_DBL& P) override { return PVTW.CalRhoW(P); }
    OCP_DBL CalXiW(const OCP_DBL& P) override { return PVTW.CalXiW(P); }
    void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void Flash(OCPMixtureVarSet& vs) override;
    void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void FlashDer(OCPMixtureVarSet& vs) override;


protected:

    ////////////////////////////////////////////////////////////
    // Water Property
    ////////////////////////////////////////////////////////////

protected:
    void CalPropertyW(const OCP_DBL& vw, OCPMixtureVarSet& vs);

protected:
    /// PVTW
    OCP_PVTW PVTW;
    /// mass density of water phase in standard condition.
    OCP_DBL  std_RhoW;
    /// index of water phase
    USI      wIdP;
    /// index of water component
    USI      wIdC;
};


////////////////////////////////////////////////////////////////
// OCPMixtureComp 
////////////////////////////////////////////////////////////////

class OCPMixtureComp : public OCPMixture
{
public:
    OCPMixtureComp() { mixtureType = OCPMIXTURE_BO_OW; }
    void Setup(const ParamReservoir& rs_param, const USI& i);
    void InitFlash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* S, const OCP_DBL& Vp) {
        SetPTS(P, T, S);
        pmMethod->InitFlash(Vp, vs);
    }
    void Flash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) {
        SetPTN(P, T, Ni);
        pmMethod->Flash(vs);
    }
    void InitFlashDer(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* S, const OCP_DBL& Vp) {
        SetPTS(P, T, S);
        pmMethod->InitFlashDer(Vp, vs);
    }
    void FlashDer(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) {
        SetPTN(P, T, Ni);
        pmMethod->FlashDer(vs);
    }
    OCP_DBL CalXi(const OCP_DBL& P, const USI& tarPhase) {
        if (tarPhase == OIL)         return pmMethod->CalXiO(P);
        else if (tarPhase == WATER)  return pmMethod->CalXiW(P);
        else                         OCP_ABORT("WRONG TarPhase");
    }
    OCP_DBL CalRho(const OCP_DBL& P, const USI& tarPhase) {
        if (tarPhase == OIL)         return pmMethod->CalRhoO(P);
        else if (tarPhase == WATER)  return pmMethod->CalRhoW(P);
        else                         OCP_ABORT("WRONG TarPhase");
    }

protected:
    void SetPTS(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* S) {
        vs.P = P;
        vs.T = T + CONV5;
        copy(S, S + vs.np, vs.S.begin());
    }
    void SetPTN(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) {
        vs.P = P;
        vs.T = T + CONV5;
        copy(Ni, Ni + vs.nc, vs.Ni.begin());
    }

protected:
    OCPMixtureCompMethod* pmMethod;
};


#endif /* end if __OCPMIXTURECOMP_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/31/2023      Create file                          */
/*----------------------------------------------------------------------------*/