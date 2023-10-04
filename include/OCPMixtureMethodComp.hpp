/*! \file    OCPMixtureMethodComp.hpp
 *  \brief   OCPMixtureMethodComp class declaration
 *  \author  Shizhe Li
 *  \date    Jul/31/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPMIXTURECOMPMETHOD_HEADER__
#define __OCPMIXTURECOMPMETHOD_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPFuncPVT.hpp"
#include "OCPMixtureVarSet.hpp"
#include "OCPPhaseEquilibrium.hpp"
#include "BulkVarSet.hpp"

#include <vector>

using namespace std;


/////////////////////////////////////////////////////
// OCPMixtureMethodComp
/////////////////////////////////////////////////////


/// OCPMixtureMethodComp is a bsaic class used in compositional model
/// all variables are about components and phases participating in 
/// Phase Equilibrium Calculations, these components and phases should be 
/// be ranked first.
class OCPMixtureMethodComp
{
public:
    OCPMixtureMethodComp() = default;
    /// Set variable set
    virtual void SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const = 0;
    /// Set variable set
    virtual void SetVarSet(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, OCPMixtureVarSet& mvs) const = 0;
    /// With P, Ni, perform flash calculations only
    virtual void Flash(OCPMixtureVarSet& vs) = 0;
    /// With P, S, Vp, perform flash calculations, and calculate VfP,Vfi only
    virtual void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    /// With P, Ni, perform flash calculations, and calculate VfP,Vfi only
    virtual void Flash(OCPMixtureVarSet& vs, const USI& ftype) = 0;
    /// With P, S, Vp, perform flash calculations, and calculate VfP,Vfi,dXsdXp
    virtual void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    /// With P, Ni, perform flash calculations, and calculate VfP,Vfi,dXsdXp
    virtual void FlashDer(OCPMixtureVarSet& vs, const USI& ftype) = 0;
    /// Calculate molar density of target phase
    virtual OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;
    /// Flash in standard conditions
    virtual void CalVStd(OCPMixtureVarSet& vs) = 0;
    /// Calculate molar density of target phase in standard conditions
    virtual OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;
    /// Calculate mass density of target phase
    virtual OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;
    /// OutPut total flash iterations during the simulation
    virtual void OutIters() const { PE.OutMixtureIters(); }
    /// if current mixture is friendly to well
    virtual OCP_BOOL IfWellFriend() const = 0;
    /// Input total number of existing phase, return the number of existing phase in PEC
    virtual USI GetNumPhasePE(const USI& np) const = 0;
    /// Get number of components in PE
    const auto& GetNC() const { return NC; }
    /// Get number if phases in PE
    const auto& GetNP() const { return NP; }
    /// Get allowable maximum number of phases in PE
    const auto& GetNPmax() const { return NPmax; }
    /// Get EoS
    const auto GetEoS() const { return &eos; }
    /// Get Ftype
    const auto& GetFtype() const { return PE.GetFtype(); }
    /// Get zi
    const auto& GetZi() const { return zi; }
    /// Get Nt
    const auto& GetNt() const { return Nt; }


////////////////////////////////////////////////////////////////
// Basic variables
////////////////////////////////////////////////////////////////

    /// Input components param from input file, allocate and setup
    void Setup(const ComponentParam& param, const USI& tarId);

protected:
    /// maximum num of allowable phases
    USI             NPmax;
    /// num of components
    USI             NC;
    /// num of existing phases
    USI             NP;
    /// total mole number of components
    OCP_DBL         Nt;
    /// molar fraction of components
    vector<OCP_DBL> zi;


////////////////////////////////////////////////////////////////
// Baisc Components Property
////////////////////////////////////////////////////////////////
    

protected:
    /// critical temperature of components
    vector<OCP_DBL> Tc;
    /// critical pressure of components
    vector<OCP_DBL> Pc;
    /// critical volume of components
    vector<OCP_DBL> Vc;
    /// Molecular Weight of components
    vector<OCP_DBL> MWC;

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
    vector<PhaseType> phaseLabel;
    /// Index of all existing phases 
    vector<USI>       epIndex;
    /// work space for reorder
    vector<OCP_DBL>   rowork;


////////////////////////////////////////////////////////////////
// Special Derivatives Calculations
////////////////////////////////////////////////////////////////

protected:
    // For Phase num : any
    /// Calculate dVf / dP, dVf / dNi (full derivatives)
    virtual void CalVfiVfp_full01(OCPMixtureVarSet& vs);
    /// Assemble Matrix for CalVfiVfp_full01
    void AssembleMatVfiVfp_full01();
    /// Assemble Rhs for CalVfiVfp_full01
    void AssembleRhsVfiVfp_full01();
    /// Calculate d (Sj xij) / d (P , Ni)
    virtual void CaldXsdXp01(OCPMixtureVarSet& vs);

    // For Phase num : <=2
    /// Calculate dVf / dP, dVf / dNi (full derivatives)
    virtual void CalVfiVfp_full02(OCPMixtureVarSet& vs);
    /// Assemble Matrix for CalVfiVfp_full02
    void AssembleMatVfiVfp_full02();
    /// Assemble Rhs for CalVfiVfp_full02
    void AssembleRhsVfiVfp_full02();
    /// Calculate d (Sj xij) / d (P , Ni)
    virtual void CaldXsdXp02(OCPMixtureVarSet& vs);

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
};


////////////////////////////////////////////////////////////////
// OCPMixtureMethodComp01
////////////////////////////////////////////////////////////////


/// water and hydrocarbon phases are immiscible
class OCPMixtureMethodComp01 : public OCPMixtureMethodComp
{
public:
    OCPMixtureMethodComp01(const ParamReservoir& rs_param, const USI& i, OCPMixtureVarSet& vs);
    void SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const override;
    void SetVarSet(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, OCPMixtureVarSet& mvs) const override;
    void Flash(OCPMixtureVarSet& vs) override;
    void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void Flash(OCPMixtureVarSet& vs, const USI& ftype) override;
    void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void FlashDer(OCPMixtureVarSet& vs, const USI& ftype) override;
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    void CalVStd(OCPMixtureVarSet& vs) override;
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    /// Water is always present
    USI GetNumPhasePE(const USI& np) const override { return np - 1; }
    OCP_BOOL IfWellFriend() const override { return OCP_FALSE; }

protected:
    void InitNtZ(OCPMixtureVarSet& vs);
    void CorrectNt(const OCP_DBL& vh, OCPMixtureVarSet& vs);

protected:

    ////////////////////////////////////////////////////////////
    // Water Property
    ////////////////////////////////////////////////////////////

protected:
    void CalPropertyW(const OCP_DBL& vw, OCPMixtureVarSet& vs);
    OCP_DBL CalXiW(const OCP_DBL& P) { return PVTW.CalXiW(P); }
    OCP_DBL CalRhoW(const OCP_DBL& P) { return PVTW.CalRhoW(P); }   

protected:
    /// PVTW table
    OCP_PVTW       PVTW;
    /// molar volume of water phase in std conditions
    const OCP_DBL  stdVw{ 1 };
    /// index of water phase
    USI            wIdP;
    /// index of water component
    USI            wIdC;


    ////////////////////////////////////////////////////////////
    // Derivatives Calculations
    ////////////////////////////////////////////////////////////

protected:
    void CalVfiVfp_full01(OCPMixtureVarSet& vs) override;
    void CaldXsdXp01(OCPMixtureVarSet& vs) override;
    void CalVfiVfp_full02(OCPMixtureVarSet& vs) override;
    void CaldXsdXp02(OCPMixtureVarSet& vs) override;
    void AddVfpVfiW(OCPMixtureVarSet& vs);
    void CaldXsdXpW(OCPMixtureVarSet& vs);
};



#endif /* end if __OCPMIXTURECOMPMETHOD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/31/2023      Create file                          */
/*----------------------------------------------------------------------------*/