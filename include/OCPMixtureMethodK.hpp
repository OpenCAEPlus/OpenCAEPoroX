/*! \file    OCPMixtureMethodK.hpp
 *  \brief   OCPMixtureMethodK class declaration
 *  \author  Shizhe Li
 *  \date    Jul/31/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPMIXTUREMETHODK_HEADER__
#define __OCPMIXTUREMETHODK_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPFuncPVT.hpp"
#include "OCPMixtureVarSet.hpp"
#include "OCPPhaseEquilibrium.hpp"
#include "BulkVarSet.hpp"

#include <vector>

using namespace std;


/////////////////////////////////////////////////////
// OCPMixtureMethodK
/////////////////////////////////////////////////////


/// OCPMixtureMethodK is a bsaic class used in non-EoS model
class OCPMixtureMethodK
{
public:
    OCPMixtureMethodK() = default;
    /// Set variable set
    virtual void SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const = 0;
    /// Set variable set
    virtual void SetVarSet(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, OCPMixtureVarSet& mvs) const = 0;
    /// With P, S, Vp, perform flash calculations, and calculate VfP,Vfi only
    virtual void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    /// With P, Ni, perform flash calculations, and calculate VfP,Vfi only
    virtual void Flash(OCPMixtureVarSet& vs) = 0;
    /// With P, S, Vp, perform flash calculations, and calculate VfP,Vfi,dXsdXp
    virtual void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    /// With P, Ni, perform flash calculations, and calculate VfP,Vfi,dXsdXp
    virtual void FlashDer(OCPMixtureVarSet& vs) = 0;
    /// Flash in standard conditions
    virtual void CalVStd(OCPMixtureVarSet& vs) = 0;
    /// Calculate molar density of target phase
    virtual OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;
    /// Calculate molar density of target phase in standard conditions
    virtual OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;
    /// Calculate mass density of target phase
    virtual OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;
    /// Output iterations
    virtual void OutIters() const { }
    /// If current method is friendly to well
    virtual OCP_BOOL IfWellFriend() const = 0;
    /// Calculate Enthalpy
    virtual OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const = 0;
};


/////////////////////////////////////////////////////
// OCPMixtureMethodK_OW01
/////////////////////////////////////////////////////


/// Use PVDO and PVTW
// Note that Vo,std, Vw,std are assumed to be 1
class OCPMixtureMethodK_OW01 : public OCPMixtureMethodK
{
public:
    OCPMixtureMethodK_OW01(const ParamReservoir& rs_param, const USI& i, OCPMixtureVarSet& vs);
    void SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const override;
    void SetVarSet(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, OCPMixtureVarSet& mvs) const override;
    void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void Flash(OCPMixtureVarSet& vs) override;
    void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void FlashDer(OCPMixtureVarSet& vs) override;
    void CalVStd(OCPMixtureVarSet& vs) override;
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    OCP_BOOL IfWellFriend() const override { return OCP_TRUE; }
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override { OCP_ABORT("Not Used!"); }

protected:
    OCP_DBL CalXiO(const OCP_DBL& P) { return PVDO->CalXiO(P); }
    OCP_DBL CalXiW(const OCP_DBL& P) { return PVTW.CalXiW(P); }
    OCP_DBL CalRhoO(const OCP_DBL& P) { return PVDO->CalRhoO(P); }
    OCP_DBL CalRhoW(const OCP_DBL& P) { return PVTW.CalRhoW(P); }


protected:
    /// PVDO table
    OCP_PVDO*     PVDO;
    /// PVTW table
    OCP_PVTW      PVTW;
    /// molar volume of oil phase in standard conditions (stb/lbmol)
    const OCP_DBL stdVo{ 1 };
    /// molar volume of water phase in standard conditions (stb/lbmol)
    const OCP_DBL stdVw{ 1 };
};


/// Oil and Water are immiscible - thermal model
class OCPMixtureMethodK_OW01T : public OCPMixtureMethodK
{
public:
    OCPMixtureMethodK_OW01T(const ComponentParam& param, const USI& tarId, OCPMixtureVarSet& vs);
    void SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const override;
    void SetVarSet(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, OCPMixtureVarSet& mvs) const override;
    void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void Flash(OCPMixtureVarSet& vs) override;
    void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void FlashDer(OCPMixtureVarSet& vs) override;
    void CalVStd(OCPMixtureVarSet& vs) override;
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override { return 1 / CalXi(P, 0, T, 0, pt); }
    OCP_BOOL IfWellFriend() const override { return OCP_FALSE; }
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override { return eC.CalEnthalpy(T + CONV4, zi); }

protected:
    OCP_DBL CalXiO(const OCP_DBL& P, const OCP_DBL& T);
    OCP_DBL CalXiW(const OCP_DBL& P, const OCP_DBL& T);
    OCP_DBL CalRhoO(const OCP_DBL& P, const OCP_DBL& T);
    OCP_DBL CalRhoW(const OCP_DBL& P, const OCP_DBL& T);

protected:
    EnthalpyCalculation  eC;
    ViscosityCalculation vC;

protected:
    /// Reference pressure
    OCP_DBL Pref;   
    /// Reference temperature
    OCP_DBL Tref; 
    /// Component molar density at reference temperature and reference pressure, lb/ft3
    vector<OCP_DBL> xi_ref;
    /// Molecular Weight of components
    vector<OCP_DBL> MWc;
    /// Molecular Weight of phase
    vector<OCP_DBL> MWp;      
 
    /// Component compressibility, 1/psi
    vector<OCP_DBL> cp;              
    /// The first thermal expansion coefficient, 1/F
    vector<OCP_DBL> ct1;  
    /// The second thermal expansion coefficient, 1/F
    vector<OCP_DBL> ct2;  
    /// The coefficient of density dependence on temperature and pressure, 1/psi-F
    vector<OCP_DBL> cpt;   
             
    /// Coefficients Ak in gas viscosity correlation formulae
    vector<OCP_DBL> avg;       
    /// Coefficients Bk in gas viscosity correlation formulae
    vector<OCP_DBL> bvg;   
};


/////////////////////////////////////////////////////
// OCPMixtureMethodK_OGW01
/////////////////////////////////////////////////////


/// Use PVDO and PVTW
// Note that Vo,std, Vg,std, Vw,std are assumed to be 1
class OCPMixtureMethodK_OGW01 : public OCPMixtureMethodK
{
public:
    OCPMixtureMethodK_OGW01(const ParamReservoir& rs_param, const USI& i, OCPMixtureVarSet& vs);
    void SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const override;
    void SetVarSet(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, OCPMixtureVarSet& mvs) const override;
    void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void Flash(OCPMixtureVarSet& vs) override;
    void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void FlashDer(OCPMixtureVarSet& vs) override;
    void CalVStd(OCPMixtureVarSet& vs) override;
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    OCP_BOOL IfWellFriend() const override { return OCP_TRUE; }
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override { OCP_ABORT("Not Used!"); }

protected:
    OCP_DBL CalXiO(const OCP_DBL& P, const OCP_DBL& Pb) { return PVCO.CalXiO(P, Pb); }
    OCP_DBL CalXiG(const OCP_DBL& P) { return PVDG.CalXiG(P); }
    OCP_DBL CalXiW(const OCP_DBL& P) { return PVTW.CalXiW(P); }
    OCP_DBL CalRhoO(const OCP_DBL& P, const OCP_DBL& Pb) { return PVCO.CalRhoO(P, Pb); }
    OCP_DBL CalRhoG(const OCP_DBL& P) { return PVDG.CalRhoG(P); }
    OCP_DBL CalRhoW(const OCP_DBL& P) { return PVTW.CalRhoW(P); }

    void CalNi(const OCP_DBL& Vp, OCPMixtureVarSet& vs);

protected:
    /// PVCO table
    OCP_PVCO        PVCO;
    /// PVDG table
    OCP_PVDG        PVDG;
    /// PVTW table
    OCP_PVTW        PVTW;
    /// 1 lbmol oil components could absorb x lbmol gas components
    OCP_DBL         x;
    /// molar volume of oil phase in standard conditions
    const OCP_DBL   stdVo{ 1 };
    /// molar volume of gas phase in standard conditions
    const OCP_DBL   stdVg{ 1 };
    /// molar volume of water phase in standard conditions
    const OCP_DBL   stdVw{ 1 };
};


/////////////////////////////////////////////////////
// OCPMixtureMethodK_GW01
/////////////////////////////////////////////////////


/// Use PVTCO2 and PVTH2O
/// Here, molar density(xi) is the mass density(rho)
class OCPMixtureMethodK_GW01 : public OCPMixtureMethodK
{
public:
    OCPMixtureMethodK_GW01(const ParamReservoir& rs_param, const USI& i, OCPMixtureVarSet& vs);
    void SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const override;
    void SetVarSet(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, OCPMixtureVarSet& mvs) const override;
    void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void Flash(OCPMixtureVarSet& vs) override;
    void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void FlashDer(OCPMixtureVarSet& vs) override;
    void CalVStd(OCPMixtureVarSet& vs) override;
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    OCP_BOOL IfWellFriend() const override { return OCP_FALSE; }
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override { OCP_ABORT("Not Used!"); }

protected:
    void CalNi(const OCP_DBL& Vp, OCPMixtureVarSet& vs);
    OCP_DBL CalXiG(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z) const { return CalRhoG(P, T, z); }
    OCP_DBL CalXiW(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z) const { return CalRhoW(P, T, z); }
    OCP_DBL CalRhoG(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z) const { return PVTCO2.CalRho(P, T); }
    OCP_DBL CalRhoW(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z) const;

protected:
    /// PVT table for CO2
    OCP_PVTCO2    PVTCO2;
    /// PVT table for H2O
    OCP_PVTH2O    PVTH2O;
    /// if use water density correction
    Garciaw       garciaw;
};


#endif /* end if __OCPMIXTUREKMETHOD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/02/2023      Create file                          */
/*----------------------------------------------------------------------------*/