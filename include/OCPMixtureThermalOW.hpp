/*! \file    OCPMixtureUnitThermalOW.hpp
 *  \brief   OCPMixtureUnitThermalOW class declaration
 *  \author  Shizhe Li
 *  \date    Jul/20/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPMIXTURETHERMALOW_HEADER__
#define __OCPMIXTURETHERMALOW_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPFuncPVT.hpp"
#include "OCPMixtureVarSet.hpp"
#include "OCPMixture.hpp"

#include <vector>

using namespace std;


/////////////////////////////////////////////////////
// OCPMixtureUnitThermalOWMethod
/////////////////////////////////////////////////////


class OCPMixtureUnitThermalOWMethod
{
public:
    OCPMixtureUnitThermalOWMethod() = default;
    virtual void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    virtual void Flash(OCPMixtureVarSet& vs) = 0;
    virtual void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    virtual void FlashDer(OCPMixtureVarSet& vs) = 0;
    virtual OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt) = 0;
    virtual OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt) = 0;
    virtual void CalVStd(OCPMixtureVarSet& vs) = 0;
    virtual OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt) = 0;
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) { return eC.CalEnthalpy(T, zi); }
    virtual OCP_BOOL IfWellFriend() const = 0;

protected:
    EnthalpyCalculation  eC;
    ViscosityCalculation vC;
};


/////////////////////////////////////////////////////
// OCPMixtureUnitThermalOWMethod01
/////////////////////////////////////////////////////


/// Oil and Water are immiscible
class OCPMixtureUnitThermalOWMethod01 : public OCPMixtureUnitThermalOWMethod
{
public:
    OCPMixtureUnitThermalOWMethod01(const ComponentParam& param, const USI& tarId, OCPMixtureVarSet& vs);
    void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void Flash(OCPMixtureVarSet& vs) override;
    void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs)override;
    void FlashDer(OCPMixtureVarSet& vs) override;
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt) override;
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt) override;
    void CalVStd(OCPMixtureVarSet& vs) override;
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt) override { return 1 / CalXi(P, T, pt); }
    OCP_BOOL IfWellFriend() const override { return OCP_FALSE; }

protected:
    OCP_DBL CalXiO(const OCP_DBL& P, const OCP_DBL& T);
    OCP_DBL CalXiW(const OCP_DBL& P, const OCP_DBL& T);
    OCP_DBL CalRhoO(const OCP_DBL& P, const OCP_DBL& T);
    OCP_DBL CalRhoW(const OCP_DBL& P, const OCP_DBL& T);

protected:
    /// Reference pressure
    OCP_DBL Pref{PRESSURE_STD};   
    /// Reference temperature
    OCP_DBL Tref{TEMPERATURE_STD}; 
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
        return pmMethod->CalXi(P, T + CONV5, pt);
    }
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& T, const PhaseType& pt) {
        return pmMethod->CalRho(P, T + CONV5, pt);
    }
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override {
        return pmMethod->CalVmStd(P, T + CONV5, pt);
    }
    void CalVStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) override {
        SetPTN(P, T, Ni);
        return pmMethod->CalVStd(vs);
    }
    
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) { return pmMethod->CalEnthalpy(T + CONV5, zi); }
    OCP_BOOL IfWellFriend() const override { return pmMethod->IfWellFriend(); }
protected:
    void SetPTN(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) {
        vs.P     = P;
        vs.T     = T + CONV5;
        vs.Ni[0] = Ni[0];
        vs.Ni[1] = Ni[1];
    }
    void SetPTS(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& Sw) {
        vs.P    = P;
        vs.T    = T + CONV5;
        vs.S[0] = 1 - Sw;
        vs.S[1] = Sw;
    }

protected:
    OCPMixtureUnitThermalOWMethod* pmMethod;
};


#endif /* end if __OCPMIXTURETHERMALOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/20/2023      Create file                          */
/*----------------------------------------------------------------------------*/