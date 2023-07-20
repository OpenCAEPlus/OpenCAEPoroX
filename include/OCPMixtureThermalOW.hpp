/*! \file    OCPMixtureThermalOW.hpp
 *  \brief   OCPMixtureThermalOW class declaration
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
#include "OCPMixture.hpp"

#include <vector>

using namespace std;


/////////////////////////////////////////////////////
// OCPMixtureThermalOWMethod
/////////////////////////////////////////////////////


class OCPMixtureThermalOWMethod
{
public:
    OCPMixtureThermalOWMethod() = default;
    virtual OCP_DBL CalRhoO(const OCP_DBL& P, const OCP_DBL& T) = 0;
    virtual OCP_DBL CalXiO(const OCP_DBL& P, const OCP_DBL& T) = 0;
    virtual OCP_DBL CalRhoW(const OCP_DBL& P, const OCP_DBL& T) = 0;
    virtual OCP_DBL CalXiW(const OCP_DBL& P, const OCP_DBL& T) = 0;
    virtual void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    virtual void Flash(OCPMixtureVarSet& vs) = 0;
    virtual void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    virtual void FlashDer(OCPMixtureVarSet& vs) = 0;

protected:
    EnthalpyCalculation eC;

};


/////////////////////////////////////////////////////
// OCPMixtureThermalOWMethod01
/////////////////////////////////////////////////////


/// Oil and Water are immiscible
class OCPMixtureThermalOWMethod01 : public OCPMixtureThermalOWMethod
{
public:
    OCPMixtureThermalOWMethod01(const ComponentParam& param, const USI& tarId, OCPMixtureVarSet& vs);
    OCP_DBL CalRhoO(const OCP_DBL& P, const OCP_DBL& T) override;
    OCP_DBL CalXiO(const OCP_DBL& P, const OCP_DBL& T) override;
    OCP_DBL CalRhoW(const OCP_DBL& P, const OCP_DBL& T) override;
    OCP_DBL CalXiW(const OCP_DBL& P, const OCP_DBL& T) override;
    void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    void Flash(OCPMixtureVarSet& vs) override;
    void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs)override;
    void FlashDer(OCPMixtureVarSet& vs) override;

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
           
    /// Coefficients in water and oil viscosity correlation formula
    vector<OCP_DBL> avisc;   
    /// Coefficients in water and oil viscosity correlation formulae
    vector<OCP_DBL> bvisc;       
    /// Coefficients Ak in gas viscosity correlation formulae
    vector<OCP_DBL> avg;       
    /// Coefficients Bk in gas viscosity correlation formulae
    vector<OCP_DBL> bvg;   
    /// If use viscosity to calculate viscosity
    OCP_BOOL        useViscTab;
    /// Viscosity-versus-temperature dependence
    OCPTable        visc; 
    /// Auxiliary variables used to calculate viscosity with interpolating tab
    vector<OCP_DBL> data;  
    ///< Auxiliary variables used to calculate viscosity with interpolating tab
    vector<OCP_DBL> cdata; 
};


/////////////////////////////////////////////////////
// OCPMixtureThermalOW 
/////////////////////////////////////////////////////

class OCPMixtureThermalOW : public OCPMixture
{
public:
    OCPMixtureThermalOW() { mixtureType = OCPMIXTURE_THERMALK_OW; }
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
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& T, const USI& tarPhase) {
        if (tarPhase == OIL)         return pmMethod->CalXiO(P, T);
        else if (tarPhase == WATER)  return pmMethod->CalXiW(P, T);
        else                         OCP_ABORT("WRONG TarPhase");
    }
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& T, const USI& tarPhase) {
        if (tarPhase == OIL)         return pmMethod->CalRhoO(P, T);
        else if (tarPhase == WATER)  return pmMethod->CalRhoW(P, T);
        else                         OCP_ABORT("WRONG TarPhase");
    }

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
        vs.S[1] = Sw;
    }

protected:
    OCPMixtureThermalOWMethod* pmMethod;
};


#endif /* end if __OCPMIXTURETHERMALOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/20/2023      Create file                          */
/*----------------------------------------------------------------------------*/