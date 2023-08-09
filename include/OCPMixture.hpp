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

#include "OCPConst.hpp"

#include <vector>

using namespace std;

const USI OCPMIXTURE_SP          = 0;
const USI OCPMIXTURE_BO_OG       = 121;
const USI OCPMIXTURE_BO_OW       = 122;
const USI OCPMIXTURE_BO_GW       = 123;
const USI OCPMIXTURE_BO_OGW      = 13;
const USI OCPMIXTURE_COMP        = 2;
const USI OCPMIXTURE_THERMALK_OW = 3;

/// mixture varset
class OCPMixtureVarSet
{
public:
    OCPMixtureVarSet() = default;
    void Init(const USI& numPhase, const USI& numCom, const OCP_BOOL& ifThermal) {
        // vars
        np = numPhase;
        nc = numCom;
        Ni.resize(nc);
        phaseExist.resize(np);
        S.resize(np);
        nj.resize(np);
        vj.resize(np);
        x.resize(np * nc);
        rho.resize(np);
        xi.resize(np);
        mu.resize(np);

        // derivatives
        vfi.resize(nc);
        vjP.resize(np);
        if (ifThermal) vjT.resize(np);
        vji.resize(np);
        for (auto& v : vji) v.resize(nc);
        rhoP.resize(np);
        if (ifThermal) rhoT.resize(np);
        rhox.resize(np * nc);
        xiP.resize(np);
        if (ifThermal) xiT.resize(np);
        xix.resize(np * nc);
        muP.resize(np);
        if (ifThermal) muT.resize(np);
        mux.resize(np * nc);
        if (ifThermal) Ufi.resize(nc);
        if (ifThermal) H.resize(np);
        if (ifThermal) HT.resize(np);
        if (ifThermal) Hx.resize(np * nc);

        if (ifThermal) dXsdXp.resize(np * (nc + 1) * (nc + 2));
        else           dXsdXp.resize(np * (nc + 1) * (nc + 1));
    }

public:
    /// Calculate total fluid volume and phase saturation with vj
    void CalVfS() {
        Vf = 0;
        for (USI j = 0; j < np; j++) {
            if (phaseExist[j]) {
                Vf += vj[j];
            }
        }
        for (USI j = 0; j < np; j++) {
            S[j] = 0;
            if (phaseExist[j]) {
                S[j] = vj[j] / Vf;
            }
        }
    }

public:
    const OCP_DBL* GetXj(const USI& j) const { return &x[j * nc]; }

public:
    /// num of phase, components
    USI                     np, nc;
    /// pressure, temperature
    OCP_DBL                 P, T;
    /// Buble point pressure 
    OCP_DBL                 Pb;
    /// total moles of components
    OCP_DBL                 Nt;
    /// total volume of components
    OCP_DBL                 Vf;
    /// mole number of components
    vector<OCP_DBL>         Ni;
    /// existence of phase
    vector<OCP_BOOL>        phaseExist;
    /// saturation of phase
    vector<OCP_DBL>         S;
    /// mole number of phases
    vector<OCP_DBL>         nj;
    /// volume of phases
    vector<OCP_DBL>         vj;
    /// molar fraction of component i in phase j
    vector<OCP_DBL>         x;
    /// mass density of phases
    vector<OCP_DBL>         rho;
    /// molar density of phases
    vector<OCP_DBL>         xi;
    /// viscosity of phases
    vector<OCP_DBL>         mu;

    // Derivatives (full derivatives)
    /// dVf / dP
    OCP_DBL                 vfP;
    /// dVf / dT
    OCP_DBL                 vfT; 
    /// dVf / dNi
    vector<OCP_DBL>         vfi;
    /// dVj / dP            
    vector<OCP_DBL>         vjP;
    /// dVj / dT            
    vector<OCP_DBL>         vjT;
    /// dVj / dNi
    vector<vector<OCP_DBL>> vji;

    // Derivatives (partial derivatives)
    /// drho / dP
    vector<OCP_DBL> rhoP;
    /// drho / dT
    vector<OCP_DBL> rhoT; 
    /// drho / dx 
    vector<OCP_DBL> rhox;
    /// dxi / dP
    vector<OCP_DBL> xiP;  
    /// dxi / dT
    vector<OCP_DBL> xiT;  
    /// dxi / dx
    vector<OCP_DBL> xix;
    /// dmu / dP
    vector<OCP_DBL> muP;  
    /// dmu / dT
    vector<OCP_DBL> muT;  
    /// dmu / dx
    vector<OCP_DBL> mux;  

    /// Internal energy of fluid
    OCP_DBL Uf;
    /// dUf / dP
    OCP_DBL UfP;
    /// dUf / dT
    OCP_DBL UfT;
    /// dUf / dNi
    vector<OCP_DBL> Ufi; 
    /// Enthalpy
    vector<OCP_DBL> H;  
    /// d Hj / d T
    vector<OCP_DBL> HT;  
    ///< d Hj / d xij
    vector<OCP_DBL> Hx;  

    /// d(Sj, xij) / d(P,Ni,(T))
    vector<OCP_DBL> dXsdXp;
};


/////////////////////////////////////////////////////
// OCPMixture
/////////////////////////////////////////////////////

class OCPMixture
{
public:
    OCPMixture() = default;
    USI MixtureType() const { return mixtureType; }
    const OCPMixtureVarSet& GetVarSet() const { return vs; }

protected:
    /// mixture type
    USI              mixtureType;
    /// mixture variables set
    OCPMixtureVarSet vs;
};


#endif /* end if __OCPMIXTURE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/12/2023      Create file                          */
/*----------------------------------------------------------------------------*/