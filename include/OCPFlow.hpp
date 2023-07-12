/*! \file    OCPFlow.hpp
 *  \brief   OCPFlow class declaration
 *  \author  Shizhe Li
 *  \date    Jul/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPFLOW_HEADER__
#define __OCPFLOW_HEADER__

#include "OCPConst.hpp"

using namespace std;

const USI OCPFLOW_OG  = 21;
const USI OCPFLOW_OW  = 22;
const USI OCPFLOW_GW  = 23;
const USI OCPFLOW_OGW = 3;

/// flow vars suite
class OCPFlowVarSet
{
    /// oil phase is reference phase
public:
    OCPFlowVarSet() { Init0(); }
    void Init0() {
		Swco     = 0;
        krocw    = 0;
		So       = 0; Sg       = 0; Sw      = 0;		          
		kro      = 0; krg      = 0; krw     = 0;
		krog     = 0; krow     = 0;
		dKrodSo  = 0; dKrodSg  = 0; dKrodSw = 0;
		dKrgdSo  = 0; dKrgdSg  = 0; dKrgdSw = 0;
		dKrwdSo  = 0; dKrwdSg  = 0; dKrwdSw = 0;
		dKrogdSo = 0; dKrogdSg = 0;
		dKrowdSo = 0; dKrowdSw = 0;
		Pcgo     = 0; Pcwo     = 0;
		dPcgodSo = 0; dPcgodSg = 0; dPcgodSw = 0;
		dPcwodSo = 0; dPcwodSg = 0; dPcwodSw = 0;
    }

public:
    /// saturaion of connate water
    OCP_DBL Swco;
    /// oil relative permeability in the presence of connate water only
    OCP_DBL krocw;
    /// oil, gas, water saturations
    OCP_DBL So, Sg, Sw;
    /// oil, gas, water relatve permeability
    OCP_DBL kro, krg, krw;
    /// the corresponding oil relative permeability when oil, gas and connate water are present
    OCP_DBL krog;
    /// the corresponding oil relative permeability when only oil and water are present
    OCP_DBL krow;
    /// the corresponding derivatives of permeability
    OCP_DBL dKrodSo, dKrodSg, dKrodSw;
    OCP_DBL dKrgdSo, dKrgdSg, dKrgdSw;
    OCP_DBL dKrwdSo, dKrwdSg, dKrwdSw;
    OCP_DBL dKrogdSo, dKrogdSg;
    OCP_DBL dKrowdSo, dKrowdSw;

    /// Capillary pressure   
    OCP_DBL Pcgo, Pcwo;               ///< Pg - Po, Pw - Po
    /// the corresponding derivatives of capillary pressure
    OCP_DBL dPcgodSo, dPcgodSg, dPcgodSw;
    OCP_DBL dPcwodSo, dPcwodSg, dPcwodSw;
};


/////////////////////////////////////////////////////
// OCPFlow
/////////////////////////////////////////////////////

class OCPFlow
{
public:
    OCPFlow() = default;
    USI FlowType() const { return flowType; }
    OCPFlowVarSet& GetVarSet() { return vs; }
    
    virtual OCP_DBL GetSwco() const = 0;
    virtual OCP_DBL GetMaxPcow() const = 0;
    virtual OCP_DBL GetMinPcow() const = 0;
    virtual OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const = 0;
    virtual OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const = 0;
protected:
    USI           flowType;
    OCPFlowVarSet vs;
};


#endif /* end if __OCPFLOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/