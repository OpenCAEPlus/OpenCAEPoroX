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


enum class OCPFlowType 
{
    /// single phase
    SP,
    /// oil and gas
    OG,
    /// oil and water
    OW,
    /// gas and water
    GW,
    /// oil, gas and water
    OGW
};

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
        Pco      = 0; Pcg      = 0; Pcw      = 0;
		dPcodSo  = 0; dPcodSg  = 0; dPcodSw  = 0;
		dPcgdSo  = 0; dPcgdSg  = 0; dPcgdSw  = 0;
		dPcwdSo  = 0; dPcwdSg  = 0; dPcwdSw  = 0;
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

    /// Capillary pressure, relative to reference phase : Po - Pr, Pg - Pr, Pw - Pr
    OCP_DBL Pco, Pcg, Pcw;
    /// the corresponding derivatives of capillary pressure
    OCP_DBL dPcodSo, dPcodSg, dPcodSw;
    OCP_DBL dPcgdSo, dPcgdSg, dPcgdSw;
    OCP_DBL dPcwdSo, dPcwdSg, dPcwdSw;
};


/////////////////////////////////////////////////////
// OCPFlow
/////////////////////////////////////////////////////

class OCPFlow
{
public:
    OCPFlow() = default;
    auto FlowType() const { return flowType; }
    OCPFlowVarSet& GetVarSet() { return vs; }
    
    virtual OCP_DBL GetSwco() const = 0;
    virtual OCP_DBL GetMaxPcow() const = 0;
    virtual OCP_DBL GetMinPcow() const = 0;
    virtual OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const = 0;
    virtual OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const = 0;
protected:
    OCPFlowType   flowType;
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