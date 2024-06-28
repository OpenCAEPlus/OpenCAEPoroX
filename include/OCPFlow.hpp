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

#include "ParamReservoir.hpp"
#include "OCPFlowMethod.hpp"


using namespace std;


/////////////////////////////////////////////////////
// OCPFlow
/////////////////////////////////////////////////////


/// There exists oil, gas and water, the reference phase is oil
class OCPFlow
{
public:
    OCPFlow(const ParamReservoir& rs_param, const USI& i);
    auto& GetVarSet() { return vs; }
    void CalKrPc(const OCP_DBL* S) {
        SetSaturation(S);
        pfMethod->CalKrPc(vs);
    }
    auto FlowType() const { return vs.flowType; }
    void CalKrPcDer(const OCP_DBL* S) {
        SetSaturation(S);
        pfMethod->CalKrPcDer(vs);
    }
    OCP_DBL GetSwco() const { return pfMethod->GetSwco(); }
    OCP_DBL GetMaxPcow() const { return pfMethod->GetMaxPcow(); }
    OCP_DBL GetMinPcow() const { return pfMethod->GetMinPcow(); }
    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const { return pfMethod->CalPcowBySw(Sw); }
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const { return pfMethod->CalSwByPcow(Pcow); }
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const { return pfMethod->CalPcgoBySg(Sg); }
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const { return pfMethod->CalSgByPcgo(Pcgo); }
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const { return pfMethod->CalSwByPcgw(Pcgw); }
    OCP_DBL CalPcgwBySw(const OCP_DBL& Sw) const { return pfMethod->CalPcgwBySw(Sw); }
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const  { return pfMethod->CalKrg(Sg, dKrgdSg); }

protected:
    void SetSaturation(const OCP_DBL* S) {
        copy(S, S + vs.np, vs.S.data());
    }

protected:
    OCPFlowVarSet vs;
    OCPFlowMethod* pfMethod;
};


#endif /* end if __OCPFLOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/