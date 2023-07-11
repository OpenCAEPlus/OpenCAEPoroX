/*! \file    OCPFlowOW.hpp
 *  \brief   OCPFlowOW class declaration
 *  \author  Shizhe Li
 *  \date    Jul/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPFLOWOW_HEADER__
#define __OCPFLOWOW_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPFuncSAT.hpp"
#include "OCPFlow.hpp"

#include <vector>

using namespace std;


/// Calculate oil, gas, water relative permeability and capillary pressure
class OCPOWFMethod
{
public:
    OCPOWFMethod() = default;
    virtual void CalKrPc() = 0;
    virtual void CalKrPcDer() = 0;

    virtual OCP_DBL GetSwco() const = 0;
    virtual OCP_DBL GetMaxPcow() const = 0;
    virtual OCP_DBL GetMinPcow() const = 0;

    virtual OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const = 0;
    virtual OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const = 0;

protected:
    OCPFlowVarSet* vs;
};


/////////////////////////////////////////////////////
// OCPOWFMethod01
/////////////////////////////////////////////////////


/// Use SWOF
class OCPOWFMethod01 : public OCPOWFMethod
{
public:
    OCPOWFMethod01(const vector<vector<OCP_DBL>>& SWOFin, OCPFlowVarSet* vsin);
    void CalKrPc() override;
    void CalKrPcDer() override;

    OCP_DBL GetSwco() const override { return SWOF.GetSwco(); }
    OCP_DBL GetMaxPcow() const override { return SWOF.GetMaxPc(); }
    OCP_DBL GetMinPcow() const override { return SWOF.GetMinPc(); }

    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { return SWOF.CalPcow(Sw); }
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const override { return SWOF.CalSw(Pcow); }

protected:
    OCP_SWOF            SWOF;
};


/////////////////////////////////////////////////////
// OCPFlowOW
/////////////////////////////////////////////////////

class OCPFlowOW : public OCPFlow
{
public:
    OCPFlowOW() { flowType = OCPFLOW_OW; }
    void Setup(const ParamReservoir& rs_param, const USI& i);
    void CalKrPc(const OCP_DBL& So, const OCP_DBL& Sw) {
        SetSaturation(So, Sw);
        pfMethod->CalKrPc();
    }
    void CalKrPcDer(const OCP_DBL& So, const OCP_DBL& Sw) {
        SetSaturation(So, Sw);
        pfMethod->CalKrPcDer();
    }

    OCP_DBL GetSwco() const override { return pfMethod->GetSwco(); }
    OCP_DBL GetMaxPcow() const override { return pfMethod->GetMaxPcow(); }
    OCP_DBL GetMinPcow() const override { return pfMethod->GetMinPcow(); }

    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { return pfMethod->CalPcowBySw(Sw); }
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const { return pfMethod->CalSwByPcow(Pcow); }

    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const override { OCP_ABORT("Wrong Usage!"); }

protected:
    void SetSaturation(const OCP_DBL& So, const OCP_DBL& Sw) {
        vs.So = So;
        vs.Sw = Sw;
    }

protected:
    OCPOWFMethod*  pfMethod;
};


#endif /* end if __OCPFLOWOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/