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

class OCPFlow
{
public:
    OCPFlow() = default;
    auto FlowType() const { return vs.flowType; }
    OCPFlowVarSet& GetVarSet() { return vs; }
    
    virtual OCP_DBL GetSwco() const = 0;
    virtual OCP_DBL GetMaxPcow() const = 0;
    virtual OCP_DBL GetMinPcow() const = 0;
    virtual OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const = 0;
    virtual OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const = 0;
protected:
    
    OCPFlowVarSet vs;
};


/////////////////////////////////////////////////////
// OCPFlowOGW
/////////////////////////////////////////////////////

/// There exists oil, gas and water, the reference phase is oil
class OCPFlowOGW : public OCPFlow
{
public:
    void Setup(const ParamReservoir& rs_param, const USI& i);
    void CalKrPc(const OCP_DBL& So, const OCP_DBL& Sg, const OCP_DBL& Sw) {
        SetSaturation(So, Sg, Sw);
        pfMethod->CalKrPc(vs);
    }
    void CalKrPcDer(const OCP_DBL& So, const OCP_DBL& Sg, const OCP_DBL& Sw) {
        SetSaturation(So, Sg, Sw);
        pfMethod->CalKrPcDer(vs);
    }

    OCP_DBL GetSwco() const override { return pfMethod->GetSwco(); }
    OCP_DBL GetMaxPcow() const override { return pfMethod->GetMaxPcow(); }
    OCP_DBL GetMinPcow() const override { return pfMethod->GetMinPcow(); }

    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { return pfMethod->CalPcowBySw(Sw); }
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const { return pfMethod->CalSwByPcow(Pcow); }
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const { return pfMethod->CalPcgoBySg(Sg); }
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const { return pfMethod->CalSgByPcgo(Pcgo); }
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const { return pfMethod->CalSwByPcgw(Pcgw); }
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const override { return pfMethod->CalKrg(Sg, dKrgdSg); }

protected:
    void SetSaturation(const OCP_DBL& So, const OCP_DBL& Sg, const OCP_DBL& Sw) {
        vs.S[vs.o] = So;
        vs.S[vs.g] = Sg;
        vs.S[vs.w] = Sw;
    }

protected:
    OCPFlowMethod* pfMethod;
};


/////////////////////////////////////////////////////
// OCPFlowOW
/////////////////////////////////////////////////////

/// There exists oil and water, the reference phase is oil
class OCPFlowOW : public OCPFlow
{
public:
    void Setup(const ParamReservoir& rs_param, const USI& i);
    void CalKrPc(const OCP_DBL& So, const OCP_DBL& Sw) {
        SetSaturation(So, Sw);
        pfMethod->CalKrPc(vs);
    }
    void CalKrPcDer(const OCP_DBL& So, const OCP_DBL& Sw) {
        SetSaturation(So, Sw);
        pfMethod->CalKrPcDer(vs);
    }

    OCP_DBL GetSwco() const override { return pfMethod->GetSwco(); }
    OCP_DBL GetMaxPcow() const override { return pfMethod->GetMaxPcow(); }
    OCP_DBL GetMinPcow() const override { return pfMethod->GetMinPcow(); }

    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { return pfMethod->CalPcowBySw(Sw); }
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const { return pfMethod->CalSwByPcow(Pcow); }

    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const override { OCP_ABORT("Wrong Usage!"); }

protected:
    void SetSaturation(const OCP_DBL& So, const OCP_DBL& Sw) {
        vs.S[vs.o] = So;
        vs.S[vs.w] = Sw;
    }

protected:
    OCPFlowMethod* pfMethod;
};


/////////////////////////////////////////////////////
// OCPFlowOG
/////////////////////////////////////////////////////

/// There exists oil and gas, the reference phase is oil
class OCPFlowOG : public OCPFlow
{
public:
    void Setup(const ParamReservoir& rs_param, const USI& i);
    void CalKrPc(const OCP_DBL& So, const OCP_DBL& Sg) {
        SetSaturation(So, Sg);
        pfMethod->CalKrPc(vs);
    }
    void CalKrPcDer(const OCP_DBL& So, const OCP_DBL& Sg) {
        SetSaturation(So, Sg);
        pfMethod->CalKrPcDer(vs);
    }

    OCP_DBL GetSwco() const override { OCP_ABORT("Wrong Usage!"); }
    OCP_DBL GetMaxPcow() const override { OCP_ABORT("Wrong Usage!"); }
    OCP_DBL GetMinPcow() const override { OCP_ABORT("Wrong Usage!"); }

    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const { return pfMethod->CalSgByPcgo(Pcgo); }
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const { return pfMethod->CalPcgoBySg(Sg); }
    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { OCP_ABORT("Wrong Usage!"); }
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const override { return pfMethod->CalKrg(Sg, dKrgdSg); }

protected:
    void SetSaturation(const OCP_DBL& So, const OCP_DBL& Sg) {
        vs.S[vs.o] = So;
        vs.S[vs.g] = Sg;
    }

protected:
    OCPFlowMethod* pfMethod;
};


/////////////////////////////////////////////////////
// OCPFlowGW
/////////////////////////////////////////////////////

class OCPFlowGW : public OCPFlow
{
public:
    void Setup(const ParamReservoir& rs_param, const USI& i);
    void CalKrPc(const OCP_DBL& Sg, const OCP_DBL& Sw) {
        SetSaturation(Sg, Sw);
        pfMethod->CalKrPc(vs);
    }
    void CalKrPcDer(const OCP_DBL& Sg, const OCP_DBL& Sw) {
        SetSaturation(Sg, Sw);
        pfMethod->CalKrPcDer(vs);
    }

    OCP_DBL GetSwco() const override { OCP_ABORT("Wrong Usage!"); }
    OCP_DBL GetMaxPcow() const override { OCP_ABORT("Wrong Usage!"); }
    OCP_DBL GetMinPcow() const override { OCP_ABORT("Wrong Usage!"); }
    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { OCP_ABORT("Wrong Usage!"); }
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const override { OCP_ABORT("Wrong Usage!"); }

protected:
    void SetSaturation(const OCP_DBL& Sg, const OCP_DBL& Sw) {
        vs.S[vs.g] = Sg;
        vs.S[vs.w] = Sw;
    }

protected:
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