/*! \file    OCPFlowOGW.hpp
 *  \brief   OCPFlowOGW class declaration
 *  \author  Shizhe Li
 *  \date    Jul/08/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPFLOWOGW_HEADER__
#define __OCPFLOWOGW_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPSATFunc.hpp"
#include "OCPFlow.hpp"

#include <vector>

using namespace std;


/// Calculate oil permeability
class OCP3POilPerMethod
{
public:
    virtual void CalOilPer(OCPFlowVarSet* vs) = 0;
    virtual void CalOilPerDer(OCPFlowVarSet* vs) = 0;
};


/////////////////////////////////////////////////////
// OCP3POilPerMethod01
/////////////////////////////////////////////////////

/// Stone 2
class OCP3POilPerMethod01 : public OCP3POilPerMethod
{
public:
    OCP3POilPerMethod01() = default;
    void CalOilPer(OCPFlowVarSet* vs) override;
    void CalOilPerDer(OCPFlowVarSet* vs) override;
};


/////////////////////////////////////////////////////
// OCP3POilPerMethod02
/////////////////////////////////////////////////////


/// Eclipse default
class OCP3POILPerMethod02 : public OCP3POilPerMethod
{
public:
    OCP3POILPerMethod02() = default;
    void CalOilPer(OCPFlowVarSet* vs) override;
    void CalOilPerDer(OCPFlowVarSet* vs) override;
};

/// Calculate oil, gas, water relative permeability and capillary pressure
class OCPOGWFMethod
{
public:
    OCPOGWFMethod() = default;
    virtual void CalKrPc() = 0;
    virtual void CalKrPcDer() = 0;

    virtual OCP_DBL GetSwco() const = 0;
    virtual OCP_DBL GetMaxPcow() const = 0;
    virtual OCP_DBL GetMinPcow() const = 0;

    virtual OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const = 0;
    virtual OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const = 0;
    virtual OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const = 0;
    virtual OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const = 0;
    virtual OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const = 0;
    virtual OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const = 0;
    

protected:
    OCPFlowVarSet*  vs;
};


/////////////////////////////////////////////////////
// OCPOGWFMethod01
/////////////////////////////////////////////////////


/// Use SGOF, SWOF
class OCPOGWFMethod01 : public OCPOGWFMethod
{
public:
    OCPOGWFMethod01(const vector<vector<OCP_DBL>>& SGOFin,
        const vector<vector<OCP_DBL>>& SWOFin,
        const USI& i, OCPFlowVarSet* vsin);
    void CalKrPc() override;
    void CalKrPcDer() override;

    OCP_DBL GetSwco() const override { return SWOF.GetSwco(); }
    OCP_DBL GetMaxPcow() const override { return SWOF.GetMaxPc(); }
    OCP_DBL GetMinPcow() const override { return SWOF.GetMinPc(); }

    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { return SWOF.CalPcow(Sw); }
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const override { return SWOF.CalSw(Pcow); }
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const override { return SGOF.CalPcgo(Sg); }
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const override { return SGOF.CalSg(Pcgo); }
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const override { return SWPCGW.Eval_Inv(1, Pcgw, 0); }
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const override { return SGOF.CalKrg(Sg, dKrgdSg); }

protected:
    void Generate_SWPCWG();

protected:
    OCP_SGOF            SGOF;
    OCP_SWOF            SWOF;
    OCP3POilPerMethod*  opMethod;
    OCPTable            SWPCGW; ///< auxiliary table: saturation of water vs. Pcgw
};


/////////////////////////////////////////////////////
// OCPOGWFMethod02
/////////////////////////////////////////////////////


/// Use SOF3, SGFN, SWFN
class OCPOGWFMethod02 : public OCPOGWFMethod
{
public:
    OCPOGWFMethod02(const vector<vector<OCP_DBL>>& SOF3in,
        const vector<vector<OCP_DBL>>& SGFNin,
        const vector<vector<OCP_DBL>>& SWFNin,
        const USI& i, OCPFlowVarSet* vsin);
    void CalKrPc() override;
    void CalKrPcDer() override;

    OCP_DBL GetSwco() const override { return SWFN.GetSwco(); }
    OCP_DBL GetMaxPcow() const override { return SWFN.GetMaxPc(); }
    OCP_DBL GetMinPcow() const override { return SWFN.GetMinPc(); }

    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { return SWFN.CalPcow(Sw); }
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const override { return SWFN.CalSw(Pcow); }
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const override { return SGFN.CalPcgo(Sg); }
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const override { return SGFN.CalSg(Pcgo); }
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const override { return SWPCGW.Eval_Inv(1, Pcgw, 0); }
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const override { return SGFN.CalKrg(Sg, dKrgdSg); }

protected:
    void Generate_SWPCWG();

protected:
    OCP_SOF3            SOF3;
    OCP_SGFN            SGFN;
    OCP_SWFN            SWFN;
    OCP3POilPerMethod*  opMethod;
    OCPTable            SWPCGW; ///< auxiliary table: saturation of water vs. Pcgw
};


/////////////////////////////////////////////////////
// OCPFlowOGW
/////////////////////////////////////////////////////

class OCPFlowOGW : public OCPFlow
{
public:
    OCPFlowOGW() = default;
    void Setup(const ParamReservoir& rs_param, const USI& i);
    OCPFlowVarSet& GetVarSet() { return vs; }
    void CalKrPc(const OCP_DBL& So, const OCP_DBL& Sg, const OCP_DBL& Sw) { 
        SetSaturation(So, Sg, Sw);
        pfMethod->CalKrPc(); 
    }
    void CalKrPcDer(const OCP_DBL& So, const OCP_DBL& Sg, const OCP_DBL& Sw) { 
        SetSaturation(So, Sg, Sw);
        pfMethod->CalKrPcDer(); 
    }

    OCP_DBL GetSwco() const { return pfMethod->GetSwco(); }
    OCP_DBL GetMaxPcow() const { return pfMethod->GetMaxPcow(); }
    OCP_DBL GetMinPcow() const { return pfMethod->GetMinPcow(); }

    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const { return pfMethod->CalPcowBySw(Sw); }
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const { return pfMethod->CalSwByPcow(Pcow); }
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const { return pfMethod->CalPcgoBySg(Sg); }
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const { return pfMethod->CalSgByPcgo(Pcgo); }
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const { return pfMethod->CalSwByPcgw(Pcgw); }
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) { return pfMethod->CalKrg(Sg, dKrgdSg); }

protected:
    void SetSaturation(const OCP_DBL& So, const OCP_DBL& Sg, const OCP_DBL& Sw) {
        vs.So = So;
        vs.Sg = Sg;
        vs.Sw = Sw;
    }

protected:
    const USI       flowType = OCPFLOW_OGW;
    OCPOGWFMethod*  pfMethod;
};


#endif /* end if __OCP3PHASEFLOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/08/2023      Create file                          */
/*----------------------------------------------------------------------------*/