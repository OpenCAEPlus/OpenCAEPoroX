/*! \file    OCPFlowOG.hpp
 *  \brief   OCPFlowOG class declaration
 *  \author  Shizhe Li
 *  \date    Jul/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPFLOWOG_HEADER__
#define __OCPFLOWOG_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPFuncSAT.hpp"
#include "OCPFlow.hpp"

#include <vector>

using namespace std;


/// Calculate oil, gas, water relative permeability and capillary pressure
class OCPOGFMethod
{
public:
    OCPOGFMethod() = default;
    virtual void CalKrPc(OCPFlowVarSet& vs) = 0;
    virtual void CalKrPcDer(OCPFlowVarSet& vs) = 0;
    virtual OCP_DBL CalPcgoBySg(const OCP_DBL& sg) const = 0;
    virtual OCP_DBL CalSgByPcgo(const OCP_DBL& pcgo) const = 0;
    virtual OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const = 0;
};


/////////////////////////////////////////////////////
// OCPOGFMethod01
/////////////////////////////////////////////////////


/// Use SGOF
class OCPOGFMethod01 : public OCPOGFMethod
{
public:
    OCPOGFMethod01(const vector<vector<OCP_DBL>>& SGOFin);
    void CalKrPc(OCPFlowVarSet& vs) override;
    void CalKrPcDer(OCPFlowVarSet& vs) override;

    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const override { return SGOF.CalPcgo(Sg); }
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const override { return SGOF.CalSg(Pcgo); }
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const override { return SGOF.CalKrg(Sg, dKrgdSg); }

protected:
    OCP_SGOF            SGOF;
};


/////////////////////////////////////////////////////
// OCPFlowOG
/////////////////////////////////////////////////////

class OCPFlowOG : public OCPFlow
{
public:
    OCPFlowOG() { flowType = OCPFlowType::OG; }
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
        vs.So = So;
        vs.Sg = Sg;
    }

protected:
    OCPOGFMethod*  pfMethod;
};


#endif /* end if __OCPFLOWOG_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/