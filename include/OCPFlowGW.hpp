/*! \file    OCPFlowGW.hpp
 *  \brief   OCPFlowGW class declaration
 *  \author  Shizhe Li
 *  \date    Sep/30/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPFLOWGW_HEADER__
#define __OCPFLOWGW_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPFuncSAT.hpp"
#include "OCPFlow.hpp"

#include <vector>

using namespace std;


/// Calculate gas, water relative permeability and capillary pressure
class OCPGWFMethod
{
public:
    OCPGWFMethod() = default;
    virtual void CalKrPc(OCPFlowVarSet& vs) = 0;
    virtual void CalKrPcDer(OCPFlowVarSet& vs) = 0;
};


/////////////////////////////////////////////////////
// OCPGWFMethod01
/////////////////////////////////////////////////////


/// Use Brooks-Corey type model
class OCPGWFMethod01 : public OCPGWFMethod
{
public:
    OCPGWFMethod01(const BrooksCoreyParam& bcp) { bc.Setup(bcp); }
    void CalKrPc(OCPFlowVarSet& vs) override;
    void CalKrPcDer(OCPFlowVarSet& vs) override;

protected:
    BrooksCorey   bc;
};


/////////////////////////////////////////////////////
// OCPFlowGW
/////////////////////////////////////////////////////

class OCPFlowGW : public OCPFlow
{
public:
    OCPFlowGW() { flowType = OCPFlowType::GW; }
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
        vs.Sg = Sg;
        vs.Sw = Sw;
    }

protected:
    OCPGWFMethod* pfMethod;
};


#endif /* end if __OCPFLOWGW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/30/2023      Create file                          */
/*----------------------------------------------------------------------------*/