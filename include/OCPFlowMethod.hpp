/*! \file    OCPFlowMethod.hpp
 *  \brief   OCPFlowMethod class declaration
 *  \author  Shizhe Li
 *  \date    Oct/04/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPFLOWMETHOD_HEADER__
#define __OCPFLOWMETHOD_HEADER__

#include "OCPFlowVarSet.hpp"
#include "OCPFuncSAT.hpp"
#include <vector>

using namespace std;


/// Used to calculate the properties of flow, for example, phase relative permeability,
/// capillary pressure and some derived properties
class OCPFlowMethod
{
public:
    OCPFlowMethod() = default;
    /// calculate relative permeability and capillary pressure
    virtual void CalKrPc(OCPFlowVarSet& vs) = 0;
    /// calculate relative permeability and capillary pressure and derivatives
    virtual void CalKrPcDer(OCPFlowVarSet& vs) = 0;
    /// get saturation of connate water
    virtual OCP_DBL GetSwco() const = 0;
    /// get maximum capillary pressure between water and oil (Po - Pw)
    virtual OCP_DBL GetMaxPcow() const = 0;
    /// get minimum capillary pressure between water and oil (Po - Pw)
    virtual OCP_DBL GetMinPcow() const = 0;
    /// calculate Pcow by Sw
    virtual OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const = 0;
    /// calculate Sw by Pcow 
    virtual OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const = 0;
    /// calculate Pcgo by Sg
    virtual OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const = 0;
    /// calculate Sg by Pcgo
    virtual OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const = 0;
    /// calculate Sw by Pcgw
    virtual OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const = 0;
    /// calculate Krg and ders by Sg
    virtual OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const = 0;
};


/////////////////////////////////////////////////////
// OCPFlowMethod_OGW01
/////////////////////////////////////////////////////


/// Use SGOF, SWOF to calculate phase flow properties in oil-gas-water situation
class OCPFlowMethod_OGW01 : public OCPFlowMethod
{
public:
    OCPFlowMethod_OGW01(const vector<vector<OCP_DBL>>& SGOFin,
        const vector<vector<OCP_DBL>>& SWOFin,
        const USI& i, OCPFlowVarSet& vs);
    void CalKrPc(OCPFlowVarSet& vs) override;
    void CalKrPcDer(OCPFlowVarSet& vs) override;
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
    /// generate table of water saturation vs. Pcgw
    void Generate_SWPCWG();

protected:
    /// table SGOF
    OCP_SGOF                SGOF;
    /// table SWOF
    OCP_SWOF                SWOF;
    /// table SWPCGW(auxiliary table: water saturation vs. Pcgw)
    OCPTable                SWPCGW;

protected:
    /// calculations for oil relative permeability in oil-gas-water flow situation
    OCP3POilPerCalculation  opC;
};


/////////////////////////////////////////////////////
// OCPFlowMethod_OGW02
/////////////////////////////////////////////////////


/// Use SOF3, SGFN, SWFN to calculate phase flow properties in oil-gas-water situation
class OCPFlowMethod_OGW02 : public OCPFlowMethod
{
public:
    OCPFlowMethod_OGW02(const vector<vector<OCP_DBL>>& SOF3in,
        const vector<vector<OCP_DBL>>& SGFNin,
        const vector<vector<OCP_DBL>>& SWFNin,
        const USI& i, OCPFlowVarSet& vs);
    void CalKrPc(OCPFlowVarSet& vs) override;
    void CalKrPcDer(OCPFlowVarSet& vs) override;

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
    /// generate table of water saturation vs. Pcgw
    void Generate_SWPCWG();

protected:
    /// table SOF3
    OCP_SOF3            SOF3;
    /// table SGFN
    OCP_SGFN            SGFN;
    /// table SWPCGW(auxiliary table: water saturation vs. Pcgw)
    OCP_SWFN            SWFN;
    /// auxiliary table: saturation of water vs. Pcgw
    OCPTable            SWPCGW; 


protected:
    /// calculations for oil relative permeability in oil-gas-water flow situation
    OCP3POilPerCalculation  opC;
};


/////////////////////////////////////////////////////
// OCPFlowMethod_OW01
/////////////////////////////////////////////////////


/// Use SWOF to calculate phase flow properties in oil-water situation
class OCPFlowMethod_OW01 : public OCPFlowMethod
{
public:
    OCPFlowMethod_OW01(const vector<vector<OCP_DBL>>& SWOFin, OCPFlowVarSet& vs);
    void CalKrPc(OCPFlowVarSet& vs) override;
    void CalKrPcDer(OCPFlowVarSet& vs) override;
    OCP_DBL GetSwco() const override { return SWOF.GetSwco(); }
    OCP_DBL GetMaxPcow() const override { return SWOF.GetMaxPc(); }
    OCP_DBL GetMinPcow() const override { return SWOF.GetMinPc(); }
    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { return SWOF.CalPcow(Sw); }
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const override { return SWOF.CalSw(Pcow); }
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const  override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const override { OCP_ABORT("Inavailable!"); }

protected:
    /// table SWOF
    OCP_SWOF            SWOF;
};


/////////////////////////////////////////////////////
// OCPFlowMethod_OG01
/////////////////////////////////////////////////////


/// Use SGOF to calculate phase flow properties in oil-gas situation
class OCPFlowMethod_OG01 : public OCPFlowMethod
{
public:
    OCPFlowMethod_OG01(const vector<vector<OCP_DBL>>& SGOFin, OCPFlowVarSet& vs);
    void CalKrPc(OCPFlowVarSet& vs) override;
    void CalKrPcDer(OCPFlowVarSet& vs) override;
    OCP_DBL GetSwco() const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL GetMaxPcow() const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL GetMinPcow() const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const override { return SGOF.CalPcgo(Sg); }
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const override { return SGOF.CalSg(Pcgo); }
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const override { return SGOF.CalKrg(Sg, dKrgdSg); }

protected:
    /// table SGOF
    OCP_SGOF            SGOF;
};


/////////////////////////////////////////////////////
// OCPFlowMethod_GW01
/////////////////////////////////////////////////////


/// Use Brooks-Corey type model to calculate phase flow properties in gas-water situation
class OCPFlowMethod_GW01 : public OCPFlowMethod
{
public:
    OCPFlowMethod_GW01(const BrooksCoreyParam& bcp, OCPFlowVarSet& vs) {
        vs.Init(OCPFlowType::GW, 2, 2);
        bc.Setup(bcp);
    }
    void CalKrPc(OCPFlowVarSet& vs) override;
    void CalKrPcDer(OCPFlowVarSet& vs) override;
    OCP_DBL GetSwco() const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL GetMaxPcow() const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL GetMinPcow() const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const override { OCP_ABORT("Inavailable!"); }

protected:
    /// Brooks-Corey type model
    BrooksCorey   bc;
};


/////////////////////////////////////////////////////
// OCPFlowMethod_GW02
/////////////////////////////////////////////////////


/// Use SWGF to calculate phase flow properties in gas-water situation
class OCPFlowMethod_GW02 : public OCPFlowMethod
{
public:
    OCPFlowMethod_GW02(const vector<vector<OCP_DBL>>& SWGFin, OCPFlowVarSet& vs) {
        vs.Init(OCPFlowType::GW, 2, 2);
        SWGF.Setup(SWGFin);
    }
    void CalKrPc(OCPFlowVarSet& vs) override;
    void CalKrPcDer(OCPFlowVarSet& vs) override;
    OCP_DBL GetSwco() const override { return SWGF.GetSwco(); }
    OCP_DBL GetMaxPcow() const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL GetMinPcow() const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const override { OCP_ABORT("Inavailable!"); }
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const override { return SWGF.CalSw(Pcgw); }
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const override { OCP_ABORT("Inavailable!"); }

protected:
    /// table SWGF
    OCP_SWGF            SWGF;
};


#endif /* end if __OCPFLOWMETHOD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/04/2023      Create file                          */
/*----------------------------------------------------------------------------*/