/*! \file    OCP3PhaseFlow.hpp
 *  \brief   OCP3PhaseFlow class declaration
 *  \author  Shizhe Li
 *  \date    Jul/08/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCP3PHASEFLOW_HEADER__
#define __OCP3PHASEFLOW_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPSATFunc.hpp"

#include <vector>

using namespace std;

/// oil-gas-water flow vars suite
class OCP3PFVarSet
{
    /// oil phase is reference phase
public:
    OCP3PFVarSet() { Init0(); }
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


/// Calculate oil permeability
class OCP3POilPerMethod
{
public:
    virtual void CalOilPer(OCP3PFVarSet* vs) = 0;
    virtual void CalOilPerDer(OCP3PFVarSet* vs) = 0;
};


/////////////////////////////////////////////////////
// OCP3POilPerMethod01
/////////////////////////////////////////////////////

/// Stone 2
class OCP3POilPerMethod01 : public OCP3POilPerMethod
{
public:
    OCP3POilPerMethod01() = default;
    void CalOilPer(OCP3PFVarSet* vs) override;
    void CalOilPerDer(OCP3PFVarSet* vs) override;
};


/////////////////////////////////////////////////////
// OCP3POilPerMethod02
/////////////////////////////////////////////////////


/// Eclipse default
class OCP3POILPerMethod02 : public OCP3POilPerMethod
{
public:
    OCP3POILPerMethod02() = default;
    void CalOilPer(OCP3PFVarSet* vs) override;
    void CalOilPerDer(OCP3PFVarSet* vs) override;
};

/// Calculate oil, gas, water relative permeability and capillary pressure
class OCP3PFMethod
{
public:
    OCP3PFMethod() = default;
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
    OCP3PFVarSet*  vs;
};


/////////////////////////////////////////////////////
// OCP3PFMethod01
/////////////////////////////////////////////////////


/// Use SGOF, SWOF
class OCP3PFMethod01 : public OCP3PFMethod
{
public:
    OCP3PFMethod01(const vector<vector<OCP_DBL>>& SGOFin,
        const vector<vector<OCP_DBL>>& SWOFin,
        const USI& i, OCP3PFVarSet* vsin);
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
// OCP3PFMethod02
/////////////////////////////////////////////////////


/// Use SOF3, SGFN, SWFN
class OCP3PFMethod02 : public OCP3PFMethod
{
public:
    OCP3PFMethod02(const vector<vector<OCP_DBL>>& SOF3in,
        const vector<vector<OCP_DBL>>& SGFNin,
        const vector<vector<OCP_DBL>>& SWFNin,
        const USI& i, OCP3PFVarSet* vsin);
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
// OCP3PhaseFlow
/////////////////////////////////////////////////////

class OCP3PhaseFlow
{
public:
    OCP3PhaseFlow() = default;
    void Setup(const ParamReservoir& rs_param, const USI& i);
    OCP3PFVarSet& GetVarSet() { return vs; }
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
    OCP3PFVarSet   vs;
    OCP3PFMethod*  pfMethod;
};


#endif /* end if __OCP3PHASEFLOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/08/2023      Create file                          */
/*----------------------------------------------------------------------------*/