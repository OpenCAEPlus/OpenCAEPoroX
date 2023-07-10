/*! \file    OCPOWFlow.hpp
 *  \brief   OCPOWFlow class declaration
 *  \author  Shizhe Li
 *  \date    Jul/10/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPOWFLOW_HEADER__
#define __OCPOWFLOW_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPSATFunc.hpp"

#include <vector>

using namespace std;

/// oil-water flow vars suite
class OCPOWFVarSet
{
    /// oil phase is reference phase
public:
    OCPOWFVarSet() { Init0(); }
    void Init0() {
        Swco     = 0;
        krocw    = 0;
        So       = 0; Sw       = 0;
        kro      = 0; krw      = 0;
        dKrodSo  = 0; dKrodSw  = 0;
        dKrwdSo  = 0; dKrwdSw  = 0;
        Pcwo     = 0;
        dPcwodSo = 0; dPcwodSw = 0;
    }

public:
    /// saturaion of connate water
    OCP_DBL Swco;
    /// oil relative permeability in the presence of connate water only
    OCP_DBL krocw;
    /// oil, gas, water saturations
    OCP_DBL So, Sw;
    /// oil, gas, water relatve permeability
    OCP_DBL kro, krw;
    /// the corresponding derivatives of permeability
    OCP_DBL dKrodSo, dKrodSw;
    OCP_DBL dKrwdSo, dKrwdSw;

    /// Capillary pressure   
    OCP_DBL Pcwo;               ///< Pw - Po
    /// the corresponding derivatives of capillary pressure
    OCP_DBL dPcwodSo, dPcwodSw;
};



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
    OCPOWFVarSet* vs;
};


/////////////////////////////////////////////////////
// OCPOWFMethod01
/////////////////////////////////////////////////////


/// Use SWOF
class OCPOWFMethod01 : public OCPOWFMethod
{
public:
    OCPOWFMethod01() = default;
    OCPOWFMethod01(const vector<vector<OCP_DBL>>& SWOFin, OCPOWFVarSet* vsin);
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
// OCPOWFlow
/////////////////////////////////////////////////////

class OCPOWFlow
{
public:
    OCPOWFlow() = default;
    void Setup(const ParamReservoir& rs_param, const USI& i);
    OCPOWFVarSet& GetVarSet() { return vs; }
    void CalKrPc(const OCP_DBL& So, const OCP_DBL& Sw) {
        SetSaturation(So, Sw);
        pfMethod->CalKrPc();
    }
    void CalKrPcDer(const OCP_DBL& So, const OCP_DBL& Sw) {
        SetSaturation(So, Sw);
        pfMethod->CalKrPcDer();
    }

    OCP_DBL GetSwco() const { return pfMethod->GetSwco(); }
    OCP_DBL GetMaxPcow() const { return pfMethod->GetMaxPcow(); }
    OCP_DBL GetMinPcow() const { return pfMethod->GetMinPcow(); }

    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const { return pfMethod->CalPcowBySw(Sw); }
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const { return pfMethod->CalSwByPcow(Pcow); }

protected:
    void SetSaturation(const OCP_DBL& So, const OCP_DBL& Sw) {
        vs.So = So;
        vs.Sw = Sw;
    }

protected:
    OCPOWFVarSet   vs;
    OCPOWFMethod* pfMethod;
};


#endif /* end if __OCPOWFLOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/