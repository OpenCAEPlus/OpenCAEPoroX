/*! \file    OCPOGFlow.hpp
 *  \brief   OCPOGFlow class declaration
 *  \author  Shizhe Li
 *  \date    Jul/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPOGFLOW_HEADER__
#define __OCPOGFLOW_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPSATFunc.hpp"

#include <vector>

using namespace std;

/// oil-water flow vars suite
class OCPOGFVarSet
{
    /// oil phase is reference phase
public:
    OCPOGFVarSet() { Init0(); }
    void Init0() {
        So       = 0; Sg       = 0;
        kro      = 0; krg      = 0;
        dKrodSo  = 0; dKrodSg  = 0;
        dKrgdSo  = 0; dKrgdSg  = 0;
        Pcgo     = 0;
        dPcgodSo = 0; dPcgodSg = 0;
    }

public:
    /// oil, gas saturations
    OCP_DBL So, Sg;
    /// oil, gas relatve permeability
    OCP_DBL kro, krg;
    /// the corresponding derivatives of permeability
    OCP_DBL dKrodSo, dKrodSg;
    OCP_DBL dKrgdSo, dKrgdSg;

    /// Capillary pressure   
    OCP_DBL Pcgo;               ///< Pg - Po
    /// the corresponding derivatives of capillary pressure
    OCP_DBL dPcgodSo, dPcgodSg;
};



/// Calculate oil, gas, water relative permeability and capillary pressure
class OCPOGFMethod
{
public:
    OCPOGFMethod() = default;
    virtual void CalKrPc() = 0;
    virtual void CalKrPcDer() = 0;
    virtual OCP_DBL CalPcgoBySg(const OCP_DBL& sg) const = 0;
    virtual OCP_DBL CalSgByPcgo(const OCP_DBL& pcgo) const = 0;

protected:
    OCPOGFVarSet* vs;
};


/////////////////////////////////////////////////////
// OCPOGFMethod01
/////////////////////////////////////////////////////


/// Use SGOF
class OCPOGFMethod01 : public OCPOGFMethod
{
public:
    OCPOGFMethod01(const vector<vector<OCP_DBL>>& SGOFin, OCPOGFVarSet* vsin);
    void CalKrPc() override;
    void CalKrPcDer() override;

    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const override { return SGOF.CalPcgo(Sg); }
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const override { return SGOF.CalSg(Pcgo); }

protected:
    OCP_SGOF            SGOF;
};


/////////////////////////////////////////////////////
// OCPOGFlow
/////////////////////////////////////////////////////

class OCPOGFlow
{
public:
    OCPOGFlow() = default;
    void Setup(const ParamReservoir& rs_param, const USI& i);
    OCPOGFVarSet& GetVarSet() { return vs; }
    void CalKrPc(const OCP_DBL& So, const OCP_DBL& Sg) {
        SetSaturation(So, Sg);
        pfMethod->CalKrPc();
    }
    void CalKrPcDer(const OCP_DBL& So, const OCP_DBL& Sg) {
        SetSaturation(So, Sg);
        pfMethod->CalKrPcDer();
    }

    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const { return pfMethod->CalSgByPcgo(Pcgo); }
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const { return pfMethod->CalPcgoBySg(Sg); }

protected:
    void SetSaturation(const OCP_DBL& So, const OCP_DBL& Sg) {
        vs.So = So;
        vs.Sg = Sg;
    }

protected:
    OCPOGFVarSet   vs;
    OCPOGFMethod*  pfMethod;
};


#endif /* end if __OCPOGFLOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/