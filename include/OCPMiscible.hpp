/*! \file    OCPMiscible.hpp
 *  \brief   OCPMiscible class declaration
 *  \author  Shizhe Li
 *  \date    Dec/26/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPMISCIBLE_HEADER__
#define __OCPMISCIBLE_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPSurfaceTension.hpp"

#include <vector>
#include <functional>

using namespace std;



class MisFacMethod
{
public:
    MisFacMethod() = default;
    virtual void CalculateMiscibleFactor(const OCP_DBL& st, OCP_DBL& fk, OCP_DBL& fp) const = 0;
};


class MisFacMethod01Params
{
    friend class MisFacMethod01;
public:
    MisFacMethod01Params() = default;
    MisFacMethod01Params(const OCP_DBL& surTenRef, const OCP_DBL& FkExp, const OCP_DBL& surTenPc) {
        ifUse   = OCP_TRUE;
        stref   = surTenRef;
        fkExp   = FkExp;
        stPc    = surTenPc;
        if (stPc < 0) stPc = stref;
    }
    OCP_BOOL IfUse() const { return ifUse; }

protected:
    OCP_BOOL ifUse{ OCP_FALSE };
    OCP_DBL  stref; ///< The reference surface tension (maximum allowed miscible surface tension)
    OCP_DBL  fkExp; ///< exponent of the surface tension calculation
    OCP_DBL  stPc;  ///< The maximum surface tension used to scale the input capillary pressure curves
};


/// Coats expression
class MisFacMethod01 : public MisFacMethod
{
public:
    MisFacMethod01() = default;
    MisFacMethod01(const MisFacMethod01Params& param) {
        stref = param.stref;
        fkExp = param.fkExp;
        stPcf = param.stPc / stref;
    }
    void CalculateMiscibleFactor(const OCP_DBL& st, OCP_DBL& fk, OCP_DBL& fp) const override;

public:
    OCP_DBL stref;   ///< reference surface tension
    OCP_DBL fkExp;   ///< exponent for Fk calcualtion
    OCP_DBL stPcf;   ///< maximum surface tension for capillary pressure calculation / reference surface tension
};


class MisFacMethodParams
{
public:
    MisFacMethodParams() = default;
    MisFacMethodParams(const MisFacMethod01Params& param01in) {
        param01 = param01in;
    }
    void SetupParams(const MisFacMethod01Params& param01in) {
        param01 = param01in;
    }

public:
    MisFacMethod01Params param01;
};

class MiscibleFcator
{
public:
    /// Default constructor
    MiscibleFcator() = default;
    USI Setup(const MisFacMethodParams& param);
    void CalculateMiscibleFactor(const USI& mIndex, const OCP_DBL& st, OCP_DBL& fk, OCP_DBL& fp) {
        mfMethod[mIndex]->CalculateMiscibleFactor(st, fk, fp);
    }
protected:
    vector<MisFacMethod*> mfMethod;
};


/// Params uesd to correct permeability and capillary pressure
class MiscibleCorrectionMethodParams
{
    friend class MiscibleCorrectionMethod01;
public:
    /// Default constructor
    MiscibleCorrectionMethodParams() = default;
    MiscibleCorrectionMethodParams(
        OCP_DBL* krin,
        OCP_DBL* Pcin,
        OCP_DBL* dPcdSin,
        const function<OCP_DBL(const OCP_DBL&)>& CalKrgin
    ) {
        kr_out    = krin;
        Pc_out    = Pcin;
        dPcdS_out = dPcdSin;
        CalKrg    = CalKrgin;
    }

protected:
    OCP_DBL*                          kr_out;    ///< permeability (will be corrected)
    OCP_DBL*                          Pc_out;    ///< capillary pressure (will be corrected)
    OCP_DBL*                          dPcdS_out; ///< dPcdS (will be corrected)
    function<OCP_DBL(const OCP_DBL&)> CalKrg;    ///< function to calculate krg
};


/// Method uesd to correct permeability and capillary pressure
class MiscibleCorrectionMethod
{
public:
    MiscibleCorrectionMethod() = default;

protected:

    OCP_DBL* kr_out;                           ///< permeability (will be corrected)
    OCP_DBL* Pc_out;                           ///< capillary pressure (will be corrected)
    OCP_DBL* dPcdS_out;                        ///< dPcdS (will be corrected)
    function<OCP_DBL(const OCP_DBL&)> CalKrg;  ///< function to calculate krg
};


/// from CMG, see *SIGMA
class MiscibleCorrectionMethod01 : public MiscibleCorrectionMethod
{
public:
    MiscibleCorrectionMethod01() = default;
    MiscibleCorrectionMethod01(const MiscibleCorrectionMethodParams& params) {
        kr_out     = params.kr_out;   
        kr_out     = params.kr_out;
        Pc_out     = params.Pc_out;
        dPcdS_out  = params.dPcdS_out;
        CalKrg     = params.CalKrg;
    }
    void CorrectOGCurve();
};


/// The class Miscible considers the effect of miscible, Now three main parts are contained:
/// surface tension calculations, miscible factor calculartions, permeability-curve correction.
class Miscible
{
public:
    /// Input param from input file
    void InputParam(const OCP_BOOL& ifmiscible);
    /// Setup for Miscible terms
    USI Setup(const OCP_USI& numBulk, const SurTenMethodParams& stparams, const MisFacMethodParams& mfparams);
    /// Calculate Misscible Factor
    void CalMiscibleFactor(const OCP_USI& bId, const USI& mIndex);
    /// Calculate Fk, Fp and return if miscible
    OCP_BOOL CalFkFp(const OCP_USI& n, OCP_DBL& fk, OCP_DBL& fp);
    /// Reset Miscible term to last time step
    void ResetTolastTimeStep() { surTen = lsurTen; }
    /// Update Miscible term at last time step
    void UpdateLastTimeStep() { lsurTen = surTen; }
    /// Return surTen
    OCP_DBL GetSurTen(const OCP_USI& n) const { return surTen[n]; }
    /// Return Fk
    OCP_DBL GetFk(const OCP_USI& n) const { return Fk[n]; }
    /// Return Fp
    OCP_DBL GetFp(const OCP_USI& n) const { return Fp[n]; }

protected:
    /// Miscible treatment of hydrocarbons, only used in compositional Model.
    OCP_BOOL ifMiscible{OCP_FALSE};
    vector<OCP_DBL> surTen; ///< Surface tensions between hydrocarbon phases.
    vector<OCP_DBL> Fk;     ///< miscibility factor for permeability
    vector<OCP_DBL> Fp;     ///< miscibility factor for capillary pressure

    // Last time step
    vector<OCP_DBL> lsurTen; ///< last surTen.

    /// Method
protected:
    /// for surface tension calculation
    SurfaceTension sT;
    /// for miscible factor calculation
    MiscibleFcator mF;

};

#endif /* end if __PHASEPERMEABILITY_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/26/2022      Create file                          */
/*  Shizhe Li           Jul/02/2022      Rename to OCPMiscible                */
/*----------------------------------------------------------------------------*/