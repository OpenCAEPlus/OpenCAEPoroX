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
#include "OCPFlowOGW.hpp"

#include <vector>


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


/// Method uesd to correct permeability and capillary pressure
class MisCurveMethod
{
public:
    MisCurveMethod() = default;
    virtual void CurveCorrect(const OCP_DBL& Fk, const OCP_DBL& Fp) = 0;
    virtual void CurveCorrectDer(const OCP_DBL& Fk, const OCP_DBL& Fp) = 0;
};


/// From CMG, see *SIGMA
class MisCurveMethod01 : public MisCurveMethod
{
public:
    MisCurveMethod01() = default;
    MisCurveMethod01(OCPFlowOGW* pf3) {
        OGWF = pf3;
    }
    void CurveCorrect(const OCP_DBL& Fk, const OCP_DBL& Fp) override;
    void CurveCorrectDer(const OCP_DBL& Fk, const OCP_DBL& Fp) override;

protected:
    OCPFlowOGW* OGWF;
};


/// For miscible curve correction
class MiscibleCurve
{
public:
    MiscibleCurve() = default;
    USI Setup(OCPFlowOGW* pf3);
    void CorrectCurve(const USI& mIndex, const OCP_DBL& Fk, const OCP_DBL& Fp) {
        mcMethod[mIndex]->CurveCorrect(Fk, Fp);
    }
    void CorrectCurveDer(const USI& mIndex, const OCP_DBL& Fk, const OCP_DBL& Fp) {
        mcMethod[mIndex]->CurveCorrectDer(Fk, Fp);
    }
protected:
    vector<MisCurveMethod*> mcMethod;
};


/// The class Miscible considers the effect of miscible, Now three main parts are contained:
/// surface tension calculations, miscible factor calculartions, permeability-curve correction.
class Miscible
{
public:
    /// Input param from input file
    void InputParam(const OCP_BOOL& ifmiscible);
    /// Setup for miscible factor calculation
    USI Setup(const OCP_USI& numBulk, const SurTenMethodParams& stparams, const MisFacMethodParams& mfparams);
    /// Setup for miscible curve correction
    USI Setup(OCPFlowOGW* pf3);
    /// Calculate Misscible Factor
    void CalMiscibleFactor(const OCP_USI& bId, const USI& mIndex);
    /// Correct miscible curve
    void CorrectCurve(const OCP_USI& bId, const USI& mIndex) {
        if (ifMiscible) mC.CorrectCurve(mIndex, Fk[bId], Fp[bId]);
    }
    void CorrectCurveDer(const OCP_USI& bId, const USI& mIndex) {
        if (ifMiscible) mC.CorrectCurveDer(mIndex, Fk[bId], Fp[bId]);
    }
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
    /// for miscible curve correction
    MiscibleCurve  mC;

};

#endif /* end if __OCPMISCIBLE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/26/2022      Create file                          */
/*  Shizhe Li           Jul/02/2022      Rename to OCPMiscible                */
/*----------------------------------------------------------------------------*/