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


/// Coats expression
class MisFacMethod01 : public MisFacMethod
{
public:
    MisFacMethod01(const Miscstr& param);
    void CalculateMiscibleFactor(const OCP_DBL& st, OCP_DBL& fk, OCP_DBL& fp) const override;

public:
    OCP_DBL stref;   ///< reference surface tension
    OCP_DBL fkExp;   ///< exponent for Fk calcualtion
    OCP_DBL stPcf;   ///< maximum surface tension for capillary pressure calculation / reference surface tension
};


class MiscibleFactor
{
public:
    /// Default constructor
    MiscibleFactor() = default;
    USI Setup(const ParamReservoir& param, const USI& i, const OCP_USI& nb, const SurfaceTension* st);
    void CalMiscibleFactor(const OCP_USI& bId, const USI& mIndex) {
        if (ifUse) {
            mfMethod[mIndex]->CalculateMiscibleFactor(surTen->GetSurfaceTension(bId), Fk[bId], Fp[bId]);
        }       
    }
    const auto IfUse()const { return ifUse; }
    const auto GetFk(const OCP_USI& bId) const { return Fk[bId]; }
    const auto GetFp(const OCP_USI& bId) const { return Fp[bId]; }
    void ResetTolastTimeStep() { }
    void UpdateLastTimeStep() { }

protected:
    /// if calculate miscible factor
    OCP_BOOL              ifUse{ OCP_FALSE };
    /// miscibility factor for permeability
    vector<OCP_DBL>       Fk;
    /// miscible factor for capillary pressure
    vector<OCP_DBL>       Fp;
    /// method of miscible factor calculation
    vector<MisFacMethod*> mfMethod;

    // Dependent modules
    /// Surface tension Calculation
    const SurfaceTension* surTen;
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
    MisCurveMethod01(OCPFlow* flowin) {
        switch (flowin->FlowType())
        {
        case OCPFlowType::OG:
        case OCPFlowType::OGW:
            break;
        default:
            OCP_ABORT("Wrong FlowType for ScalePcow!");
            break;
        }

        flow = flowin;
    }
    void CurveCorrect(const OCP_DBL& Fk, const OCP_DBL& Fp) override;
    void CurveCorrectDer(const OCP_DBL& Fk, const OCP_DBL& Fp) override;

protected:
    OCPFlow* flow;
};


/// For miscible curve correction
class MiscibleCurve
{
public:
    MiscibleCurve() = default;
    USI Setup(OCPFlow* flowin, const MiscibleFactor* misfactor);
    void CorrectCurve(const OCP_USI& bId, const USI& mIndex) {
        if (ifUse) {
            mcMethod[mIndex]->CurveCorrect(misFac->GetFk(bId), misFac->GetFp(bId));
        }
    }
    void CorrectCurveDer(const OCP_USI& bId, const USI& mIndex) {
        if (ifUse) {
            mcMethod[mIndex]->CurveCorrectDer(misFac->GetFk(bId), misFac->GetFp(bId));
        }
    }
    void ResetTolastTimeStep() { }
    void UpdateLastTimeStep() { }
protected:
    /// If use Miscible curve correction
    OCP_BOOL                ifUse{ OCP_FALSE };
    /// Methods for Miscible Curve 
    vector<MisCurveMethod*> mcMethod;

    // Dependent modules
    const MiscibleFactor*   misFac;
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