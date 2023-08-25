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


class MisFacVarSet
{
public:
    void SetNb(const OCP_USI& nbin) { nb = nbin; }
public:
    /// number of bulks
    OCP_USI         nb;
    /// miscibility factor for permeability
    vector<OCP_DBL> Fk;
    /// miscible factor for capillary pressure
    vector<OCP_DBL> Fp;
};


class MisFacMethod
{
public:
    MisFacMethod() = default;
    virtual void CalculateMiscibleFactor(const OCP_USI& bId, const SurTenVarSet& stvs, MisFacVarSet& mfvs) const = 0;
};


/// Coats expression
class MisFacMethod01 : public MisFacMethod
{
public:
    MisFacMethod01(const Miscstr& param, MisFacVarSet& mfvs);
    void CalculateMiscibleFactor(const OCP_USI& bId, const SurTenVarSet& stvs, MisFacVarSet& mfvs) const override;

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
            mfMethod[mIndex]->CalculateMiscibleFactor(bId, surTen->GetVS(), vs);
        }       
    }
    const auto& GetVarSet() const { return vs; }
    const auto IfUse()const { return ifUse; }
    const auto GetFk(const OCP_USI& bId) const { return vs.Fk[bId]; }
    const auto GetFp(const OCP_USI& bId) const { return vs.Fp[bId]; }
    void ResetTolastTimeStep() { }
    void UpdateLastTimeStep() { }

protected:
    /// if calculate miscible factor
    OCP_BOOL              ifUse{ OCP_FALSE };
    /// miscible factor calculation variable set
    MisFacVarSet          vs;
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
    virtual void CurveCorrect(const OCP_USI& bId, const MisFacVarSet& mfvs) = 0;
    virtual void CurveCorrectDer(const OCP_USI& bId, const MisFacVarSet& mfvs) = 0;
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
    void CurveCorrect(const OCP_USI& bId, const MisFacVarSet& mfvs) override;
    void CurveCorrectDer(const OCP_USI& bId, const MisFacVarSet& mfvs) override;

protected:
    // Dependent modules
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
            mcMethod[mIndex]->CurveCorrect(bId, misFac->GetVarSet());
        }
    }
    void CorrectCurveDer(const OCP_USI& bId, const USI& mIndex) {
        if (ifUse) {
            mcMethod[mIndex]->CurveCorrectDer(bId, misFac->GetVarSet());
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