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

#include <vector>

using namespace std;

/////////////////////////////////////////////////////////////////////
// Miscible For Compositional Model
/////////////////////////////////////////////////////////////////////

class SurfaceTensionMethodParams
{
    friend class SurfaceTensionMethod01;
public:
    /// Default constructor
    SurfaceTensionMethodParams() = default;
    SurfaceTensionMethodParams(
        const vector<OCP_DBL>& parachorin,
        const OCP_DBL*         xvin,
        const OCP_DBL*         xlin,
        const OCP_DBL*         bvin,
        const OCP_DBL*         blin,
        const USI*             NPin
    ) 
    {
        parachor = parachorin;
        xv       = xvin;
        xl       = xlin;
        bv       = bvin;
        bl       = blin;
        NP       = NPin;
    }
protected:
    vector<OCP_DBL> parachor;  ///< used to calculate oil-gas surface tension by Macleod-Sugden correlation
    const OCP_DBL*  xv;        ///< molar fractions of vapour phase
    const OCP_DBL*  xl;        ///< molar fractions of liquid phase
    const OCP_DBL*  bv;        ///< molar density of vapour phase
    const OCP_DBL*  bl;        ///< molar density of liquid phase
    const USI*      NP;        ///< num of present phases
};


class SurfaceTensionMethod
{
public:
    /// Default constructor
    SurfaceTensionMethod() = default;
    /// Calculate surface tensions
    virtual OCP_DBL CalSurfaceTension() = 0;

protected:
    const USI* NP;     ///< num of present phases
};

/// Macleod - Sugden correlation
class SurfaceTensionMethod01 : public SurfaceTensionMethod
{
public:
    /// Default constructor
    SurfaceTensionMethod01() = default;
    /// Construct with parachor
    SurfaceTensionMethod01(const SurfaceTensionMethodParams& params)
    {
        parachor = params.parachor;
        xv       = params.xv;
        xl       = params.xl;
        bv       = params.bv;
        bl       = params.bl;
        NP       = params.NP;

        NC       = parachor.size();
    }
    OCP_DBL CalSurfaceTension();
    vector<OCP_DBL> parachor; ///< used to calculate oil-gas surface tension by Macleod-Sugden correlation
    const OCP_DBL* xv;        ///< molar fractions of vapour phase
    const OCP_DBL* xl;        ///< molar fractions of liquid phase
    const OCP_DBL* bv;        ///< molar density of vapour phase
    const OCP_DBL* bl;        ///< molar density of liquid phase
    USI            NC;        ///< num of components
};


class MiscibleMethod
{
public:
    MiscibleMethod() = default;
};

class Miscible
{
public:
    /// Input param from input file
    void InputParam(const Miscstr& misterm);
    /// Allocate memory for Miscible term
    USI Setup(const OCP_USI& numBulk, const SurfaceTensionMethodParams& params);
    /// Calculate SurfaceTension
    void CalSurfaceTension(const OCP_USI& bId, const USI& mIndex);
    /// Assign value to surTen
    void AssignValue(const OCP_USI& n, const OCP_DBL& v) { surTen[n] = v; }
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
    /// The reference surface tension - flow is immiscible when the surface tension
    //  is greater than or equal to this value
    OCP_DBL surTenRef;
    OCP_DBL surTenPc; ///< Maximum surface tension for capillary pressure / surTenRef
    OCP_DBL Fkexp;    ///< Exponent set used to calculate Fk
    vector<OCP_DBL> surTen; ///< Surface tensions between hydrocarbon phases.
    vector<OCP_DBL> Fk;     ///< The relative permeability interpolation parameter
    vector<OCP_DBL> Fp;     ///< The capillary pressure interpolation parameter.

    // Last time step
    vector<OCP_DBL> lsurTen; ///< last surTen.

    /// Method
protected:
    vector<SurfaceTensionMethod*> surfaceTensionMethod;
};

#endif /* end if __PHASEPERMEABILITY_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/26/2022      Create file                          */
/*  Shizhe Li           Dec/26/2022      Rename to OCPMiscible                */
/*----------------------------------------------------------------------------*/