/*! \file    Well.hpp
 *  \brief   Well class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __WELL_HEADER__
#define __WELL_HEADER__

// Standard header files
#include <cassert>

// OpenCAEPoroX header files
#include "Bulk.hpp"
#include "DenseMat.hpp"
#include "LinearSystem.hpp"
#include "OCPConst.hpp"
#include "OCPStructure.hpp"
#include "ParamWell.hpp"
#include "WellOpt.hpp"
#include "WellPerf.hpp"
#include "OCPMixture.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////
// General Well Error Type
/////////////////////////////////////////////////////////////////////

const int WELL_SUCCESS           = 10;
const int WELL_NEGATIVE_PRESSURE = -11;
const int WELL_SWITCH_TO_BHPMODE = -12;
const int WELL_CROSSFLOW         = -13;

/// Well class defines well, and any operations referred to wells are in it.
/// Well connects to the bulks by perforations, which serve as source and sink.
/// Due to practical difficulties in production, a good treatment for well is important,
/// excellent treatment will make the flow rate in well more stable.
class Well
{
    friend class AllWells;
    friend class Out4RPT;

public:
    Well() = default;


    /////////////////////////////////////////////////////////////////////
    // Input Param and Setup
    /////////////////////////////////////////////////////////////////////

public:
    /// Input the param of perforations.
    void InputPerfo(const WellParam& well, const Domain& domain, const USI& wId);
    /// Setup the well after Grid and Bulk finish setup.
    void Setup(const Bulk& bk, const vector<SolventINJ>& sols);

protected:
    void SetupUnit();
    /// reservoir unit -> surface unit
    OCP_DBL UnitConvertR2S(const PhaseType& pt, const OCP_DBL& val) const;
    OCP_DBL UnitConvertR2S(const USI& j, const OCP_DBL& val) const;
protected:
    /// Unit Convert (m3,ft3 -> m3,stb,Mscf)
    vector<OCP_DBL> unitConvert;

protected:
    /// Setup well operations
    void SetupOpts(const vector<SolventINJ>& sols);
    void SetupOptsInj(WellOpt& opt, const vector<SolventINJ>& sols);
    void SetupOptsProd(WellOpt& opt);
    

    /////////////////////////////////////////////////////////////////////
    // Basic Well information
    /////////////////////////////////////////////////////////////////////

protected:
    string  name;  ///< well name
    string  group; ///< group well belongs to, it should be moved to opt if necessary!!!
    USI     I;     ///< I-index of the well header.
    USI     J;     ///< J-index of the well header.
    OCP_DBL depth; ///< reference depth of well.
    WellOpt opt;   ///< well control parameters, contains current control parameters.
    vector<WellOpt> optSet;      ///< well control parameters set, contains control
                                 ///< parameters in all critical time.
    USI                 numPerf; ///< num of perforations belonging to this well.
    vector<Perforation> perf;    ///< information of perforation belonging to this well.

    /// num of phases
    USI np;
    /// num of components
    USI nc;

    /// Well surface pressure, psia
    OCP_DBL Psurf{PRESSURE_STD};
    /// Well surface temperature, F
    OCP_DBL Tsurf{TEMPERATURE_STD};

    USI wOId; ///< well index in allWells, closed well is excluded, it's the well index
              ///< in equations

    /// fluid model(where is it from?)
    OCPMixture* mixture;

public:
    /// Initialize the Well BHP
    void InitBHP(const Bulk& myBulk) { bhp = myBulk.P[perf[0].location]; }
    /// Calculate Well Index with Peaceman model.
    void CalWI_Peaceman(const Bulk& myBulk);
    /// Calculate transmissibility for each phase in perforations.
    void CalTrans(const Bulk& myBulk);
    /// Calculate the flux for each perforations.
    void CalFlux(const Bulk& myBulk, const OCP_BOOL ReCalXi = OCP_FALSE);
    /// calculate flow rate of moles of phases for injection well with maxBHP.
    OCP_DBL CalInjRateMaxBHP(const Bulk& myBulk);
    /// calculate flow rate of moles of phases for production well with minBHP.
    OCP_DBL CalProdRateMinBHP(const Bulk& myBulk);
    /// Calculate flow rate of moles of phases for injection well with calculated
    /// qi_lbmol.
    void CalInjQj(const Bulk& myBulk, const OCP_DBL& dt);
    /// Calculate flow rate of moles of phases for production well with calculated
    /// qi_lbmol.
    void CalProdQj(const Bulk& myBulk, const OCP_DBL& dt);
    /// Calculate pressure difference between well and perforations.
    void CaldG(const Bulk& myBulk);
    /// Calculate pressure difference between well and perforations for Injection.
    void CalInjdG(const Bulk& myBulk);
    /// Calculate pressure difference between well and perforations for Production.
    void CalProddG(const Bulk& myBulk);
    /// Calculate pressure difference between well and perforations for Production.
    void CalProddG01(const Bulk& myBulk);
    /// Calculate pressure difference between well and perforations for Production.
    void CalProddG02(const Bulk& myBulk);
    /// Calculate the production weight
    void CalFactor(const Bulk& myBulk) const;
    /// Correct BHP if opt mode is BHPMode
    void CorrectBHP();
    /// Check if well operation mode would be changed.
    void CheckOptMode(const Bulk& myBulk);
    /// Check if abnormal Pressure occurs.
    OCP_INT CheckP(const Bulk& myBulk);
    /// Check if cross flow happens.
    OCP_INT CheckCrossFlow(const Bulk& myBulk);
    /// Update pressure in Perforation after well pressure updates.
    void CalPerfP()
    {
        for (USI p = 0; p < numPerf; p++) perf[p].P = bhp + dG[p];
    }
    /// Display operation mode of well and state of perforations.
    void ShowPerfStatus(const Bulk& myBulk) const;

    USI     PerfNum() const { return numPerf; }
    void    SetBHP(const OCP_DBL& p) { bhp = p; }
    OCP_DBL BHP() const { return bhp; }

    /// Return the state of the well, Open or Close.
    OCP_BOOL IsOpen() const { return (opt.state == WellState::open); }
    /// Return the type of well, Inj or Prod.
    auto WellType() const { return opt.type; }
    auto PerfState(const USI& p) const { return perf[p].state; }
    USI      PerfLocation(const USI& p) const { return perf[p].location; }
    OCP_DBL  PerfQi_lbmol(const USI& p, const USI& i) const
    {
        return perf[p].qi_lbmol[i];
    }
    OCP_DBL PerfProdQj_ft3(const USI& p, const USI& j) const
    {
        OCP_ASSERT(opt.type == WellType::productor, "Wrong Call");
        return perf[p].qj_ft3[j];
    }

protected:
    
    /////////////////////////////////////////////////////////////////////
    // Well Physical information
    /////////////////////////////////////////////////////////////////////

    mutable OCP_DBL bhp; ///< Well pressure in reference depth.
    vector<OCP_DBL>
        dG; ///< difference of pressure between well and perforation: numPerf.

    // Last time step
    OCP_DBL         lbhp; ///< Last BHP
    vector<OCP_DBL> ldG;  ///< Last dG

    // PROD/INJ Rate
    vector<OCP_DBL> qi_lbmol; ///< flow rate of moles of component inflowing/outflowing
    /// components mole number -> target phase volume
    mutable vector<OCP_DBL> factor;


    OCP_DBL WOPR{0}; ///< well oil production rate.
    OCP_DBL WOPT{0}; ///< well total oil production.
    OCP_DBL WGPR{0}; ///< well gas production rate.
    OCP_DBL WGPT{0}; ///< well total gas production.
    OCP_DBL WWPR{0}; ///< well water production rate.
    OCP_DBL WWPT{0}; ///< well total water production.
    OCP_DBL WGIR{0}; ///< well gas injection rate.
    OCP_DBL WGIT{0}; ///< well total gas injection.
    OCP_DBL WWIR{0}; ///< well water injection rate.
    OCP_DBL WWIT{0}; ///< well total water injection.

    OCP_BOOL ifUseUnweight{OCP_FALSE};
    /////////////////////////////////////////////////////////////////////
    // IMPEC
    /////////////////////////////////////////////////////////////////////

public:

    /////////////////////////////////////////////////////////////////////
    // FIM
    /////////////////////////////////////////////////////////////////////

public:
    void AssembleMatFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
    void CalResFIM(OCP_USI& wId, OCPRes& res, const Bulk& bk, const OCP_DBL& dt) const;
protected:
    void AssembleMatInjFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
    void AssembleMatProdFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;

public:
    void AssembleMatFIM_T(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
    void CalResFIM_T(OCP_USI& wId, OCPRes& res, const Bulk& bk, const OCP_DBL& dt) const;
protected:
    void AssembleMatInjFIM_T(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
    void AssembleMatProdFIM_T(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;

public:
    void AssembleMatIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
protected:
    void AssembleMatInjIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
    void AssembleMatProdIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;

    // for output
    //void SetPolyhedronWell(const Grid& myGrid, OCPpolyhedron& mypol);
};

#endif /* end if __WELL_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/