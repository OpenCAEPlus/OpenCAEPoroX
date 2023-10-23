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
#include "ParamWell.hpp"
#include "WellOpt.hpp"
#include "WellPerf.hpp"
#include "OCPMixture.hpp"
#include "OCPNRresidual.hpp"

using namespace std;


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
    virtual void InputPerfo(const WellParam& well, const Domain& domain, const USI& wId) = 0;
    /// Setup the well after Grid and Bulk finish setup.
    virtual void Setup(const Bulk& bk, const vector<SolventINJ>& sols) = 0;
    /// Apply ith operations
    OCP_BOOL ApplyOpt(const USI& i);

protected:
    /// Setup well operations
    void SetupOpts(const vector<SolventINJ>& sols);
    /// Setup well operations for Injection well
    void SetupOptsInj(WellOpt& opt, const vector<SolventINJ>& sols);
    /// Setup well operations for Production well
    void SetupOptsProd(WellOpt& opt);


public:
    /// Initialize the Well Pressure
    virtual void InitWellP(const Bulk& bk) = 0;
    /// Check if well operation mode would be changed.
    virtual void CheckOptMode(const Bulk& bk) = 0;
    /// Calculate Flux and initialize values which will not be changed during this time step
    virtual void CalFluxInit(const Bulk& bk) = 0;
    /// Calculate Flux
    virtual void CalFlux(const Bulk& bk) = 0;
    /// Check if abnormal Pressure occurs.
    virtual ReservoirState CheckP(const Bulk& bk) = 0;
    /// Calculate flow rate of moles of phases for injection well and production well
    virtual void CalIPRate(const Bulk& bk, const OCP_DBL& dt) = 0;
    /// Calculate max change of well pressure between two time step
    virtual OCP_DBL CalMaxChangeTime() const = 0;
    /// Calculate max change of well pressure between two NR step
    virtual OCP_DBL CalMaxChangeNR() = 0;
    /// Reset to last time step
    virtual void ResetToLastTimeStep(const Bulk& bk) = 0;
    /// Update last time step
    virtual void UpdateLastTimeStep() = 0;


public:
    /// Display operation mode of well and state of perforations.
    void     ShowPerfStatus(const Bulk& bk) const;
    /// Get num of peferations
    USI      PerfNum() const { return numPerf; }
    /// Return the state of the well, Open or Close.
    OCP_BOOL IsOpen() const { return (opt.state == WellState::open); }
    /// Return the type of well, Inj or Prod.
    auto     WellType() const { return opt.type; }
    /// Return the state of perforations
    auto     PerfState(const USI& p) const { return perf[p].state; }
    /// Return the location(bulk index) of perforations
    USI      PerfLocation(const USI& p) const { return perf[p].location; }
    /// Return mole flow rate of component i from perforation j 
    OCP_DBL  PerfQi_lbmol(const USI& p, const USI& i) const { return perf[p].qi_lbmol[i]; }
    /// Return volume flow rate of component i from perforation j 
    OCP_DBL  PerfProdQj_ft3(const USI& p, const USI& j) const { return perf[p].qj_ft3[j]; }

protected:
    /// well name
    string              name;
    /// belonging group of well, it would be moved to opt if necessary.
    string              group;
    /// reference depth of well.
    OCP_DBL             depth;
    /// Well pressure in reference depth.
    OCP_DBL             bhp;
    /// well control parameters, contains current control parameters.
    WellOpt             opt;
    /// well control parameters set.
    vector<WellOpt>     optSet;
    /// num of perforations.
    USI                 numPerf; 
    /// information of perforations.
    vector<Perforation> perf;    
    /// Well surface pressure, psia
    OCP_DBL             Psurf;
    /// Well surface temperature, F
    OCP_DBL             Tsurf;
    /// mixture model
    OCPMixture*         mixture;
    /// num of phases
    USI                 np;
    /// num of components
    USI                 nc;
    /// reservoir temperature
    OCP_DBL             rsTemp;
    /// flow rate of moles of component inflowing/outflowing
    vector<OCP_DBL>     qi_lbmol;
    /// ifUse unweighted permeability.
    OCP_BOOL            ifUseUnweight{ OCP_FALSE };
    /// bhp at last time step
    OCP_DBL             lbhp;
    /// bhp at last NR step
    OCP_DBL             NRbhp;

    /// well oil production rate.
    OCP_DBL             WOPR{0};
    /// well total oil production.
    OCP_DBL             WOPT{0}; 
    /// well gas production rate.
    OCP_DBL             WGPR{0};
    /// well total gas production.
    OCP_DBL             WGPT{0}; 
    /// well water production rate.
    OCP_DBL             WWPR{0};
    /// well total water production.
    OCP_DBL             WWPT{0}; 
    /// well gas injection rate.
    OCP_DBL             WGIR{0};
    /// well total gas injection.
    OCP_DBL             WGIT{0}; 
    /// well water injection rate.
    OCP_DBL             WWIR{0}; 
    /// well total water injection.
    OCP_DBL             WWIT{0}; 



    /////////////////////////////////////////////////////////////////////
    // Unit
    /////////////////////////////////////////////////////////////////////

protected:
    void SetupUnit();
    /// reservoir unit -> surface unit for specific phase
    OCP_DBL UnitConvertR2S(const PhaseType& pt, const OCP_DBL& val) const;
    /// reservoir unit -> surface unit for specific phase index
    OCP_DBL UnitConvertR2S(const USI& j, const OCP_DBL& val) const;
protected:
    /// Unit Convert (m3,ft3 -> m3,stb,Mscf)
    vector<OCP_DBL> unitConvert;


    /////////////////////////////////////////////////////////////////////
    // Method
    /////////////////////////////////////////////////////////////////////

public:
    /// Calculate resiual for FIM method
    virtual void CalResFIM(OCP_USI& wId, OCPNRresidual& res, const Bulk& bk, const OCP_DBL& dt) const = 0;
    /// Assemble matrix for FIM method
    virtual void AssembleMatFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const = 0;
    /// Get solution for FIM method
    virtual void GetSolutionFIM(const vector<OCP_DBL>& u, OCP_USI& wId) = 0;
    /// Assemble matrix for IMPEC method
    virtual void AssembleMatIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const = 0;
    /// Get solution for IMPEC method
    virtual void GetSolutionIMPEC(const vector<OCP_DBL>& u, OCP_USI& wId) = 0;
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