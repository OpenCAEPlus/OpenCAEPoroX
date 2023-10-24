/*! \file    WellOpt.hpp
 *  \brief   WellOpt class declaration
 *  \author  Shizhe Li
 *  \date    Nov/22/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __WELLOPT_HEADER__
#define __WELLOPT_HEADER__

// Standard header files
#include <cmath>

// OpenCAEPoroX header files
#include "ParamWell.hpp"
#include "OCPMixture.hpp"

using namespace std;

enum class WellType : USI
{
    /// injection well
    injector,
    /// production well
    productor
};

enum class WellState : USI
{
    /// Open
    open,
    /// Close
    close
};

enum class WellOptMode : USI
{
    /// injection well rate control
    irate,
    /// production well oil rate control
    orate,
    /// production well gas rate control
    grate,
    /// production well water rate control
    wrate,
    /// production well liquid rate control
    lrate,
    /// injection/production well bhp control
    bhp,
};

/// WellOpt describes the operation mode of a well.
/// usually it changes over time, specifically, each attributes could be changed
/// including the well type.
class WellOpt
{
    friend class Well;
    friend class AllWells;
    friend class Out4RPT;

public:
    /// Default constructor.
    WellOpt() = default;

    /// Constructor well operation mode using params.
    WellOpt(const WellOptParam& Optparam);

    /// overload inequality
    OCP_BOOL operator!=(const WellOpt& Opt) const;

public:
    /// type of well, Injection well or Production well.
    WellType type;
    /// indicate which type of fluids will be injected, water, gas, or other solvent.
    /// it's decided by users and only useful for injection well.
    string   injFluidName;
    /// state of well
    WellState state{ WellState::close };
    /// Well control mode
    WellOptMode mode;
    /// Initial control mode
    WellOptMode initMode;
    /// it gives the upper limit of flow rate of specified fluids if the well is under
    /// the control of constant pressure. it gives the flow rate of specified fluids if
    /// the well is under the control of constant flow rate.
    OCP_DBL maxRate;
    /// it gives the upper limit of injection well pressure if the well is under the control of
    /// constant flow rate. it gives the pressure of well if the well is under the
    /// control of constant pressure.
    OCP_DBL maxBHP;
    /// it gives the lower limit of production well pressure if the well is under the control of
    /// constant flow rate. it gives the pressure of well if the well is under the
    /// control of constant pressure.
    OCP_DBL minBHP;
    /// target injection/production well rate
    OCP_DBL tarRate;
    /// target injection/production well BHP
    OCP_DBL tarBHP;

    /// for injection well, it describes the components of injected fluids.
    vector<OCP_DBL> injZi;
    /// type of injecting fluid
    PhaseType       injPhase;
    /// production weight of phase (1.0 or 0.0)
    vector<OCP_DBL> prodPhaseWeight;
    /// temperature of injected fluid (unit : F)
    OCP_DBL injTemp;
};

/// Describe the molar fraction of components of fluid injected to reservoir from INJ.
class SolventINJ
{
public:
    SolventINJ() = default;
    SolventINJ(const Solvent& other)
    {
        name = other.name;
        data = other.comRatio;
    };
    string          name; ///< name of solvent
    vector<OCP_DBL> data; ///< molar fraction of components
};

#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           NOV/22/2022      Create file                          */
/*----------------------------------------------------------------------------*/
