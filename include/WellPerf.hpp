/*! \file    WellPerf.hpp
 *  \brief   WellPerf class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PERFORATION_HEADER__
#define __PERFORATION_HEADER__

// Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "WellOpt.hpp"

using namespace std;


enum class PerfDirection : USI
{
    /// x-direction
    x,
    /// y-direction
    y,
    /// z-direction
    z,
    /// unstructral grid
    usg
};


/// Perforation describe the connections between wells and bulks.
class Perforation
{

public:
    /// Default constructor.
    Perforation() = default;

    /// Set state of perf
    void SetState(const WellState& flag) { state = flag; };
    /// Return the location of perf: index of bulk
    OCP_USI Location() const { return location; }

public:
    /// state of perforation
    WellState state{ WellState::close };
    OCP_USI  location; ///< Index of bulks which connects to current perforation.
    OCP_DBL  depth;    ///< Depth of bulks which connects to current perforation.
    OCP_DBL  P;        ///< Pressure in current perforation.

    OCP_DBL WI;     ///< Connection transmissibility factor, it can be provided directly
                    ///< from the users.
    OCP_DBL       radius; ///< Well radius.
    OCP_DBL       kh;     ///< Effective permeability times net thickness of the connection.
    OCP_DBL       skinFactor; ///< Skin factor.
    PerfDirection direction;  ///< Direction of the well penetrating the grid block

    /// Multiplier factor for transmissibility of current perforation.
    /// It equals to 0 (close) or 1 (open) now.
    OCP_DBL multiplier;
    /// Molar density of fluid in current perforation. It's used in injection well,
    /// where the fluid consists only single phase.
    mutable OCP_DBL xi;
    vector<OCP_DBL> qi_lbmol; ///< Flow rate of moles of components from into/out
                              ///< current perforation.
    vector<OCP_DBL> transj;   ///< Transmissibility of phase in current perforation.
    OCP_DBL         transINJ;
    vector<OCP_DBL> qj_ft3; ///< Flow rate of volume of phase from into/out current
                            ///< perforation.
    OCP_DBL qt_ft3;         ///< Flow rate of volume of fluids from into/out current
                            ///< perforation.
};

#endif /* end if __PERFORATION_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/