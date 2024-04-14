/*! \file    Solver.hpp
 *  \brief   Solver class declaration
 *  \author  Shizhe Li
 *  \date    Oct/21/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// OpenCAEPoroX header files
#include "IsoThermalSolver.hpp"
#include "ThermalSolver.hpp"
#include "OCPOutput.hpp"

#ifndef __SOLVER_HEADER__
#define __SOLVER_HEADER__

/// Solver class for overall solution methods.
class Solver
{
public:
    /// Setup Solver
    void Setup(Reservoir& rs, const OCPControl& ctrl);
    /// Initialize the reservoir.
    void InitReservoir(Reservoir& rs);
    /// General API
    const OCPNRsuite& GoOneStep(Reservoir& rs, OCPControl& ctrl);

protected:
    /// Model
    OCPModel         model{ OCPModel::none };
    /// Solver for isothermal models with fixed T
    IsothermalSolver IsoTSolver;
    /// Solver for ifThermal models with varied T
    ThermalSolver    TSolver;
};

#endif /* end if __SOLVER_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Shizhe Li           Oct/21/2021      Change from OCPMethod to Solver      */
/*  Chensong Zhang      Oct/27/2021      Rearrange and add comments           */
/*----------------------------------------------------------------------------*/