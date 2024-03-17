/*! \file    ThermalSolver.cpp
 *  \brief   ThermalSolver class declaration
 *  \author  Shizhe Li
 *  \date    Nov/10/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// OpenCAEPoroX header files
#include "ThermalSolver.hpp"

void ThermalSolver::SetupMethod(Reservoir& rs, const OCPControl& ctrl)
{
    fim.Setup(rs, ctrl);
    fim.SetWorkLS(LSolver.Setup(ctrl.SM.GetModel(), ctrl.SM.GetWorkDir(), ctrl.SM.GetLsFile(0), rs.GetDomain(), rs.GetComNum() + 2), 0);
}

void ThermalSolver::InitReservoir(Reservoir& rs) { fim.InitReservoir(rs); }

void ThermalSolver::Prepare(Reservoir& rs, const OCPControl& ctrl)
{
    fim.Prepare(rs, ctrl);
}


const OCPNRsuite& ThermalSolver::GoOneStep(Reservoir& rs, OCPControl& ctrl)
{
    // Prepare for time marching
    Prepare(rs, ctrl);

    // Time marching with adaptive time stepsize
    while (OCP_TRUE) {
        if (ctrl.time.GetCurrentDt() < MIN_TIME_CURSTEP) {
            if (CURRENT_RANK == MASTER_PROCESS)
                OCP_WARNING("Time stepsize is too small: " + to_string(ctrl.time.GetCurrentDt()) + TIMEUNIT);
            ctrl.StopSim = OCP_TRUE;
            break;
        }
        // Assemble linear system
        AssembleMat(rs, ctrl);
        // Solve linear system
        SolveLinearSystem(rs, ctrl);
        if (!UpdateProperty(rs, ctrl)) {
            continue;
        }
        if (FinishNR(rs, ctrl)) break;
    }

    // Finish current time step
    FinishStep(rs, ctrl);

    return fim.GetNRsuite();
}


void ThermalSolver::AssembleMat(const Reservoir& rs, OCPControl& ctrl)
{
    GetWallTime timer;
    timer.Start();

    fim.AssembleMat(LSolver, rs, ctrl.time.GetCurrentDt());

    OCPTIME_ASSEMBLE_MAT += timer.Stop() / TIME_S2MS;
}

void ThermalSolver::SolveLinearSystem(Reservoir& rs, OCPControl& ctrl)
{
    fim.SolveLinearSystem(LSolver, rs, ctrl);
}

/// Update properties of fluid.
OCP_BOOL ThermalSolver::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    GetWallTime timer;
    timer.Start();
    OCP_BOOL flag = fim.UpdateProperty(rs, ctrl);

    OCPTIME_UPDATE_GRID += timer.Stop() / TIME_S2MS;
    return flag;
}

/// Finish the Newton-Raphson iteration.
OCP_BOOL ThermalSolver::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    return fim.FinishNR(rs, ctrl);
}

/// Finish the current time step.
void ThermalSolver::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    fim.FinishStep(rs, ctrl);
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/