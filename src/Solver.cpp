/*! \file    Solver.cpp
 *  \brief   Solver class definition
 *  \author  Shizhe Li
 *  \date    Oct/21/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// OpenCAEPoroX header files
#include "Solver.hpp"

void Solver::Setup(Reservoir& rs, const OCPControl& ctrl)
{

    model = ctrl.SM.GetModel();
    switch (model) 
    {
        case OCPModel::isothermal:
            IsoTSolver.SetupMethod(rs, ctrl);
            break;
        case OCPModel::thermal:
            TSolver.SetupMethod(rs, ctrl);
            break;
        default:
            OCP_ABORT("Wrong model type specified!");
            break;
    }
}


/// Initialize the reservoir setting for different solution methods.
void Solver::InitReservoir(Reservoir& rs)
{
    switch (model) 
    {
        case OCPModel::isothermal:
            IsoTSolver.InitReservoir(rs);
            break;
        case OCPModel::thermal:
            TSolver.InitReservoir(rs);
            break;
        default:
            OCP_ABORT("Wrong model type specified!");
            break;
    }
}


/// This is one time step of dynamic simulation in an abstract setting.
const OCPNRsuite& Solver::GoOneStep(Reservoir& rs, OCPControl& ctrl)
{
    switch (model) 
    {
        case OCPModel::isothermal:
            return IsoTSolver.GoOneStep(rs, ctrl);
        case OCPModel::thermal:
            return TSolver.GoOneStep(rs, ctrl);
            break;
        default:
            OCP_ABORT("Wrong model type specified!");
            break;
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/21/2021      Create file                          */
/*----------------------------------------------------------------------------*/