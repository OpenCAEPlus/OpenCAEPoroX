/*! \file    IsothermalSolver.cpp
 *  \brief   IsothermalSolver class definition
 *  \author  Shizhe Li
 *  \date    Oct/21/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "IsothermalSolver.hpp"

/// Setup solution methods, including IMPEC and FIM.
void IsothermalSolver::SetupMethod(Reservoir& rs, const OCPControl& ctrl)
{
    method = ctrl.GetMethod();

    switch (method) {
        case OCPNLMethod::IMPEC:
            impec.Setup(rs, LSolver, ctrl);
            break;
        case OCPNLMethod::FIM:
            fim.Setup(rs, LSolver, ctrl);
            break;
        case OCPNLMethod::AIMc:
            aimc.Setup(rs, LSolver, ctrl);
            break;
        case OCPNLMethod::FIMddm:
            fim_ddm.IsoT_FIM::Setup(rs, LSolver, ctrl);
            break;
        default:
            OCP_ABORT("Wrong method type!");
    }
}


/// Setup solution methods, including IMPEC and FIM.
void IsothermalSolver::InitReservoir(Reservoir& rs)
{
    switch (method) {
        case OCPNLMethod::IMPEC:
            impec.InitReservoir(rs);
            break;
        case OCPNLMethod::FIM:
            fim.InitReservoir(rs);
            break;
        case OCPNLMethod::AIMc:
            aimc.InitReservoir(rs);
            break;
        case OCPNLMethod::FIMddm:
            fim_ddm.IsoT_FIM::InitReservoir(rs);
            break;
        default:
            OCP_ABORT("Wrong method type!");
    }
}


/// Prepare solution methods, including IMPEC and FIM.
void IsothermalSolver::Prepare(Reservoir& rs, OCPControl& ctrl)
{
    switch (method) {
        case OCPNLMethod::IMPEC:
            impec.Prepare(rs, ctrl);
            break;
        case OCPNLMethod::FIM:
            fim.Prepare(rs, ctrl.time.GetCurrentDt());
            break;
        case OCPNLMethod::AIMc:
            aimc.Prepare(rs, ctrl.time.GetCurrentDt());
            break;
        case OCPNLMethod::FIMddm:
            fim_ddm.IsoT_FIM::Prepare(rs, ctrl.time.GetCurrentDt());
            break;
        default:
            OCP_ABORT("Wrong method type!");
    }
}


/// Assemble linear systems for IMPEC and FIM.
void IsothermalSolver::AssembleMat(const Reservoir& rs, OCPControl& ctrl)
{
    // Assemble
    const OCP_DBL dt = ctrl.time.GetCurrentDt();
    GetWallTime timer;
    timer.Start();

    switch (method) {
        case OCPNLMethod::IMPEC:
            impec.AssembleMat(LSolver, rs, dt);
            break;
        case OCPNLMethod::FIM:
            fim.AssembleMat(LSolver, rs, dt);
            break;
        case OCPNLMethod::AIMc:
            aimc.AssembleMat(LSolver, rs, dt);
            break;
        case OCPNLMethod::FIMddm:
            fim_ddm.IsoT_FIM::AssembleMat(LSolver, rs, dt);
            break;
        default:
            OCP_ABORT("Wrong method type!");
    }
   
    OCPTIME_ASSEMBLE_MAT += timer.Stop() / TIME_S2MS;
}


/// Solve linear systems for IMPEC and FIM.
void IsothermalSolver::SolveLinearSystem(Reservoir& rs, OCPControl& ctrl)
{ 
    switch (method) {
        case OCPNLMethod::IMPEC:
            impec.SolveLinearSystem(LSolver, rs, ctrl);
            break;
        case OCPNLMethod::FIM:
            fim.SolveLinearSystem(LSolver, rs, ctrl);
            break;
        case OCPNLMethod::AIMc:
            aimc.SolveLinearSystem(LSolver, rs, ctrl);
            break;
        case OCPNLMethod::FIMddm:
            fim_ddm.IsoT_FIM::SolveLinearSystem(LSolver, rs, ctrl);
            break;
        default:
            OCP_ABORT("Wrong method type!");
    }
}


/// Update physical properties for IMPEC and FIM.
OCP_BOOL IsothermalSolver::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_BOOL flag;

    GetWallTime timer;
    timer.Start();

    switch (method) {
        case OCPNLMethod::IMPEC:
            flag = impec.UpdateProperty(rs, ctrl);
            break;
        case OCPNLMethod::FIM:
            flag = fim.UpdateProperty(rs, ctrl);
            break;
        case OCPNLMethod::AIMc:
            flag = aimc.UpdateProperty(rs, ctrl);
            break;
        case OCPNLMethod::FIMddm:
            flag = fim_ddm.IsoT_FIM::UpdateProperty(rs, ctrl);
            break;
        default:
            OCP_ABORT("Wrong method type!");
    }

    OCPTIME_UPDATE_GRID += timer.Stop() / TIME_S2MS;

    return flag;
}


/// Finish up Newton-Raphson iteration for IMPEC and FIM.
OCP_BOOL IsothermalSolver::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    switch (method) {
        case OCPNLMethod::IMPEC:
            return impec.FinishNR(rs);
        case OCPNLMethod::FIM:
            return fim.FinishNR(rs, ctrl);
        case OCPNLMethod::AIMc:
            return aimc.FinishNR(rs, ctrl);
        case OCPNLMethod::FIMddm:
            return fim_ddm.FinishNR(rs, ctrl);
        default:
            OCP_ABORT("Wrong method type!");
    }
}


/// Finish up time step for IMPEC and FIM.
void IsothermalSolver::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    switch (method) {
        case OCPNLMethod::IMPEC:
            impec.FinishStep(rs, ctrl);
            break;
        case OCPNLMethod::FIM:
            fim.FinishStep(rs, ctrl);
            break;
        case OCPNLMethod::AIMc:
            aimc.FinishStep(rs, ctrl);
            break;
        case OCPNLMethod::FIMddm:
            fim_ddm.IsoT_FIM::FinishStep(rs, ctrl);
            break;
        default:
            OCP_ABORT("Wrong method type!");
    }
}


const OCPNRsuite& IsothermalSolver::GetNRsuite() const
{
    switch (method) {
    case OCPNLMethod::IMPEC:
        return impec.GetNRsuite();
        break;
    case OCPNLMethod::FIM:
        return fim.GetNRsuite();
        break;
    case OCPNLMethod::AIMc:
        return aimc.GetNRsuite();
        break;
    case OCPNLMethod::FIMddm:
        return fim_ddm.GetNRsuite();
        break;
    default:
        OCP_ABORT("Wrong method type!");
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/21/2021      Create file                          */
/*  Chensong Zhang      Jan/08/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/