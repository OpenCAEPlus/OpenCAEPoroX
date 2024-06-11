/*! \file    IsoThermalSolver.cpp
 *  \brief   IsoThermalSolver class definition
 *  \author  Shizhe Li
 *  \date    Oct/21/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "IsoThermalSolver.hpp"

/// Setup solution methods used including solver and preconditioner
void IsothermalSolver::SetupMethod(Reservoir& rs, const OCPControl& ctrl)
{
    const auto& methods = ctrl.SM.GetMethod();

    for (USI i = 0; i < methods.size(); i++) {
        switch (methods[i]) {
        case OCPNLMethod::IMPEC:        
            impec.Setup(rs, ctrl);
            impec.SetWorkLS(LSolver.Setup(ctrl.SM.GetModel(), ctrl.SM.GetWorkDir(), ctrl.SM.GetLsFile(i), rs.GetDomain(), 1), i);
            break;
        case OCPNLMethod::FIM:
            fim.Setup(rs, ctrl);
            fim.SetWorkLS(LSolver.Setup(ctrl.SM.GetModel(), ctrl.SM.GetWorkDir(), ctrl.SM.GetLsFile(i), rs.GetDomain(), rs.GetComNum() + 1), i);
            break;
        case OCPNLMethod::AIMc:
            aimc.Setup(rs, ctrl);
            aimc.SetWorkLS(LSolver.Setup(ctrl.SM.GetModel(), ctrl.SM.GetWorkDir(), ctrl.SM.GetLsFile(i), rs.GetDomain(), rs.GetComNum() + 1), i);
            break;
        case OCPNLMethod::FIMddm:
            fim_ddm.Setup(rs, ctrl);
            fim_ddm.SetWorkLS(LSolver.Setup(ctrl.SM.GetModel(), ctrl.SM.GetWorkDir(), ctrl.SM.GetLsFile(i), rs.GetDomain(), rs.GetComNum() + 1), i);
            break;
        default:
            OCP_ABORT("Wrong method type!");
        }
    }

    mainMethod = methods[0];
    preMethod  = ctrl.SM.InitMethod();
    curMethod  = preMethod;
}


/// Setup solution methods, including IMPEC and FIM.
void IsothermalSolver::InitReservoir(Reservoir& rs)
{
    switch (curMethod) {
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
            fim_ddm.InitReservoir(rs);
            break;
        default:
            OCP_ABORT("Wrong method type!");
    }
}


const OCPNRsuite& IsothermalSolver::GoOneStep(Reservoir& rs, OCPControl& ctrl)
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
        if (!SolveLinearSystem(rs, ctrl)) {
            continue;
        }
        if (!UpdateProperty(rs, ctrl)) {
            continue;
        }
        if (FinishNR(rs, ctrl)) break;
    }

    // Finish current time step
    FinishStep(rs, ctrl);

    return GetNRsuite();
}


/// Prepare solution methods, including IMPEC and FIM.
void IsothermalSolver::Prepare(Reservoir& rs, OCPControl& ctrl)
{
    curMethod = preMethod;
	switch (curMethod) {
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
		fim_ddm.Prepare(rs, ctrl.time.GetCurrentDt());
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

    switch (curMethod) {
        case OCPNLMethod::IMPEC:
            LSolver.SetWorkLS(impec.GetWorkLS());
            impec.AssembleMat(LSolver, rs, dt);
            break;
        case OCPNLMethod::FIM:
            LSolver.SetWorkLS(fim.GetWorkLS());
            fim.AssembleMat(LSolver, rs, dt);
            break;
        case OCPNLMethod::AIMc:
            LSolver.SetWorkLS(aimc.GetWorkLS());
            aimc.AssembleMat(LSolver, rs, dt);
            break;
        case OCPNLMethod::FIMddm:
            LSolver.SetWorkLS(fim_ddm.GetWorkLS());
            fim_ddm.AssembleMat(LSolver, rs, dt);
            break;
        default:
            OCP_ABORT("Wrong method type!");
    }
   
    OCPTIME_ASSEMBLE_MAT += timer.Stop();
}


/// Solve linear systems for IMPEC and FIM.
OCP_BOOL IsothermalSolver::SolveLinearSystem(Reservoir& rs, OCPControl& ctrl)
{ 
    switch (curMethod) {
        case OCPNLMethod::IMPEC:
            return impec.SolveLinearSystem(LSolver, rs, ctrl);
            break;
        case OCPNLMethod::FIM:
            return fim.SolveLinearSystem(LSolver, rs, ctrl);
            break;
        case OCPNLMethod::AIMc:
            return aimc.SolveLinearSystem(LSolver, rs, ctrl);
            break;
        case OCPNLMethod::FIMddm:
            return fim_ddm.SolveLinearSystem(LSolver, rs, ctrl);
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

    switch (curMethod) {
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
            flag = fim_ddm.UpdateProperty(rs, ctrl);
            break;
        default:
            OCP_ABORT("Wrong method type!");
    }

    OCPTIME_UPDATE_GRID += timer.Stop();

    return flag;
}


/// Finish up Newton-Raphson iteration for IMPEC and FIM.
OCP_BOOL IsothermalSolver::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    OCP_BOOL conFlag = OCP_FALSE;

    switch (curMethod) {
        case OCPNLMethod::IMPEC:
            conFlag = impec.FinishNR(rs);
            break;
        case OCPNLMethod::FIM:
            conFlag = fim.FinishNR(rs, ctrl);
            break;
        case OCPNLMethod::AIMc:
            conFlag = aimc.FinishNR(rs, ctrl);
            break;
        case OCPNLMethod::FIMddm:
            conFlag = fim_ddm.FinishNR(rs, ctrl);
            break;
        default:
            OCP_ABORT("Wrong method type!");
    }

    if (conFlag) {
        if (curMethod == mainMethod) {
            return OCP_TRUE;
        }
        curMethod = mainMethod;
        if (curMethod == OCPNLMethod::FIM) {
            if (CURRENT_RANK == 0)
                cout << "FIMddm iters = " << fim_ddm.NR.GetIterNR() << "  " << endl;
            if (fim.TransferToFIM(fim_ddm.global_res0, rs, ctrl) == OCP_TRUE) {
                return OCP_FALSE;
            }
            else {
                return OCP_TRUE;
            }
        }
        return OCP_TRUE;
    }
    else {
        return OCP_FALSE;
    }
}


/// Finish up time step for IMPEC and FIM.
void IsothermalSolver::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    if (ctrl.SM.GetMethod().size() == 1) {
        switch (curMethod) {
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
            fim_ddm.FinishStep(rs, ctrl);
            break;
        default:
            OCP_ABORT("Wrong method type!");
        }
    }
    else {
        FinishStepComb(rs, ctrl);
    }
}

void IsothermalSolver::FinishStepComb(Reservoir& rs, OCPControl& ctrl)
{
    if (mainMethod == OCPNLMethod::FIM && preMethod == OCPNLMethod::FIMddm) {

        rs.CalIPRT(ctrl.time.GetCurrentDt());
        fim.NR.CalMaxChangeTime(rs);
        ctrl.CalNextTimeStep(fim.NR, { "dP", "dS", "iter" });

        fim_ddm.SetStarBulkSet(rs.GetBulk(), rs.GetDomain());

        OCPITER_NR_DDM  += fim_ddm.NR.GetIterNR();
        OCPITER_NRW_DDM += fim_ddm.NR.GetIterNRw();
        OCPITER_LS_DDM  += fim_ddm.NR.GetIterLS();
    }
}


const OCPNRsuite& IsothermalSolver::GetNRsuite() const
{
    switch (curMethod) {
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