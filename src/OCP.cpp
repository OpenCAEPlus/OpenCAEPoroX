/*! \file    OCP.cpp
 *  \brief   OpenCAEPoroX class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCP.hpp"

/// Read Param from input file
void OpenCAEPoroX::InputDistParam(const string& filename, PreProcess& prepro, const OCP_INT& myRank)
{

    GetWallTime timer;
    timer.Start();

    OCP_BOOL disable_grid = OCP_FALSE;
    if (myRank != MASTER_PROCESS)
        disable_grid = OCP_TRUE;

    ParamRead rp(disable_grid);
    rp.ReadInputFile(filename);

    reservoir.InputParam(prepro, rp);
    control.InputParam(rp.paramControl);
    output.InputParam(rp.paramOutput);

    OCPTIME_READPARAM = timer.Stop() / 1000;
    OCPTIME_TOTAL     += OCPTIME_READPARAM;
}


/// Call setup procedures for reservoir, output, and linear solver.
void OpenCAEPoroX::SetupSimulator(const USI& argc, const char* options[])
{
    const Domain& domain = reservoir.GetDomain();

    GetWallTime timer;
    timer.Start();

    control.SetupFastControl(argc, options); // Read Fast control

    if (CURRENT_RANK == MASTER_PROCESS) {
        switch (control.GetModel()) {
        case ISOTHERMALMODEL:
            if (control.printLevel >= PRINT_MIN) {
                cout << endl << "Dynamic simulation for isothermal models" << endl;
            }
            break;
        case THERMALMODEL:
            if (control.printLevel >= PRINT_MIN) {
                cout << endl << "Dynamic simulation for thermal models" << endl;
            }
            break;
        default:
            OCP_ABORT("Wrong model type specified!");
        }
    }

    solver.Setup(reservoir, control); // Setup static info for solver

    output.Setup(reservoir, control, domain); // Setup output for dynamic simulation

    control.Setup(domain);

    double finalTime = timer.Stop() / 1000;
    if (control.printLevel >= PRINT_MIN && CURRENT_RANK == MASTER_PROCESS) {
        cout << endl
             << "Setup simulation done. Wall time : " << fixed << setprecision(3)
             << finalTime << " Sec" << endl;
    }

    OCPTIME_SETUP_SIM = finalTime;
    OCPTIME_TOTAL     += OCPTIME_SETUP_SIM;
}


/// Initialize the reservoir class.
void OpenCAEPoroX::InitReservoir()
{
    GetWallTime timer;
    timer.Start();

    solver.InitReservoir(reservoir);

    double finalTime = timer.Stop() / 1000;
    if (control.printLevel >= PRINT_MIN && CURRENT_RANK == MASTER_PROCESS) {
        cout << endl
             << "Initialization done. Wall time : " << fixed << setprecision(3)
             << finalTime << " Sec" << endl;
    }

    OCPTIME_INIT_RESERVOIR = finalTime;
    OCPTIME_TOTAL          += OCPTIME_INIT_RESERVOIR;
}

// Call IMPEC, FIM, AIM, etc for dynamic simulation.
void OpenCAEPoroX::RunSimulation()
{
    if (CURRENT_RANK == MASTER_PROCESS) {
        switch (control.GetMethod()) {
        case IMPEC:
            if (control.printLevel >= PRINT_MIN) {
                cout << "\nDynamic simulation with IMPEC\n" << endl;
            }
            break;
        case FIM:
            if (control.printLevel >= PRINT_MIN) {
                cout << "\nDynamic simulation with FIM\n" << endl;
            }
            break;
        case FIMn:
            if (control.printLevel >= PRINT_MIN) {
                cout << "\nDynamic simulation with FIMn\n" << endl;
            }
            break;
        case AIMc:
            if (control.printLevel >= PRINT_MIN) {
                cout << "\nDynamic simulation with AIMc\n" << endl;
            }
            break;
        default:
            OCP_ABORT("Wrong method type is used!");
        }
    }

    solver.RunSimulation(reservoir, control, output);
}

/// Print summary information on screen and SUMMARY.out file.
void OpenCAEPoroX::OutputResults() const
{
    GetWallTime timer;
    timer.Start();
    output.PrintInfo();
    output.PostProcess();
    OCPTIME_TOTAL += timer.Stop() / 1000;
    // find an appropriate size for printing times
    int fixWidth = OCP_MAX(log10(control.current_time), log10(OCP_MAX(OCPTIME_TOTAL, 1.0))) + 6;
    if (CURRENT_RANK == MASTER_PROCESS) {
        
        cout << "==================================================" << endl;

        // print numbers of steps
        cout << "Final time:                  " << right << fixed << setprecision(3)
            << setw(fixWidth) << control.current_time << " (Days)" << endl;
        cout << " - Avg time step size ......." << setw(fixWidth)
            << control.current_time / control.numTstep << " (" << control.numTstep
            << " steps)" << endl;
        cout << " - Avg Newton steps ........." << setw(fixWidth)
            << static_cast<double>(control.iterNR_total) / control.numTstep << " ("
            << control.iterNR_total << " succeeded + " << control.wastedIterNR
            << " wasted)" << endl;
        cout << " - Avg linear steps ........." << setw(fixWidth)
            << static_cast<double>(control.iterLS_total) / control.iterNR_total << " ("
            << control.iterLS_total << " succeeded + " << control.wastedIterLS
            << " wasted)" << endl;

        // print time usages
        cout << "Simulation time:             " << setw(fixWidth) << OCPTIME_TOTAL
            << " (Seconds)" << endl;
        cout << " - % Input & Partition ......" << setw(fixWidth)
            << 100.0 * OCPTIME_PARTITION / OCPTIME_TOTAL << " (" << OCPTIME_PARTITION
            << "s)" << endl;
        cout << " - % Input Reservoir ........" << setw(fixWidth)
            << 100.0 * OCPTIME_READPARAM / OCPTIME_TOTAL << " (" << OCPTIME_READPARAM
            << "s)" << endl;
        cout << " - % Setup Simulator ........" << setw(fixWidth)
            << 100.0 * OCPTIME_SETUP_SIM / OCPTIME_TOTAL << " (" << OCPTIME_SETUP_SIM
            << "s)" << endl;
        cout << " - % Initialization ........." << setw(fixWidth)
            << 100.0 * OCPTIME_INIT_RESERVOIR / OCPTIME_TOTAL << " (" << OCPTIME_INIT_RESERVOIR
            << "s)" << endl;
        cout << " - % Assembling ............." << setw(fixWidth)
            << 100.0 * OCPTIME_ASSEMBLE_MAT / OCPTIME_TOTAL << " ("
            << OCPTIME_ASSEMBLE_MAT << "s)" << endl;
        cout << " - % Assembling for LS ......" << setw(fixWidth)
            << 100.0 * OCPTIME_ASSEMBLE_MAT_FOR_LS / OCPTIME_TOTAL << " ("
            << OCPTIME_ASSEMBLE_MAT_FOR_LS << "s)" << endl;
        cout << " - % Linear Solver .........." << setw(fixWidth)
            << 100.0 * OCPTIME_LSOLVER / OCPTIME_TOTAL << " ("
            << OCPTIME_LSOLVER << "s)" << endl;
        cout << " - % Newton Step ............" << setw(fixWidth)
            << 100.0 * OCPTIME_NRSTEP / OCPTIME_TOTAL << " ("
            << OCPTIME_NRSTEP << "s)" << endl;
        cout << " - % Updating Properties ...." << setw(fixWidth)
            << 100.0 * OCPTIME_UPDATE_GRID / OCPTIME_TOTAL << " ("
            << OCPTIME_UPDATE_GRID << "s)" << endl;
        cout << " - % Output ................." << setw(fixWidth)
            << 100.0 * OCPTIME_OUTPUT / OCPTIME_TOTAL << " ("
            << OCPTIME_OUTPUT << "s)" << endl;

        cout << "==================================================" << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    cout << fixed << setprecision(3)
        << "Rank " << setw(4) << CURRENT_RANK 
        << ": Updating properties : " << OCPTIME_UPDATE_GRID << "s   "
        << "Assembling : " << OCPTIME_ASSEMBLE_MAT << "s   " 
        << "NumGrid =  " << reservoir.GetInteriorBulkNum() << endl;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Dec/05/2021      Format file                          */
/*----------------------------------------------------------------------------*/