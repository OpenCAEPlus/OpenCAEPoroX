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
      
    OutputTimeMain(cout.rdbuf());
    OutputTimeProcess();
}


void OpenCAEPoroX::OutputTimeMain(streambuf* mysb) const
{
    if (CURRENT_RANK == MASTER_PROCESS) {

        streambuf* oldcout = cout.rdbuf(mysb);
        // find an appropriate size for printing times
        int fixWidth = OCP_MAX(log10(control.current_time), log10(OCP_MAX(OCPTIME_TOTAL, 1.0))) + 6;
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
        cout << "Simulation time:                " << setw(fixWidth) << OCPTIME_TOTAL
            << " (Seconds)" << endl;
        cout << " - % Input & Partition ........." << setw(fixWidth)
            << 100.0 * OCPTIME_PARTITION / OCPTIME_TOTAL << " (" << OCPTIME_PARTITION
            << "s)" << endl;
        cout << " - % Input Reservoir ..........." << setw(fixWidth)
            << 100.0 * OCPTIME_READPARAM / OCPTIME_TOTAL << " (" << OCPTIME_READPARAM
            << "s)" << endl;
        cout << " - % Setup Simulator ..........." << setw(fixWidth)
            << 100.0 * OCPTIME_SETUP_SIM / OCPTIME_TOTAL << " (" << OCPTIME_SETUP_SIM
            << "s)" << endl;
        cout << " - % Initialization ............" << setw(fixWidth)
            << 100.0 * OCPTIME_INIT_RESERVOIR / OCPTIME_TOTAL << " (" << OCPTIME_INIT_RESERVOIR
            << "s)" << endl;
        cout << " - % Assembling ................" << setw(fixWidth)
            << 100.0 * OCPTIME_ASSEMBLE_MAT / OCPTIME_TOTAL << " ("
            << OCPTIME_ASSEMBLE_MAT << "s)" << endl;
        cout << " - % Assembling for LS ........." << setw(fixWidth)
            << 100.0 * OCPTIME_ASSEMBLE_MAT_FOR_LS / OCPTIME_TOTAL << " ("
            << OCPTIME_ASSEMBLE_MAT_FOR_LS << "s)" << endl;
        cout << " - % Linear Solver ............." << setw(fixWidth)
            << 100.0 * OCPTIME_LSOLVER / OCPTIME_TOTAL << " ("
            << OCPTIME_LSOLVER << "s)" << endl;
        cout << " - % Newton Step ..............." << setw(fixWidth)
            << 100.0 * OCPTIME_NRSTEP / OCPTIME_TOTAL << " ("
            << OCPTIME_NRSTEP << "s)" << endl;
        cout << " - % Updating Properties ......." << setw(fixWidth)
            << 100.0 * OCPTIME_UPDATE_GRID / OCPTIME_TOTAL << " ("
            << OCPTIME_UPDATE_GRID << "s)" << endl;
        cout << " - % Output ...................." << setw(fixWidth)
            << 100.0 * OCPTIME_OUTPUT / OCPTIME_TOTAL << " ("
            << OCPTIME_OUTPUT << "s)" << endl;
        cout << " - % Communication(collect) ...." << setw(fixWidth)
            << 100.0 * OCPTIME_COMM_COLLECTIVE / OCPTIME_TOTAL << " ("
            << OCPTIME_COMM_COLLECTIVE << "s)" << endl;
        cout << " - % Communication(P2P) ........" << setw(fixWidth)
            << 100.0 * OCPTIME_COMM_P2P / OCPTIME_TOTAL << " ("
            << OCPTIME_COMM_P2P << "s)" << endl;

        cout << "==================================================" << endl;

        cout.rdbuf(oldcout);
    }
}


void OpenCAEPoroX::OutputTimeProcess() const
{
    // Record information of each process
    const Domain& domain = reservoir.GetDomain();
    if (domain.numproc > 1) {
        const OCP_INT record_var_num = 4;
        vector<OCP_DBL> record_local{ OCPTIME_UPDATE_GRID , OCPTIME_ASSEMBLE_MAT,
                OCPTIME_COMM_P2P, static_cast<OCP_DBL>(reservoir.GetInteriorBulkNum()) };
        vector<OCP_DBL> record_total;

        if (CURRENT_RANK == MASTER_PROCESS) {
            record_total.resize(record_var_num * domain.numproc);
            MPI_Gather(record_local.data(), record_var_num, MPI_DOUBLE, record_total.data(), record_var_num, MPI_DOUBLE, MASTER_PROCESS, domain.myComm);

            // Calculate averge num           
            OCP_DBL       aveTimeUG  = 0;
            OCP_DBL       aveTimeAM  = 0;
            OCP_DBL       aveTimeP2P = 0;
            const OCP_DBL aveGrid    = 1.0 * domain.GetNumGridTotal() / domain.numproc;
            OCP_DBL       minTimeUG{ record_local[0] }, maxTimeUG{ record_local[0] };
            OCP_DBL       minTimeAM{ record_local[1] }, maxTimeAM{ record_local[1] };
            OCP_DBL       minTimeP2P{ record_local[2] }, maxTimeP2P{ record_local[2] };
            OCP_USI       minNG{ record_local[3] }, maxNG{ record_local[3] };

            for (OCP_USI p = 0; p < domain.numproc; p++) {
                OCP_DBL tmpUG  = record_total[p * record_var_num + 0];
                OCP_DBL tmpAM  = record_total[p * record_var_num + 1];
                OCP_DBL tmpP2P = record_total[p * record_var_num + 2];
                OCP_DBL tmpNG  = record_total[p * record_var_num + 3];

                aveTimeUG  += tmpUG;
                aveTimeAM  += tmpAM;
                aveTimeP2P += tmpP2P;

                minTimeUG  = minTimeUG < tmpUG ? minTimeUG : tmpUG;
                maxTimeUG  = maxTimeUG > tmpUG ? maxTimeUG : tmpUG;
                minTimeAM  = minTimeAM < tmpAM ? minTimeAM : tmpAM;
                maxTimeAM  = maxTimeAM > tmpAM ? maxTimeAM : tmpAM;
                minTimeP2P = minTimeP2P < tmpP2P ? minTimeP2P : tmpP2P;
                maxTimeP2P = maxTimeP2P > tmpP2P ? maxTimeP2P : tmpP2P;
                minNG      = minNG < tmpNG ? minNG : tmpNG;
                maxNG      = maxNG > tmpNG ? maxNG : tmpNG;
            }

            aveTimeUG  /= domain.numproc;
            aveTimeAM  /= domain.numproc;
            aveTimeP2P /= domain.numproc;

            // Calculate standard variance
            OCP_DBL varTimeUG  = 0;
            OCP_DBL varTimeAM  = 0;
            OCP_DBL varTimeP2P = 0;
            OCP_DBL varGrid    = 0;
            for (OCP_USI p = 0; p < domain.numproc; p++) {
                varTimeUG  += pow((record_total[p * record_var_num + 0] - aveTimeUG), 2);
                varTimeAM  += pow((record_total[p * record_var_num + 1] - aveTimeAM), 2);
                varTimeP2P += pow((record_total[p * record_var_num + 2] - aveTimeP2P), 2);
                varGrid    += pow((record_total[p * record_var_num + 3] - aveGrid), 2);
            }
            varTimeUG  = sqrt(varTimeUG) / domain.numproc;
            varTimeAM  = sqrt(varTimeAM) / domain.numproc;
            varTimeP2P = sqrt(varTimeP2P) / domain.numproc;
            varGrid    = sqrt(varGrid) / domain.numproc;

            // output general information to screen
            cout << fixed << setprecision(3);
            cout << "Item                 " << setw(12) << " Average Time " << setw(12) << "Varance" << setw(12) << "Max" << setw(12) << "Min" << " \n";
            cout << "Updating Properties  " << setw(12) << aveTimeUG << "s" << setw(12) << varTimeUG << "s" << setw(12) << maxTimeUG << "s" << setw(12) << minTimeUG << "s\n";
            cout << "Assembling           " << setw(12) << aveTimeAM << "s" << setw(12) << varTimeAM << "s" << setw(12) << maxTimeAM << "s" << setw(12) << minTimeAM << "s\n";
            cout << "Communication(P2P)   " << setw(12) << aveTimeP2P << "s" << setw(12) << varTimeP2P << "s" << setw(12) << maxTimeP2P << "s" << setw(12) << minTimeP2P << "s\n";
            cout << "Grid Num             " << setw(12) << aveGrid << " " << setw(12) << varGrid << " " << setw(12) << maxNG << " " << setw(12) << minNG << " \n";
            cout << "==================================================" << endl;
            // output detailed inforamtion to files
            if (true) {
                ofstream myFile;
                myFile.open(control.workDir + "statistics.out");

                ios::sync_with_stdio(false);
                myFile.tie(0);

                myFile << fixed << setprecision(3);
                myFile << setw(6) << "Rank"
                    << setw(30) << "Updating Properties (s)"
                    << setw(30) << "Assembling (s)"
                    << setw(30) << "Communication(P2P) (s)"
                    << setw(30) << "Grid Num" << "\n";
                for (OCP_USI p = 0; p < domain.numproc; p++) {
                    myFile << setw(6) << p
                        << setprecision(3) << setw(30) << record_total[p * record_var_num + 0]
                        << setprecision(3) << setw(30) << record_total[p * record_var_num + 1]
                        << setprecision(3) << setw(30) << record_total[p * record_var_num + 2]
                        << setprecision(0) << setw(30) << record_total[p * record_var_num + 3] << "\n";
                }
                myFile << fixed << setprecision(3);
                myFile << "\n==================================================\n";
                myFile << "Item                 " << setw(12) << " Average Time " << setw(12) << "Varance" << setw(12) << "Max" << setw(12) << "Min" << " \n";
                myFile << "Updating Properties  " << setw(12) << aveTimeUG << "s" << setw(12) << varTimeUG << "s" << setw(12) << maxTimeUG << "s" << setw(12) << minTimeUG << "s\n";
                myFile << "Assembling           " << setw(12) << aveTimeAM << "s" << setw(12) << varTimeAM << "s" << setw(12) << maxTimeAM << "s" << setw(12) << minTimeAM << "s\n";
                myFile << "Communication(P2P)   " << setw(12) << aveTimeP2P << "s" << setw(12) << varTimeP2P << "s" << setw(12) << maxTimeP2P << "s" << setw(12) << minTimeP2P << "s\n";
                myFile << "Grid Num             " << setw(12) << aveGrid << " " << setw(12) << varGrid << " " << setw(12) << maxNG << " " << setw(12) << minNG<< " \n";

                OutputTimeMain(myFile.rdbuf());

                myFile.close();
            }
        }
        else {
            MPI_Gather(record_local.data(), record_var_num, MPI_DOUBLE, record_total.data(), record_var_num, MPI_DOUBLE, MASTER_PROCESS, domain.myComm);
        }
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Dec/05/2021      Format file                          */
/*----------------------------------------------------------------------------*/