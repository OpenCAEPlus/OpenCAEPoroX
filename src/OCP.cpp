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

    ParamRead rp;
    rp.ReadInputFile(filename);

    reservoir.InputParam(prepro, rp);
    control.InputParam(rp.paramControl);
    output.InputParam(rp.paramOutput);

    OCPTIME_READPARAM = timer.Stop() / TIME_S2MS;
    OCPTIME_TOTAL     += OCPTIME_READPARAM;
}


/// Call setup procedures for reservoir, output, and linear solver.
void OpenCAEPoroX::SetupSimulator(const USI& argc, const char* options[])
{
    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Setup Simulator -- begin");
    }


    const Domain& domain = reservoir.GetDomain();

    GetWallTime timer;
    timer.Start();

    control.SetupFastControl(argc, options); // Read Fast control

    if (CURRENT_RANK == MASTER_PROCESS) {
        switch (control.GetModel()) {
        case OCPModel::isothermal:
            if (control.printLevel >= PRINT_MIN) {
                cout << endl << "Dynamic simulation for isothermal models" << endl;
            }
            break;
        case OCPModel::thermal:
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

    double finalTime = timer.Stop() / TIME_S2MS;
    if (control.printLevel >= PRINT_MIN && CURRENT_RANK == MASTER_PROCESS) {
        cout << endl
             << "Setup simulation done. Wall time : " << fixed << setprecision(3)
             << finalTime << " Sec" << endl;
    }

    OCPTIME_SETUP_SIM = finalTime;
    OCPTIME_TOTAL     += OCPTIME_SETUP_SIM;


    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Setup Simulator -- end");
    }
}


/// Initialize the reservoir class.
void OpenCAEPoroX::InitReservoir()
{
    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Initialize Resevoir -- begin");
    }

    GetWallTime timer;
    timer.Start();

    solver.InitReservoir(reservoir);

    double finalTime = timer.Stop() / TIME_S2MS;
    if (control.printLevel >= PRINT_MIN && CURRENT_RANK == MASTER_PROCESS) {
        std::ostringstream initTime;
        initTime << fixed << setprecision(3) << finalTime;
        OCP_INFO("Initialize Resevoir -- end. Wall time : " + to_string(finalTime) + " Sec");
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
    OCPTIME_TOTAL += timer.Stop() / TIME_S2MS;
      
    OutputTimeMain(cout.rdbuf());
    OutputTimeProcess();
}


void OpenCAEPoroX::OutputTimeMain(streambuf* mysb) const
{
    if (CURRENT_RANK == MASTER_PROCESS) {

        streambuf* oldcout = cout.rdbuf(mysb);
        // find an appropriate size for printing times
        int fixWidth = OCP_MAX(log10(control.time.GetCurrentTime()), log10(OCP_MAX(OCPTIME_TOTAL, 1.0))) + 6;
        cout << "==================================================" << endl;

        // print numbers of steps
        cout << "Final time:                  " << right << fixed << setprecision(3)
            << setw(fixWidth) << control.time.GetCurrentTime() << "(" + TIMEUNIT + ")" << endl;
        cout << " - Avg time step size ......." << setw(fixWidth)
            << control.time.GetCurrentTime() / output.iters.GetNumTimeStep() << " (" << output.iters.GetNumTimeStep()
            << " steps)" << endl;
        cout << " - Avg Newton steps ........." << setw(fixWidth)
            << static_cast<double>(output.iters.GetNRt()) / output.iters.GetNumTimeStep() << " ("
            << output.iters.GetNRt() << " succeeded + " << output.iters.GetNRwt()
            << " wasted)" << endl;
        cout << " - Avg linear steps ........." << setw(fixWidth)
            << static_cast<double>(output.iters.GetLSt()) / output.iters.GetNRt() << " ("
            << output.iters.GetLSt() << " succeeded + " << output.iters.GetLSwt()
            << " wasted)" << endl;

        // print time usages
        cout << "Simulation time:                " << setw(fixWidth) << OCPTIME_TOTAL
            << " (Seconds)" << endl;
        cout << " - % Input & Partition ........." << setw(fixWidth)
            << 100.0 * OCPTIME_PARTITION / OCPTIME_TOTAL << " (" << OCPTIME_PARTITION
            << "s)" << endl;
        cout << " - % Partition - ParMetis ......" << setw(fixWidth)
            << 100.0 * OCPTIME_PARMETIS / OCPTIME_TOTAL << " (" << OCPTIME_PARMETIS
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
        cout << " - % Converting for LS interface" << setw(fixWidth)
            << 100.0 * OCPTIME_CONVERT_MAT_FOR_LS_IF / OCPTIME_TOTAL << " ("
            << OCPTIME_CONVERT_MAT_FOR_LS_IF << "s)" << endl;
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
        cout << "==================================================" << endl;

        cout.rdbuf(oldcout);
    }
}


void OpenCAEPoroX::OutputTimeProcess() const
{
    // Record information of each process
    const Domain& domain = reservoir.GetDomain();
    if (domain.numproc > 1) {
     
        const vector<OCP_DBL> record_local{ 
            static_cast<OCP_DBL>(reservoir.GetInteriorBulkNum()),
            OCPTIME_UPDATE_GRID,
            OCPTIME_ASSEMBLE_MAT,
            OCPTIME_COMM_COLLECTIVE,
            OCPTIME_COMM_P2P,
            OCPTIME_COMM_1ALLREDUCE,
            OCPTIME_NRSTEP,
            OCPTIME_NRSTEPC};
        vector<OCP_DBL> record_total;
        const OCP_INT record_var_num = record_local.size();

        if (CURRENT_RANK == MASTER_PROCESS) {

            class statisticsVar
            {
            public:
                statisticsVar(const string& name, const USI& n, const USI& pre) {
                    itemName = name;
                    val.resize(n, 0);
                    precision = pre;
                    clens = itemName.size() + 5;
                }
                string          itemName;  ///< name of item            
                vector<OCP_DBL> val;       ///< averge, varance, max, min
                USI             precision; ///< precision for output
                USI             clens;     ///< length of itemname + 5
            };
            const vector<string> myItems{
                "Item", "Average","Varance","Max","Min"
            };
            const USI len = myItems.size() - 1;

            vector<statisticsVar> staVar{
                statisticsVar("Grid Num", len, 0),
                statisticsVar("Updating Properties(s)", len, 3),
                statisticsVar("Assembling(s)", len, 3),
                statisticsVar("Communication(collect)(s)", len, 3),
                statisticsVar("Communication(P2P)(s)", len, 3),                
                statisticsVar("1AllReduce(OCPCheck)(s)", len, 3),
                statisticsVar("Newton Step(s)", len, 3),
                statisticsVar("Newton Step(c)(s)", len, 3),
            };

            OCP_ASSERT(record_var_num == staVar.size(), "wrong staVar");

            record_total.resize(record_var_num * domain.numproc);
            MPI_Gather(record_local.data(), record_var_num, OCPMPI_DBL, record_total.data(), record_var_num, OCPMPI_DBL, MASTER_PROCESS, domain.myComm);

            // Calculate average, max, min
            for (USI n = 0; n < record_var_num; n++) {
                staVar[n].val[2] = record_local[n];
                staVar[n].val[3] = record_local[n];
            }

            for (OCP_USI p = 0; p < domain.numproc; p++) {
                for (USI n = 0; n < record_var_num; n++) {
                    const OCP_DBL tmp = record_total[p * record_var_num + n];
                    staVar[n].val[0] += tmp;
                    staVar[n].val[2] = staVar[n].val[2] < tmp ? tmp : staVar[n].val[2];
                    staVar[n].val[3] = staVar[n].val[3] > tmp ? tmp : staVar[n].val[3];
                }
            }

            for (USI n = 0; n < record_var_num; n++) {
                staVar[n].val[0] /= domain.numproc;
            }

            // Calculate standard variance
            for (OCP_USI p = 0; p < domain.numproc; p++) {
                for (USI n = 0; n < record_var_num; n++) {
                    staVar[n].val[1] += pow((record_total[p * record_var_num + n] - staVar[n].val[0]), static_cast<OCP_DBL>(2));
                }
            }
            for (USI n = 0; n < record_var_num; n++) {
                staVar[n].val[1] = pow(staVar[n].val[1] / domain.numproc, static_cast<OCP_DBL>(0.5));
            }

            // output general information to screen
            cout << std::left << setw(25) << myItems[0];
            for (USI n = 1; n < myItems.size(); n++) {
                cout << std::right << setw(12) << myItems[n];
            }
            cout << endl;
            for (USI n = 0; n < record_var_num; n++) {
                cout << std::left << setw(25) << staVar[n].itemName;
                for (USI i = 0; i < len; i++)
                    cout << std::right << setw(12) << setprecision(staVar[n].precision) << staVar[n].val[i];
                cout << endl;
            }
            cout << "==================================================" << endl;
            // output detailed inforamtion to files
            if (true) {
                ofstream myFile;
                myFile.open(control.GetWorkDir() + "statistics.out");

                ios::sync_with_stdio(false);
                myFile.tie(0);

                myFile << fixed << setprecision(3);
                myFile << setw(6) << "Rank";
                for (USI i = 0; i < record_var_num; i++) {
                    myFile << setw(staVar[i].clens) << staVar[i].itemName;
                }
                myFile << "\n";                   
                for (OCP_USI p = 0; p < domain.numproc; p++) {
                    myFile << setw(6) << p;
                    for (USI i = 0; i < record_var_num; i++) {
                        myFile << setprecision(staVar[i].precision) << setw(staVar[i].clens) << record_total[p * record_var_num + i];
                    }
                    myFile << "\n";
                }
                myFile << "\n==================================================\n";
                myFile << std::left << setw(25) << myItems[0];
                for (USI n = 1; n < myItems.size(); n++) {
                    myFile << std::right << setw(12) << myItems[n];
                }
                myFile << endl;
                for (USI n = 0; n < record_var_num; n++) {
                    myFile << std::left << setw(25) << staVar[n].itemName;
                    for (USI i = 0; i < len; i++)
                        myFile << std::right << setw(12) << setprecision(staVar[n].precision) << staVar[n].val[i];
                    myFile << endl;
                }

                OutputTimeMain(myFile.rdbuf());

                myFile.close();
            }
        }
        else {
            MPI_Gather(record_local.data(), record_var_num, OCPMPI_DBL, record_total.data(), record_var_num, OCPMPI_DBL, MASTER_PROCESS, domain.myComm);
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