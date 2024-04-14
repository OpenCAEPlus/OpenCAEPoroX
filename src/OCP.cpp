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
void OpenCAEPoroX::SetupDistParam(const USI& argc, const char* argv[], PreProcess& prepro, const OCP_INT& myRank)
{

    GetWallTime timer;
    timer.Start();

    ParamRead rp;
    rp.ReadInputFile(argv[1]);

    reservoir.Setup(prepro, rp);
    control.Setup(argc, argv, rp.paramControl, reservoir.GetDomain());
    output.Setup(rp.paramOutput, control, reservoir);

    OCPTIME_READPARAM = timer.Stop();
    OCPTIME_TOTAL     += OCPTIME_READPARAM;
}


/// Call setup procedures for reservoir, output, and linear solver.
void OpenCAEPoroX::SetupSolver()
{
    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Setup Simulator -- begin");
    }

    GetWallTime timer;
    timer.Start();

    solver.Setup(reservoir, control); // Setup static info for solver

    control.OutputModelMethodInfo();

    double finalTime = timer.Stop();
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

    double finalTime = timer.Stop();
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
    control.simTime.Initialize();

    GetWallTime timer;
    timer.Start();
    output.PrintInfoSched(reservoir, control, timer.Stop());

    while (OCP_TRUE) {

        INT d = control.ApplyControl();
        if (d < 0)  break;

        reservoir.ApplyControl(d);   
        control.time.CalInitTimeStep4TSTEP(reservoir.GetWellOptChange());

        /// calculations between TSTEPS
        while (!control.time.IfEndTSTEP()) {
            output.PrintCurrentTimeIter(control);
            const OCPNRsuite& NR = solver.GoOneStep(reservoir, control);
            output.SetValAtTimeStep(reservoir, control, NR);
            if (control.printLevel >= PRINT_ALL) {
                // Print Summary and critical information at every time step
                output.PrintAtTimeStep();
            }

            if (control.simTime.IfStop())  control.StopSim = OCP_TRUE;
            if (control.StopSim)  break;
        }

        output.PrintInfoSched(reservoir, control, timer.Stop());
        // reservoir.allWells.ShowWellStatus(reservoir.bulk);     
        if (control.StopSim) break;
    }
    OCPTIME_TOTAL += timer.Stop();
}

/// Print summary information on screen and SUMMARY.out file.
void OpenCAEPoroX::OutputResults() const
{
    GetWallTime timer;
    timer.Start();
    output.PrintAtTimeStep();
    output.PostProcess();
    OCPTIME_TOTAL += timer.Stop();
      
    OutputTimeMain(cout.rdbuf());
    OutputTimeProcess();
    reservoir.OutInfoFinal();
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
        cout << "CPU time:                       " << setw(fixWidth) << OCPTIME_TOTAL
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
        cout << " - % Wating time ..............." << setw(fixWidth) 
            << 100.0 * control.simTime.GetWaitingTime() / OCPTIME_TOTAL << " (" << control.simTime.GetWaitingTime()
            << "s)" << endl;
        cout << " - % Simulation time:           " << setw(fixWidth)
            << 100.0 * control.simTime.GetTotalSimTime() / OCPTIME_TOTAL << " (" << control.simTime.GetTotalSimTime()
            << "s)" << endl;
        cout << "     - % Assembling ................" << setw(fixWidth)
            << 100.0 * OCPTIME_ASSEMBLE_MAT / control.simTime.GetTotalSimTime() << " ("
            << OCPTIME_ASSEMBLE_MAT << "s)" << endl;
        cout << "     - % Converting for LS interface" << setw(fixWidth)
            << 100.0 * OCPTIME_CONVERT_MAT_FOR_LS_IF / control.simTime.GetTotalSimTime() << " ("
            << OCPTIME_CONVERT_MAT_FOR_LS_IF << "s)" << endl;
        cout << "     - % Linear Solver ............." << setw(fixWidth)
            << 100.0 * OCPTIME_LSOLVER / control.simTime.GetTotalSimTime() << " ("
            << OCPTIME_LSOLVER << "s)" << endl;
        cout << "     - % Newton Step ..............." << setw(fixWidth)
            << 100.0 * OCPTIME_NRSTEP / control.simTime.GetTotalSimTime() << " ("
            << OCPTIME_NRSTEP << "s)" << endl;
        cout << "     - % Updating Properties ......." << setw(fixWidth)
            << 100.0 * OCPTIME_UPDATE_GRID / control.simTime.GetTotalSimTime() << " ("
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
            MPI_Gather(record_local.data(), record_var_num, OCPMPI_DBL, record_total.data(), record_var_num, OCPMPI_DBL, MASTER_PROCESS, domain.global_comm);

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
            MPI_Gather(record_local.data(), record_var_num, OCPMPI_DBL, record_total.data(), record_var_num, OCPMPI_DBL, MASTER_PROCESS, domain.global_comm);
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