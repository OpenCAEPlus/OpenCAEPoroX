/*! \file    Main.cpp
 *  \brief   An example to demonstrate main steps of the OCP simulator
 *  \author  Shizhe Li
 *  \date    Feb/15/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// Standard header files
#include <cstdio>
#include <iostream>
#include <string>
#include <mpi.h>

// OpenCAEPoroX header files
#include "OCP.hpp"
#include "PreProcess.hpp"

using namespace std;

/// The main() function performs dynamic simulation in five steps.
int main(int argc, char* argv[])
{

    OCP_INT myRank, commSize;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    CURRENT_RANK = myRank;

    OpenCAEPoroX simulator;

    // Step 0. Print simulator version information.   
    if (myRank == MASTER_PROCESS) {
        if (argc < 2) {
            simulator.PrintUsage(argv[0]);
            return OCP_ERROR_NUM_INPUT; // Need at least one parameter
        }
        else {
            simulator.PrintVersion();
            if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
                simulator.PrintUsage(argv[0]);
                return OCP_SUCCESS;
            }
        }
    }

    {
        // Step 1. Input and generate Grid infomation and partition
        PreProcess preProcess(argv[1], myRank, MPI_COMM_WORLD);

        // Step 2. Input reservoir information and distribute
        simulator.InputDistParam(argv[1], preProcess, myRank);
    }

    // Step 3. Setup params
    simulator.SetupSimulator(argc, const_cast<const char**>(argv));

    // Step 4. Initialize the reservoir
    simulator.InitReservoir();

    // Step 4. Run dynamic simulation using methods like IMPEC, AIM, and FIM.
    simulator.RunSimulation();

    // Step 5. Output the results according to control params.
    simulator.OutputResults();

    MPI_Finalize();

    return OCP_SUCCESS;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Feb/15/2023      Create file                          */
/*----------------------------------------------------------------------------*/
