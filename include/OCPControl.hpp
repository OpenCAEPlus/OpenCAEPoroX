/*! \file    OCPControl.hpp
 *  \brief   OCPControl class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPCONTROL_HEADER__
#define __OCPCONTROL_HEADER__

// Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "ParamControl.hpp"
#include "OCPControlTime.hpp"
#include "OCPControlNR.hpp"
#include "OCPControlMethod.hpp"
#include "OCPControlSimTime.hpp"

using namespace std;


/// All parameters used for solution control
class OCPControl
{
public:
    /// Input parameters for control.
    void Setup(const USI& argc, const char* argv[], const ParamControl& CtrlParam, const Domain& domain);
    /// Apply control for time step i.
    INT ApplyControl();
    /// Check if converge
    OCPNRStateC CheckConverge(const OCPNRsuite& NRs, const initializer_list<string>& il) {
        return NR.CheckConverge(NRs, il);
    }
    // Calculate next time step
    void CalNextTimeStep(const OCPNRsuite& NRs, const initializer_list<string>& il) {
        time.CalNextTimeStep(NRs, il);
    }

protected:
    /// Setup communicator
    void SetupComm(const Domain& domain);

public:
    MPI_Comm         myComm;
    OCP_INT          numproc, myrank;

public:
    /// Print level
    USI            printLevel{0};
    /// Time control 
    ControlTime    time;
    /// NR control    
    ControlNR      NR;
    /// Solver Method control
    ControlMethod  SM;
    /// simluation time control
    ControlSimTime simTime;
    /// Stop simulation
    OCP_BOOL    StopSim{ OCP_FALSE };


public:  
    /// Get work dir name.
    auto GetWorkDir() const { return workDir; }
    /// Get OCP file name.
    auto GetOCPFile() const { return ocpFile; }
    /// Output model information
    void OutputModelMethodInfo() const;

protected:
    /// Setup fast Control.
    void SetFastControl(const USI& argc, const char* optset[]);

protected:
    /// Current work directory
    string              workDir;
    /// Current file name
    string              ocpFile;

};

#endif /* end if __OCPControl_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/08/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/