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

using namespace std;


/// All parameters used for solution control
class OCPControl
{
public:
    /// Input parameters for control.
    void InputParam(const ParamControl& CtrlParam);
    /// Setup Comm
    void Setup(const Domain& domain);
    /// Apply control for time step i.
    void ApplyControl(const USI& i, const OCP_BOOL& wellOptChange_loc);
    /// Check if converge
    OCPNRStateC CheckConverge(const OCPNRsuite& NRs, const initializer_list<string>& il) {
        return NR.CheckConverge(NRs, il);
    }
    // Calculate next time step
    void CalNextTimeStep(const OCPNRsuite& NRs, const initializer_list<string>& il) {
        time.CalNextTimeStep(NRs, il);
    }

public:
    MPI_Comm         myComm;
    OCP_INT          numproc, myrank;

public:
    /// Print level
    USI           printLevel{0};
    /// Time control 
    ControlTime   time;
    /// NR control    
    ControlNR     NR;
    /// Solver Method control
    ControlMethod SM;
    /// Stop simulation
    OCP_BOOL    StopSim{ OCP_FALSE };
    /// Stop time
    OCP_DBL     MaxSimTime{ 1E20 };


public:  
    /// Get model
    auto GetModel() const { return model; }
    /// Get type of the solution method.
    auto GetMethod() const { return method; }
    /// Get work dir name.
    auto GetWorkDir() const { return workDir; }
    /// Get OCP file name.
    auto GetOCPFile() const { return ocpFile; }
    /// Get linear solver file name.
    auto GetLsFile() const { return lsFile; }
    /// Setup fast Control.
    void SetupFastControl(const USI& argc, const char* optset[]);

protected:
    /// model: isothermal, thermal
    OCPModel            model{ OCPModel::none };
    /// Discrete method
    vector<OCPNLMethod> method;
    /// File name of linear Solver
    vector<string>      lsFile;
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