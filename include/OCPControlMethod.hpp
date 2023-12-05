/*! \file    OCPControlMethod.hpp
 *  \brief   OCPControlMethod class declaration
 *  \author  Shizhe Li
 *  \date    Dec/05/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPCONTROLMETHOD_HEADER__
#define __OCPCONTROLMETHOD_HEADER__

 // Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"


using namespace std;


/// All parameters used for solution control
class OCPControlMethod
{
    friend class OCPControl;
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

#endif /* end if __OCPControlMethod_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/05/2023      Create file                          */
/*----------------------------------------------------------------------------*/