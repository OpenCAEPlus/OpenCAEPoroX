/*! \file    ControlMethod.hpp
 *  \brief   ControlMethod class declaration
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
#include "ParamControl.hpp"
#include "OCPControlFast.hpp"


using namespace std;


/// control the usage of solver method
class ControlMethod
{
    friend class OCPControl;
public:
    /// Set control param
    void SetCtrlParam(const ParamControl& CtrlParam);
    /// Set fast control
    void SetFastControl(const FastControl& fCtrl);
    /// Switch to main method
    auto SwitchMethod() const;
    /// Output model and method info
    void OutputInfo(const USI& pl) const;
    /// Get model
    auto GetModel() const { return model; }
    /// Get type of the solution method.
    auto GetMethod() const { return method; }
    /// Get ith methode
    auto GetMethod(const USI& i) const { return method[i]; }
    /// Get type of the main method.
    auto GetMainMethod() const { return mainMethod; }
    /// Get type of the preconditioner method.
    auto GetPreMethod() const { return preMethod; }
    /// Get linear solver file name.
    auto GetLsFile() const { return lsFile; }
    /// Get ith ls file
    auto GetLsFile(const USI& i) const { return lsFile[i]; }
    /// Get work dir name.
    auto GetWorkDir() const { return workDir; }

protected:
    /// work directory
    string              workDir;
    /// model: isothermal, thermal
    OCPModel            model{ OCPModel::none };
    /// Discrete method
    vector<OCPNLMethod> method;
    /// File name of linear Solver
    vector<string>      lsFile;
    /// main method
    OCPNLMethod         mainMethod{ OCPNLMethod::none };
    /// pre method
    OCPNLMethod         preMethod{ OCPNLMethod::none };
};

#endif /* end if __OCPControlMethod_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/05/2023      Create file                          */
/*----------------------------------------------------------------------------*/