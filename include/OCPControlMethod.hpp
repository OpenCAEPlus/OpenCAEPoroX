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
    /// Initialize calling sequence of methods
    OCPNLMethod InitMethod() const;
    /// Switch to main method
    OCPNLMethod SwitchMethod() const;
    /// Output model and method info
    void OutputInfo(const USI& pl) const;
    /// Get model
    auto GetModel() const { return model; }
    /// Get type of the solution method.
    auto GetMethod() const { return method; }
    /// Get ith ls file
    auto GetLsFile(const USI& i) const { return lsFile[i]; }
    /// Get work dir name.
    auto GetWorkDir() const { return workDir; }

protected:
    /// work directory
    string              workDir;
    /// model: isothermal, thermal
    OCPModel            model{ OCPModel::none };
    /// Index of working method
    mutable USI         wIndex;
    /// main method is placed in the front following by the preconditioner method sequentially,
    /// only at most two methods' cooperation is supported now
    vector<OCPNLMethod> method;
    /// File name of linear Solver
    vector<string>      lsFile;
};

#endif /* end if __OCPControlMethod_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/05/2023      Create file                          */
/*----------------------------------------------------------------------------*/