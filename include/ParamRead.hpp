/*! \file    ParamRead.hpp
 *  \brief   ParamRead class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *  \note    The params used in OpenCAEPoroX is mostly compatible with Eclipse by SLB,
 *           but it has some own rules for easy to use. It is extensible and friendly.
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARAMREAD_HEADER__
#define __PARAMREAD_HEADER__

// Standard header files
#include <fstream>
#include <iostream>
#include <string>

// OpenCAEPoroX header files
#include "ParamControl.hpp"
#include "ParamOutput.hpp"
#include "ParamReservoir.hpp"
#include "ParamWell.hpp"

using namespace std;

/// Pre-processing unit for OpenCAEPoroX for reading params from input files.
class ParamRead
{
public:
    string   inputFile; ///< Input file with its path (absolute or relative).
    string   workDir;   ///< Current work directory.
    string   fileName;  ///< File name of input file.
    OCP_BOOL disable_grid{ OCP_FALSE }; ///< grid-based params are ignored

    // Main workloads for ParamRead: read reservoir params, read well params,
    // read control params, and read output param. These workloads are assigned
    // to the following classes, respectively.
    ParamReservoir paramRs;      ///< Read the reservoir params.
    ParamWell      paramWell;    ///< Read the well params.
    ParamControl   paramControl; ///< Read the control params.
    ParamOutput    paramOutput;  ///< Read the output params.

public:

    ParamRead(const OCP_BOOL& disableGrid = OCP_FALSE) {
        disable_grid         = disableGrid;
        paramRs.disable_grid = disable_grid;
    };

    /// Get current work dir and input file name from the full file path.
    void GetDirAndName();

    /// Initialize the param reading process.
    void Init();

    /// General interface for reading input data.
    void ReadInputFile(const string& file);

    /// Read the input file.
    void ReadFile(const string& file);

    /// Handle the INCLUDE keyword, which contains other input files.
    void ReadINCLUDE(ifstream& ifs);

    /// Post-processing
    void PostPrecess();

    /// Check whether the params contain error.
    void CheckParam();
};

#endif /* end if __PARAMREAD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/16/2021      Format file                          */
/*----------------------------------------------------------------------------*/