/*! \file    ParamOutput.hpp
 *  \brief   ParamOutput class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARAMOUTPUT_HEADER__
#define __PARAMOUTPUT_HEADER__

// Standard header files
#include <fstream>
#include <vector>

// OpenCAEPoroX header files
#include "UtilInput.hpp"

/// A structure of three-dimensional coordinates.
class COOIJK
{
public:
    COOIJK() = default;
    COOIJK(const USI& i, const USI& j, const USI& k)
        : I(i)
        , J(j)
        , K(k){};
    USI I;
    USI J;
    USI K;
};

/// Used to store the contents of keyword whose contents are in form of coordinates.
class Type_B_o
{
public:
    OCP_BOOL       activity{OCP_FALSE};
    vector<COOIJK> obj;
};

/// Used to store the contents of keyword whose contents are in form of string.
class Type_A_o
{
public:
    OCP_BOOL       activity{OCP_FALSE};
    vector<string> obj;
};

/// OutputSummary contains all information about SUMMARY, these information tells
/// which results will be print to the summary output file.
class OutputSummary
{

public:
    OCP_BOOL FPR{OCP_FALSE};  ///< Field average Pressure.
    OCP_BOOL FTR{OCP_FALSE};  ///< Field average Temperature.
    OCP_BOOL FOPR{OCP_FALSE}; ///< Field oil production rate.
    OCP_BOOL FOPT{OCP_FALSE}; ///< Field total oil production.
    OCP_BOOL FGPR{OCP_FALSE}; ///< Field gas production rate.
    OCP_BOOL FGPt{OCP_FALSE}; ///< Field total gas production.
    OCP_BOOL FWPR{OCP_FALSE}; ///< Field water production rate.
    OCP_BOOL FWPT{OCP_FALSE}; ///< Field total water production.
    OCP_BOOL FGIR{OCP_FALSE}; ///< Field gas injection rate.
    OCP_BOOL FGIT{OCP_FALSE}; ///< Field total gas injection.
    OCP_BOOL FWIR{OCP_FALSE}; ///< Field water injection rate.
    OCP_BOOL FWIT{OCP_FALSE}; ///< Field total water injection.

    Type_A_o WOPR; ///< Well oil production rate.
    Type_A_o WOPT; ///< Well total oil production rate.
    Type_A_o WGPR; ///< Well gas production rate.
    Type_A_o WGPT; ///< Well total gas production.
    Type_A_o WWPR; ///< Well water production rate.
    Type_A_o WWPT; ///< Well total water production.
    Type_A_o WGIR; ///< Well gas injection rate.
    Type_A_o WGIT; ///< Well total gas injection.
    Type_A_o WWIR; ///< Well water injection rate.
    Type_A_o WWIT; ///< Well total water injection.
    Type_A_o WBHP; ///< Well pressure.
    Type_A_o DG;   ///< Pressure difference between wells and perforations.

    Type_B_o BPR;  ///< Bulk pressure.
    Type_B_o SOIL; ///< Oil saturation of bulk.
    Type_B_o SGAS; ///< Gas saturation of bulk.
    Type_B_o SWAT; ///< Water saturation of bulk.
};

class OutGridParam
{
public:
    OCP_BOOL ASCII{OCP_FALSE};   ///< If output using ASCII format
    OCP_DBL  DOUBLE{OCP_FALSE};  ///< If use double precision for float number
    OCP_BOOL DEPTH{OCP_FALSE};   ///< depth of grids.
    OCP_BOOL PRE{OCP_FALSE};     ///< Pressure of grids.
    OCP_BOOL PHASEP{OCP_FALSE};  ///< Pressure of phases.
    OCP_BOOL COMPM{OCP_FALSE};   ///< Moles of components of grids.
    OCP_BOOL SOIL{OCP_FALSE};    ///< Oil saturation of grids.
    OCP_BOOL SGAS{OCP_FALSE};    ///< Gas saturation of grids.
    OCP_BOOL SWAT{OCP_FALSE};    ///< Water saturation of grids.
    OCP_BOOL DENO{OCP_FALSE};    ///< Oil density of grids.
    OCP_BOOL DENG{OCP_FALSE};    ///< Gas density of grids.
    OCP_BOOL DENW{OCP_FALSE};    ///< Water density of grids.
    OCP_BOOL KRO{OCP_FALSE};     ///< Oil relative permeability of grids.
    OCP_BOOL KRG{OCP_FALSE};     ///< Gas relative permeability of grids.
    OCP_BOOL KRW{OCP_FALSE};     ///< Water relative permeability of grids.
    OCP_BOOL BOIL{OCP_FALSE};    ///< Oil reservoir molar densities of grids.
    OCP_BOOL BGAS{OCP_FALSE};    ///< Gas reservoir molar densities of grids.
    OCP_BOOL BWAT{OCP_FALSE};    ///< Water reservoir molar densities of grids.
    OCP_BOOL VOIL{OCP_FALSE};    ///< Oil viscosity of grids.
    OCP_BOOL VGAS{OCP_FALSE};    ///< Gas viscosity of grids.
    OCP_BOOL VWAT{OCP_FALSE};    ///< Water viscosity of grids.
    OCP_BOOL XMF{OCP_FALSE};     ///< liquid component mole fractions.
    OCP_BOOL YMF{OCP_FALSE};     ///< gas component mole fractions.
    OCP_BOOL PCW{OCP_FALSE};     ///< capillary pressure: Po - Pw.
    OCP_BOOL CO2{OCP_FALSE};     ///< CO2 Concentration
    OCP_BOOL SATNUM{OCP_FALSE};  ///< SAT region
    OCP_BOOL PERMX{OCP_FALSE};   ///< permeability of rock in x-direction
    OCP_BOOL PERMY{OCP_FALSE};   ///< permeability of rock in y-direction
    OCP_BOOL PERMZ{OCP_FALSE};   ///< permeability of rock in z-direction
    OCP_BOOL DSAT{OCP_FALSE};    ///< saturation change of phases
    OCP_BOOL DP{ OCP_FALSE };    ///< pressure change of phases
    OCP_BOOL CSFLAG{OCP_FALSE};  ///< flag of the process group of coupled solution
    OCP_BOOL ITERNRDDM{ OCP_FALSE }; ///< Accumulated iters for NR of DDM
    OCP_BOOL ITERLSDDM{ OCP_FALSE }; ///< Accumulated iters for LS of DDM
    OCP_BOOL TIMELSDDM{ OCP_FALSE }; ///< Accumulated time for LS of DDM
};

/// OutputRPTParam is a part of ParamOutput, it's used to control the output of detailed
/// information of reservoir. For example, pressure in every grid, due to excessive data
/// volume, this option is usually used in Debug Model.
class OutputRPTParam
{
public:
    OCP_BOOL               useRPT{OCP_FALSE};
    OutGridParam bgp;
};

class OutputVTKParam
{
public:
    OCP_BOOL               useVTK{OCP_FALSE};
    OutGridParam bgp;
};

/// ParamOutput is an internal structure used to stores the information of outputting
/// from input files. It is an intermediate interface and independent of the main
/// simulator. After all file inputting finishes, the params in it will pass to
/// corresponding modules.
class ParamOutput
{
public:
    OutputSummary  summary;     ///< See OutputSummary.
    OutputRPTParam outRPTParam; ///< See OutputRPTParam.
    OutputVTKParam outVTKParam; ///< See OutputVTKParam

    /// Input the keyword SUMMARY, which contains many sub-keyword, indicating which
    /// results are interested by user. After the simulation, these results will be
    /// output into a summary file.
    void InputSUMMARY(ifstream& ifs);
    /// Input the sub-keyword in SUMMARY, the contents in these keyword is in the form
    /// of string.
    void InputType_A(ifstream& ifs, Type_A_o& obj);
    /// Input the sub-keyword in SUMMARY, the contents in these keyword is in the form
    /// of coordinates.
    void InputType_B(ifstream& ifs, Type_B_o& obj);

    /// Input the keyword RPTSCHED, which tells which detailed information will be
    /// output to the RPTfile.
    void InputRPTSCHED(ifstream& ifs, const string& keyword);
};

#endif /* end if __PARAMOUTPUT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/