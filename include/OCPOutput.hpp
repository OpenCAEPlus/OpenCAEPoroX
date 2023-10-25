/*! \file    OCPOutput.hpp
 *  \brief   OCPOutput class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCP_OUTPUT_HEADER__
#define __OCP_OUTPUT_HEADER__

// Standard header files
#include <iomanip>
#include <iostream>
#include <cstdio>

// OpenCAEPoroX header files
#include "OCPControl.hpp"
#include "Output4Vtk.hpp"
#include "ParamOutput.hpp"
#include "Reservoir.hpp"
#include "UtilOutput.hpp"
#include "UtilTiming.hpp"


using namespace std;

/// 3D coordinate representation in OpenCAEPoroX
class OCPIJK
{
public:
    OCPIJK() = default;
    OCPIJK(const USI& i, const USI& j, const USI& k)
        : I(i)
        , J(j)
        , K(k){};
    OCPIJK(const COOIJK& src)
    {
        I = src.I;
        J = src.J;
        K = src.K;
    };
    OCPIJK& operator=(const COOIJK& src)
    {
        I = src.I;
        J = src.J;
        K = src.K;
        return *this;
    }
    USI I, J, K;
};

/// TODO: Add Doxygen
template <typename T>
class OCPType_Sum
{
public:
    OCPType_Sum& operator=(const Type_A_o& src)
    {
        activity = src.activity;
        obj      = src.obj;
        return *this;
    }
    OCPType_Sum& operator=(const Type_B_o& src)
    {
        activity = src.activity;
        obj.assign(src.obj.begin(), src.obj.end());
        return *this;
    }
    OCP_BOOL  activity{OCP_FALSE};
    vector<T> obj;
    vector<USI>
        index; ///< Records the index of bulk or well, whose properties will be printed
};


/// Record iteration information
class ItersInfo
{

public:
    /// Update all Iters
    void Update(const OCPNRsuite& NRs);
    /// Get num of time steps
    auto GetNumTimeStep() const { return numTstep; }
    /// Return the total number of Newton iterations.
    auto GetNRt() const { return NRt; }
    /// Return the total number of wasted Newton iterations.
    auto GetNRwt() const { return NRwt; }
    /// Return the total number of linear iterations.
    auto GetLSt() const { return LSt; }
    /// Return the total number of wasted linear iterations.
    auto GetLSwt() const { return LSwt; }

protected:
    /// number of time step
    USI numTstep{ 0 };
    /// total number of Newton iterations
    USI NRt{ 0 };
    /// total number of wasted Newton iterations
    USI NRwt{ 0 };
    /// total number of iterations of linear solver
    USI LSt{ 0 };
    /// totalnumber of wasted linear iterations
    USI LSwt{ 0 };
};


/// The SumItem class is an auxiliary structure storing summary data to output.
class SumItem
{
public:
    SumItem() = default;
    SumItem(const string& item, const string& obj) :Item(item), Obj(obj) {};
    SumItem(const string& item, const USI& n) : Item(item) { val.reserve(n); };
    SumItem(const string& item,
            const string& obj,
            const string& unit,
            const string& type,
            const OCP_USI& n)
        : Item(item)
        , Obj(obj)
        , Unit(unit)
        , Type(type) { val.reserve(n); };
    OCP_BOOL operator==(const SumItem& other) {
        if (this->Item == other.Item && this->Obj == other.Obj)    return OCP_TRUE;
        else                                                       return OCP_FALSE;
    }
    string          Item;
    string          Obj;
    string          Unit;
    string          Type;
    vector<OCP_DBL> val;
};

/// The Summary class manages the output in the summary file.
//  Note: It contains the most interested information in each time step, which usually
//  will be convert to figures for later analysis.
class Summary
{
public:
    /// Input var names those will be printed
    void InputParam(const OutputSummary& summary_param);

    /// Setup outputting structure
    void Setup(const Reservoir& reservoir);

    /// Set value for vars
    void SetVal(const Reservoir& reservoir, const OCPControl& ctrl, const ItersInfo& iters);

    /// Write output information to a file.
    void PrintInfo(const string& dir, const string& filename, const OCP_INT& rank) const;

    /// Combine all files into 1 by Master process
    void PostProcess(const string& dir, const string& filename, const OCP_INT& numproc) const;

protected:
    vector<SumItem> Sumdata; ///< Contains all information to be printed.

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

    OCPType_Sum<string> WOPR; ///< Well oil production rate.
    OCPType_Sum<string> WOPT; ///< Well total oil production.
    OCPType_Sum<string> WGPR; ///< Well gas production rate.
    OCPType_Sum<string> WGPT; ///< Well total gas production.
    OCPType_Sum<string> WWPR; ///< Well water production rate.
    OCPType_Sum<string> WWPT; ///< Well total water production.
    OCPType_Sum<string> WGIR; ///< Well gas injection rate.
    OCPType_Sum<string> WGIT; ///< Well total gas injection.
    OCPType_Sum<string> WWIR; ///< Well water injection rate.
    OCPType_Sum<string> WWIT; ///< Well total water injection.
    OCPType_Sum<string> WBHP; ///< Well pressure.
    OCPType_Sum<string> DG;   ///< Pressure difference between wells and perforations.

    OCPType_Sum<OCPIJK> BPR;  ///< Bulk pressure.
    OCPType_Sum<OCPIJK> SOIL; ///< Oil saturation of bulk.
    OCPType_Sum<OCPIJK> SGAS; ///< Gas saturation of bulk.
    OCPType_Sum<OCPIJK> SWAT; ///< Water saturation of bulk.
};

/// Collect important information of each time step for fast review.
class CriticalInfo
{
public:
    /// Setup vars to print and reserve memory
    void Setup();

    /// Set values for vars
    void SetVal(const OCPControl& ctrl, const OCPNRsuite& NR);

    /// Print to file
    void PrintFastReview(const string& dir, const string& filename, const OCP_INT& rank, const ItersInfo& iters) const;

    /// Combine all files into 1 by Master process
    void PostProcess(const string& dir, const string& filename, const OCP_INT& numproc, const ItersInfo& iters) const;

private:
    vector<SumItem> Sumdata;
};

/// Basic grid properties for output
class BasicGridProperty
{
    friend class Out4RPT;
    friend class Out4VTK;

public:
    void SetBasicGridProperty(const BasicGridPropertyParam& param);
    void Check(const Reservoir& rs);

private:
    OCP_BOOL PRE{OCP_FALSE};  ///< Pressure of grids.
    OCP_BOOL PGAS{OCP_FALSE}; ///< Gas pressure of grids.
    OCP_BOOL PWAT{OCP_FALSE}; ///< Water pressure of grids.
    OCP_BOOL SOIL{OCP_FALSE}; ///< Oil saturation of grids.
    OCP_BOOL SGAS{OCP_FALSE}; ///< Gas saturation of grids.
    OCP_BOOL SWAT{OCP_FALSE}; ///< Water saturation of grids.
    OCP_BOOL DENO{OCP_FALSE}; ///< Oil density of grids.
    OCP_BOOL DENG{OCP_FALSE}; ///< Gas density of grids.
    OCP_BOOL DENW{OCP_FALSE}; ///< Water density of grids.
    OCP_BOOL KRO{OCP_FALSE};  ///< Oil relative permeability of grids.
    OCP_BOOL KRG{OCP_FALSE};  ///< Gas relative permeability of grids.
    OCP_BOOL KRW{OCP_FALSE};  ///< Water relative permeability of grids.
    OCP_BOOL BOIL{OCP_FALSE}; ///< Oil reservoir molar densities of grids.
    OCP_BOOL BGAS{OCP_FALSE}; ///< Gas reservoir molar densities of grids.
    OCP_BOOL BWAT{OCP_FALSE}; ///< Water reservoir molar densities of grids.
    OCP_BOOL VOIL{OCP_FALSE}; ///< Oil viscosity of grids.
    OCP_BOOL VGAS{OCP_FALSE}; ///< Gas viscosity of grids.
    OCP_BOOL VWAT{OCP_FALSE}; ///< Water viscosity of grids.
    OCP_BOOL XMF{OCP_FALSE};  ///< liquid component mole fractions.
    OCP_BOOL YMF{OCP_FALSE};  ///< gas component mole fractions.
    OCP_BOOL PCW{OCP_FALSE};  ///< capillary pressure: Po - Pw.

    USI             bgpnum;   ///< num of Basic grid information to be printed
};

/// Collect more detailed information of each time step.
class Out4RPT
{
public:
    void InputParam(const OutputRPTParam& RPTparam);
    void Setup(const string& dir, const Reservoir& reservoir);
    void PrintRPT(const string& dir, const Reservoir& rs, const OCP_DBL& days) const;
    template <typename T>
    void PrintRPT_Scalar(ofstream&              ifs,
                         const string&          dataName,
                         const OCP_DBL&         days,
                         const T*               gridVal,
                         const USI&             gap,
                         const vector<GB_Pair>& gbPair,
                         const bool&            useActive,
                         const OCP_DBL&         alpha = 1.0) const;
    void GetIJKGrid(USI& i, USI& j, USI& k, const OCP_ULL& n) const;

private:
    OCP_BOOL          useRPT{OCP_FALSE};
    OCP_ULL           numGrid;
    OCP_USI           nx;
    OCP_USI           ny;
    USI               IJKspace;
    BasicGridProperty bgp;
};

template <typename T>
void Out4RPT::PrintRPT_Scalar(ofstream&              myRPT,
                              const string&          dataName,
                              const OCP_DBL&         days,
                              const T*               gridVal,
                              const USI&             gap,
                              const vector<GB_Pair>& gbPair,
                              const bool&            useActive,
                              const OCP_DBL&         alpha) const
{
    USI     I, J, K;
    OCP_ULL bId;

    myRPT << OCP_SEP01(50) << "\n";
    myRPT << dataName << "                   " << fixed << setprecision(3) << days
          << "  DAYS";

    if (useActive) {
        for (OCP_ULL n = 0; n < numGrid; n++) {
            if (n % nx == 0) myRPT << "\n";
            if (n % (nx * ny) == 0) myRPT << "\n\n";

            if (n % nx == 0) {
                GetIJKGrid(I, J, K, n);
                myRPT << GetIJKformat("*", to_string(J), to_string(K), IJKspace);
            }

            if (gbPair[n].IsAct()) {
                bId = gbPair[n].GetId();
                myRPT << setw(10) << fixed << setprecision(3)
                      << gridVal[bId * gap] * alpha;
            } else {
                myRPT << setw(10) << " --- ";
            }
        }
    } else {
        for (OCP_ULL n = 0; n < numGrid; n++) {
            if (n % nx == 0) myRPT << "\n";
            if (n % (nx * ny) == 0) myRPT << "\n\n";

            if (n % nx == 0) {
                GetIJKGrid(I, J, K, n);
                myRPT << GetIJKformat("*", to_string(J), to_string(K), IJKspace);
            }
            myRPT << setw(10) << fixed << setprecision(3) << gridVal[n * gap] * alpha;
        }
    }

    myRPT << "\n\n\n";
}

class Out4VTK
{
public:
    /// Input Params about vtk output
    void InputParam(const OutputVTKParam& VTKParam);
    /// Setup up vtk output
    void Setup(const string& dir, const Reservoir& rs);
    /// Output info with the vtk format
    void PrintVTK(const Reservoir& rs) const;
    /// Combine all files into 1 by Master process
    void PostProcess(const string& dir, const string& filename, const OCP_INT& numproc) const;
    void PostProcessP(const string& dir, const string& filename, const OCP_INT& numproc) const;
    void PostProcessS(const string& dir, const string& filename) const;

private:
    OCP_BOOL          useVTK{OCP_FALSE}; ///< If use vtk
    BasicGridProperty bgp;               ///< Basic grid information
    Output4Vtk        out4vtk;           ///< Output for vtk
    string            myFile;            ///< output file name

    /// total number of grids
    mutable OCP_ULL   numGrid;
};

/// The OCPOutput class manages different kinds of ways to output information.
//  Note: The most commonly used is the summary file, which usually gives the
//  information of bulks and wells in each time step, such as average pressure, oil
//  production rate of wells. If other information at critical dates is of interest, you
//  can chose the PRT file (TODO). Also, some infomation will be printed on the screen
//  at the critical dates to make sure the program is at the right way.
class OCPOutput
{
    friend class OpenCAEPoroX;

public:
    /// Input Params about output
    void InputParam(const ParamOutput& paramOutput);
    /// Setup all kinds of output
    void Setup(const Reservoir& reservoir, const OCPControl& ctrl, const Domain& domain);
    /// Assign values to be output in PrintInfo()
    void SetVal(const Reservoir& reservoir, const OCPControl& ctrl, const OCPNRsuite& NR);
    /// Output info which is each time step based
    void PrintInfo() const;
    /// Output info which is Keyword TSTEP based
    void PrintInfoSched(const Reservoir&  rs, const OCPControl& ctrl, const OCP_DBL& time) const;
    /// Post-process the output file
    void PostProcess() const;
    /// output current time step and iterations to screen
    void PrintCurrentTimeIter(const OCPControl& ctrl) const;

protected:
    MPI_Comm  myComm;
    OCP_INT   numproc, myrank;


protected:
    string       workDir;
    string       fileName;
    ItersInfo    iters;
    Summary      summary;
    CriticalInfo crtInfo;
    Out4RPT      out4RPT;
    Out4VTK      out4VTK;
};

#endif /* end if __OCPOUTPUT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/08/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/