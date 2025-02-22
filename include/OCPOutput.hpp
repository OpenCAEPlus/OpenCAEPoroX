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
    /// index along x,y,z direction
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
    /// if use
    OCP_BOOL  activity{OCP_FALSE};
    /// object under current item
    vector<T> obj;
    /// Records the index of bulk or well, whose properties will be printed
    vector<USI> index; 
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
    /// Setup outputting structure
    void Setup(const OutputSummary& summary_param, const Reservoir& reservoir);

    /// Set value for vars
    void SetVal(const Reservoir& reservoir, const OCPControl& ctrl, const ItersInfo& iters, GetWallTime& timer);

    /// Write output information to a file.
    void PrintInfo(const string& dir, const string& filename, const OCP_INT& rank) const;

    /// Combine all files into 1 by Master process
    void PostProcess(const string& dir, const string& filename, const OCP_INT& numproc) const;

protected:
    void SetupOutputTerm(const OutputSummary& summary_param);

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

/// ToDo
class OutGridVar
{
    friend class OutGridVarSet;

protected:
    string  name;
         
    USI     gap;
    USI     offset;
};


/// Basic grid properties for output
class OutGridVarSet
{
    friend class Out4RPT;
    friend class Out4VTK;

public:
    void Setup(const OutGridParam& param, const Bulk& bk);

protected:
    USI      nc;                  ///< num of components
    USI      np;                  ///< num of phases
     
protected:
    OCP_BOOL DEPTH{OCP_FALSE};    ///< Depth
    OCP_BOOL PRE{OCP_FALSE};      ///< Pressure of grids.
    OCP_BOOL PHASEP{OCP_FALSE};   ///< Pressure of phases.
    OCP_BOOL COMPM{OCP_FALSE};    ///< Moles of components of grids.
    OCP_BOOL SOIL{OCP_FALSE};     ///< Oil saturation of grids.
    OCP_BOOL SGAS{OCP_FALSE};     ///< Gas saturation of grids.
    OCP_BOOL SWAT{OCP_FALSE};     ///< Water saturation of grids.
    OCP_BOOL DENO{OCP_FALSE};     ///< Oil density of grids.
    OCP_BOOL DENG{OCP_FALSE};     ///< Gas density of grids.
    OCP_BOOL DENW{OCP_FALSE};     ///< Water density of grids.
    OCP_BOOL KRO{OCP_FALSE};      ///< Oil relative permeability of grids.
    OCP_BOOL KRG{OCP_FALSE};      ///< Gas relative permeability of grids.
    OCP_BOOL KRW{OCP_FALSE};      ///< Water relative permeability of grids.
    OCP_BOOL BOIL{OCP_FALSE};     ///< Oil reservoir molar densities of grids.
    OCP_BOOL BGAS{OCP_FALSE};     ///< Gas reservoir molar densities of grids.
    OCP_BOOL BWAT{OCP_FALSE};     ///< Water reservoir molar densities of grids.
    OCP_BOOL VOIL{OCP_FALSE};     ///< Oil viscosity of grids.
    OCP_BOOL VGAS{OCP_FALSE};     ///< Gas viscosity of grids.
    OCP_BOOL VWAT{OCP_FALSE};     ///< Water viscosity of grids.
    OCP_BOOL XMF{OCP_FALSE};      ///< liquid component mole fractions.
    OCP_BOOL YMF{OCP_FALSE};      ///< gas component mole fractions.
    OCP_BOOL PCW{OCP_FALSE};      ///< capillary pressure: Po - Pw.
    OCP_BOOL CO2{ OCP_FALSE };    ///< CO2 Concentration in water phase
    OCP_BOOL SATNUM{ OCP_FALSE }; ///< SAT region
    OCP_BOOL PERMX{ OCP_FALSE };  ///< permeability of rock in x-direction
    OCP_BOOL PERMY{ OCP_FALSE };  ///< permeability of rock in y-direction
    OCP_BOOL PERMZ{ OCP_FALSE };  ///< permeability of rock in z-direction
    OCP_BOOL DSAT{ OCP_FALSE };   ///< saturation change of phases
    OCP_BOOL DP{ OCP_FALSE };    ///< pressure change of phases
    OCP_BOOL CSFLAG{ OCP_FALSE }; ///< flag of the process group of coupled solution
    OCP_BOOL ITERNRDDM{ OCP_FALSE }; ///< Accumulated iters for NR of DDM
    OCP_BOOL ITERLSDDM{ OCP_FALSE }; ///< Accumulated iters for LS of DDM
    OCP_BOOL TIMELSDDM{ OCP_FALSE }; ///< Accumulated time for LS of DDM

    USI      bgpnum;              ///< num of Basic grid information to be printed
};

/// Collect more detailed information of each time step.
//class Out4RPT
//{
//public:
//    void InputParam(const OutputRPTParam& RPTparam);
//    void Setup(const string& dir, const Reservoir& reservoir);
//    void PrintRPT(const string& dir, const Reservoir& rs, const OCP_DBL& days) const;
//    void GetIJKGrid(USI& i, USI& j, USI& k, const OCP_ULL& n) const;
//
//private:
//    OCP_BOOL          useRPT{OCP_FALSE};
//    OCP_ULL           numGrid;
//    OCP_USI           nx;
//    OCP_USI           ny;
//    USI               IJKspace;
//    OutGridVarSet     bgp;
//};


class Out4VTK
{
public:
    /// Input Params about vtk output
    void Setup(const OutputVTKParam& VTKParam, const string& dir, const Reservoir& rs);
    /// Output info with the vtk format
    void PrintVTK(const Reservoir& rs, const OCPControl& ctrl) const;
    /// Combine all files into 1 by Master process
    void PostProcess(const string& dir, const string& filename, const OCP_INT& numproc) const;

protected:
    void Initialize(const string& dir, const Reservoir& rs);
    void PostProcessP(const string& dir, const string& filename, const OCP_INT& numproc) const;
    void PostProcessS(const string& dir, const string& filename) const;

protected:
    OCP_BOOL          useVTK{OCP_FALSE}; ///< If use vtk
    OutGridVarSet     bgp;               ///< Basic grid information
    Output4Vtk        out4vtk;           ///< Output for vtk
    string            myFile;            ///< output file name
    mutable USI       countPrint{ 0 };   ///< record the count printed
    
    mutable OCP_CHAR* timeInfo;

    /// total number of grids
    mutable OCP_ULL   numGrid;
};

/// The OCPOutput class manages different kinds of ways to output information.
//  Note: The most commonly used is the summary file, which usually gives the
//  information of bulks and wells in each time step, such as average pressure, oil
//  production rate of wells. If other information at critical dates is of interest, you
//  can chose the PRT/VTK file. Also, some infomation will be printed on the screen
//  at the critical dates to make sure the program is at the right way.
class OCPOutput
{
    friend class OpenCAEPoroX;

public:
    /// Input Params about output
    void Setup(const ParamOutput& paramOutput, const OCPControl& ctrl, const Reservoir& rs);
    /// Assign values to be output in PrintAtTimeStep()
    void SetValAtTimeStep(const Reservoir& rs, const OCPControl& ctrl, const OCPNRsuite& NR, GetWallTime& timer_total);
    /// Output info which is each time step based
    void PrintAtTimeStep() const;
    /// Output info which is Keyword TSTEP based
    void PrintInfoSched(const Reservoir& rs, const OCPControl& ctrl, const OCP_DBL& time) const;
    /// Post-process the output file
    void PostProcess() const;
    /// output current time step and iterations to screen
    void PrintCurrentTimeIter(const OCPControl& ctrl) const;

protected:
    /// Setup communicator
    void SetupComm(const Domain& domain);

protected:
    MPI_Comm  myComm{ MPI_COMM_NULL };
    OCP_INT   numproc, myrank;


protected:
    string       workDir;
    string       fileName;
    ItersInfo    iters;
    Summary      summary;
    CriticalInfo crtInfo;
    Out4VTK      out4VTK;
    // Out4RPT      out4RPT;
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