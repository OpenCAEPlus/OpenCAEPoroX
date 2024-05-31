/*! \file    ParamWell.hpp
 *  \brief   ParamWell class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARAMWELL_HEADER__
#define __PARAMWELL_HEADER__

// Standard header files
#include <cassert>
#include <fstream>
#include <vector>
#include <map>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "UtilInput.hpp"
#include "UtilOutput.hpp"

using namespace std;

/// WellOptParam contains all the param used for well operation. for a well, all the
/// params in it can be changed over the time.
class WellOptParam
{
public:
    WellOptParam(string intype, vector<string>& vbuf);
    /// INJE well
    WellOptParam(string fluid_type, string state_, string mode_, string max_rate, OCP_DBL max_bhp);
    /// PROD well
    WellOptParam(string state_, string mode_, string max_rate, OCP_DBL min_bhp);
    /// 拷贝构造函数
    WellOptParam(const WellOptParam& opt);

    enum Type
    {
        INVALID=0,
        PROD,
        INJ,
    };

    // WCONINJE & WCONPROD
    string type;         ///< Type of well, injection or production
    string fluidType;    ///< Type of fluid into the injection well. (injection well only)
    string state;        ///< State of well, open or close?
    string mode;         ///< Mode of well, Rate or BHP?

    OCP_DBL maxRate;     ///< Maximum allowable flow rate into/out the well.
    OCP_DBL maxBHP;      ///< Maximum allowable pressure in the injection well.
    OCP_DBL minBHP{0.0};      ///< Minimum allowable pressure in the production well.
    OCP_DBL injTemp;     ///< Temperature of injected fluid.
};

static std::map<std::string, WellOptParam::Type> TypeMap = {
        {"WIR", WellOptParam::INJ},
        {"GIR", WellOptParam::INJ},
        {"WIBHP", WellOptParam::INJ},
        {"GIBHP", WellOptParam::INJ},
        {"WITHP", WellOptParam::INJ},
        {"GITHP", WellOptParam::INJ},
        //
        {"ORAT", WellOptParam::PROD},
        {"GRAT", WellOptParam::PROD},
        {"WRAT", WellOptParam::PROD},
        {"LRAT", WellOptParam::PROD},
        {"BHP", WellOptParam::PROD},
        {"THP", WellOptParam::PROD},
};

/// WellOptPair contains two parts, one is the operation mode of well, the other is the
/// beginning time at which the operation is applied. The beginning time is represented
/// by an index in critical time.
class WellOptPair
{
public:
    WellOptPair(USI i, string type, vector<string>& vbuf)
        : d(i)
        , opt(type, vbuf){};
    WellOptPair(USI i, WellOptParam opt_): d(i), opt(opt_) {}

    USI          d;
    WellOptParam opt;
};

/// TODO: Add Doxygen
class WellParam
{
public:
    /// Input well in structured grid
    WellParam(vector<string>& info);
    /// Input well in unstructured grid
    WellParam(vector<string>& info, const string& unstructured);
    WellParam(string name_, USI i, USI j, OCP_DBL depth_=-1.0); /// fff
//    WellParam(string name_, OCP_DBL x, OCP_DBL y, OCP_DBL z);
    /// Input perforations
    void InputCOMPDAT(vector<string>& vbuf);
    /// Input perforations in structured grid
    void InputCOMPDATS(vector<string>& vbuf);
    /// Input perforations in unstructured grid
    void InputCOMPDATUS(vector<string>& vbuf);
    /// Get number of perforations
    USI GetPerfNum() const { return max(I_perf.size(), X_perf.size()); }
    // static infomation
    // WELSPECS
    /// Grid type
    GridType gridType;
    // WELSPECS
    /// Name of Well
    string   name;
    /// Group the well belongs to.
    string   group{ "FEILD" };
    /// I index of well header, for structured grid
    USI      I;
    /// J index of well header, for structured grid
    USI      J;
    /// Depth of well header
    OCP_DBL  depth{ -1.0 };
    /// x-coordinate of well header, for unstructured grid
    OCP_DBL  X;
    /// y-coordinate of well header, for unstructured grid
    OCP_DBL  Y;
    /// z-coordinate of well header, for unstructured grid
    OCP_DBL  Z;

    // COMPDAT ---- for all perforation.
    /// I-index of perforation, for structured grid
    vector<USI>     I_perf;
    /// J-index of perforation, for structured grid
    vector<USI>     J_perf;
    /// K-index of perforation, for structured grid
    vector<USI>     K_perf;
    /// x-coordinate of perforation, for unstructured grid
    vector<OCP_DBL> X_perf;
    /// y-coordinate of perforation, for unstructured grid
    vector<OCP_DBL> Y_perf;
    /// z-coordinate of perforation, for unstructured grid
    vector<OCP_DBL> Z_perf;
    /// Transmissibility connection factor
    vector<OCP_DBL> WI;
    /// Diameter of perforations
    vector<OCP_DBL> diameter;
    /// Effective Kh
    vector<OCP_DBL> kh;
    /// Skin factor
    vector<OCP_DBL> skinFactor;
    /// Direction of perforations, for structured grid
    vector<string>  direction;
    /// If use unweighted injector.
    OCP_BOOL        ifUseUnweight{OCP_FALSE}; 
    /// initial BHP
    OCP_DBL initP{ -1 }; 
    // dynamic infomation
    vector<WellOptPair> optParam;
    void SetWellParams(USI i, USI j, USI k, OCP_DBL diam);

    static std::vector<std::string> Templates; /// Only for HiSim
};

/// Describe the molar fraction of components of fluid injected to reservoir from INJ.
class Solvent
{
public:
    Solvent() = default;
    Solvent(const vector<string>& vbuf);
    Solvent(const string& name_, const vector<double>& ratios);
    string          name;
    vector<OCP_DBL> comRatio;
};

/// ParamWell is an internal structure used to stores the information of wells from
/// input files. It is an intermediate interface and independent of the main simulator.
/// After all file inputting finishes, the params in it will pass to corresponding
/// modules.
class ParamWell
{
public:
    /// if thermal model is used
    OCP_BOOL          thermal{ OCP_FALSE };
    /// Contains all the information of wells.
    vector<WellParam> well;    
    /// Records the critical time given by users.
    vector<OCP_DBL>   criticalTime; 
    /// Sets of Solvent.
    vector<Solvent>   solSet;
    /// Pressure in surface condition.
    OCP_DBL           Psurf; 
    /// Temperature in surface condition.
    OCP_DBL           Tsurf;

    /// Construct a new well for HiSim, return the index of the newly constructed well
    struct WellOperation
    {
        double tstep;
        string name;
        string type; /// PROD, INJE
        WellOptParam opt;

        WellOperation(double ts, string name_, string type_, WellOptParam opt_)
                    : tstep(ts), name(name_), type(type_), opt(opt_) {}
    }; /// 从HiSim中读取的数据
    int ConstructNewWell(string name_, USI i, USI j);
    void SetSTREAM(const vector<double>& ratios);
    void AddWellOpts(int idx, const WellOptPair& opt);
    vector<WellOperation> well_oper_list;

    /// Initialize the inputting the params of wells.
    void Init();
    /// Initialize the critical time.
    void InitTime() { criticalTime.push_back(0); };
    /// Input the well keyword WELSPECS. WELSPECS defines wells including well name,
    /// well location, well depth and son on.
    void InputWELSPECS(ifstream& ifs);
    /// Input the well keyword COMPDAT. COMPDAT contains the information of perforations
    /// of wells, for example the location, the trans or directions.
    void InputCOMPDAT(ifstream& ifs);
    /// Input the well keyword WCONINJE. WCONINJE describes the initial operation mode
    /// of injection well.
    void InputWCONINJE(ifstream& ifs);
    /// Input the well keyword WCONPROD. WCONPROD describes the initial operation mode
    /// of production well.
    void InputWCONPROD(ifstream& ifs);
    /// Input the keyword: TSTEP. TSTEP is used to divide the simulation time and
    /// ususally the time point is critical, at which for example, the operation mode of
    /// well will change. So the params of solving equations could be adjusted
    /// correspondingly.
    void InputTSTEP(ifstream& ifs);
    /// Input the well keyword: WELTARG. WELTARG is used to change the operation mode of
    /// well anytime. For example, the oil production rate is changed from 1000 stb/day
    /// to 1500 stb/day at the 100th day.
    void InputWELTARG(ifstream& ifs);
    /// Input the temperature of injected fluid
    void InputWTEMP(ifstream& ifs);
    /// Input injector type -- MOBWEIGHT(defaulted) or UNWEIGHT
    void InputUNWEIGHT(ifstream& ifs);
    /// Input well keyword: Solvent. It describes the molar fraction of components of
    /// fluid injected to reservoir from INJ.
    void InputWELLSTRE(ifstream& ifs);
    /// Input surface pressure
    void InputPSURF(ifstream& ifs);
    /// Input surface temperature
    void InputTSURF(ifstream& ifs);
    /// Input initial BHP
    void InputWELINITP(ifstream& ifs);
    // check
    /// Check if wrong params are input.
    void CheckParam() const;
    /// Check if params of Perforation is wrong.
    void CheckPerf() const;
};

#endif /* end if __PARAMWELL_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/