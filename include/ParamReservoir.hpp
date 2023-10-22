/*! \file    ParamReservoir.hpp
 *  \brief   ParamReservoir class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARAMRESERVOIR_HEADER__
#define __PARAMRESERVOIR_HEADER__

// Standard header files
#include <fstream>
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "UtilInput.hpp"
#include "UtilOutput.hpp"

using namespace std;

/// TableSet is used to store a series of tables.
class TableSet
{
public:
    void DisplayTable() const; ///< Print table

public:
    /// Name of table.
    string                          name;  
    /// Number of columns of table.
    USI                             colNum; 
    /// All table with the same name.
    vector<vector<vector<OCP_DBL>>> data;    
};

/// 2D table, which means the table values usually depend on two variable, 
/// eg. pressure and temperature
class Table2
{
public:
    Table2(const USI& n) { data.resize(n); }
    void SetColNum() { colNum = data[0].size(); }
public:
    /// name of refData
    string                          refName;
    /// refData for each sub table
    vector<OCP_DBL>                 refData;
    /// Number of columns of table.
    USI                             colNum;  
    /// data
    vector<vector<vector<OCP_DBL>>> data;
};


/// a series of Table2
class Table2Set
{
public:
    /// Name of table
    string           name;    
    /// all table with the same name.
    vector<Table2>   data;    
};


class HLoss
{
public:
    OCP_BOOL ifHLoss{OCP_FALSE}; ///< If use Heat loss
    OCP_BOOL obUse{ OCP_FALSE }; ///< If use heat loss in overburden
    OCP_DBL  obC{ -1 };          ///< Volumetric heat capacity of overburden rock
    OCP_DBL  obK{ -1 };          ///< Thermal conductivity of overburden rock
    OCP_BOOL ubUse{ OCP_FALSE }; ///< If use heat loss in underburden
    OCP_DBL  ubC{ -1 };          ///< Volumetric heat capacity of underburden rock
    OCP_DBL  ubK{ -1 };          ///< Thermal conductivity of underburden rock
};

/// RockParam class contains information about the keyword ROCK.
class RockParam
{
public:
    string   type{"LINEAR"}; ///< LINEAR or EXPONENT for porosity model
    OCP_DBL  Pref{14.7};     ///< Reference pressure at initial porosity.
    OCP_DBL  Tref{60};       ///< Reference temperature at initial porosity.
    OCP_DBL  cp1{3.406E-6};  ///< Compressibility factor of rock in reservoir.
    OCP_DBL  cp2{0};         ///< 2 order Compressibility factor of rock in reservoir.
    OCP_DBL  ct{0};          ///< Expansion factor of rock in reservoir, ifThermal only
    OCP_DBL  cpt{0};         ///< cross items, ifThermal only
    OCP_BOOL ConstRock{OCP_TRUE}; ///< if true, rock volume remains const, else, bulk
                                  ///< volume remains const
    OCP_DBL HCP1{35}; ///< coefficients of the rock enthalpy formula, Btu/ft^3 - F
    OCP_DBL HCP2{0};  ///< coefficients of the rock enthalpy formula, Btu/ft^3 - F
};


/// Initial reservoir infomation for calculating initial equilibration.
class EQUILParam
{
public:
    vector<OCP_DBL> data;
};


/// Brooks-Corey type relative permeability and capillary pressure
class BrooksCoreyParam
{
public:
    /// Immobile wetting phase saturation
    OCP_DBL sw_imm;
    /// Immobile nonwetting phase(gas) saturation
    OCP_DBL sn_imm;
    /// Gas entry pressure
    OCP_DBL Pentry;
    /// Max capillary pressure
    OCP_DBL Pcmax;
    /// Shape exponents for relative permeability of wetting phase
    OCP_DBL Cw_kr;
    /// Shape exponents for relative permeability of nonwetting phase
    OCP_DBL Cn_kr;
    /// Shape exponents for capillary pressure
    OCP_DBL C_pc;
};


/// External boundary condition params
class BoundaryParam
{
public:
    BoundaryParam(const string& Name) : name(Name) {}

    /// boundary name
    string   name;
    /// if use const pressure
    OCP_BOOL constP{ OCP_FALSE };
    /// pressure
    OCP_DBL  P;
};


/// A internal structure used to store some params for reservoir, it can tell
/// if these params are given by users.
template <typename T>
class Type_A_r
{
public:
    OCP_BOOL  activity{OCP_FALSE}; ///< If OCP_FALSE, this param is not given.
    vector<T> data;                ///< Data of param.
};

/// ComponentParam contains information of components
class ComponentParam
{
public:
    // Basic params
    /// Init Params
    void Init();
    /// Input the information of hydrocarbon components
    void InputCOMPONENTS(ifstream& ifs, const string& keyword);
    /// Find corresponding variable according to the name of variable.
    /// It is used for the basic properties of hydrocarbon components such as TCRIT
    Type_A_r<vector<OCP_DBL>>* FindPtr01(const string& varName);
    /// input reference pressure, temperature
    void InputRefPR(ifstream& ifs, const string& keyword);
    /// Find corresponding variable according to the name of variable
    vector<OCP_DBL>* FindPtr02(const string& varName);
    /// Input the names of hydrocarbon components
    void InputCNAMES(ifstream& ifs);
    /// Input LBC coefficients for viscosity calculation
    void InputLBCCOEF(ifstream& ifs);
    /// Input the Binary interaction of components
    void InputBIC(ifstream& ifs);
    // Method params
    void InputSSMSTA(ifstream& ifs);
    void InputNRSTA(ifstream& ifs);
    void InputSSMSP(ifstream& ifs);
    void InputNRSP(ifstream& ifs);
    void InputRR(ifstream& ifs);

public:
    /// num of EOS region.
    USI            NTPVT;   
    /// num of components who will be in EOS calculations
    USI            numCom{0}; 
    /// num of phase who will be in EOS calculations
    USI            numPhase{2}; 
    /// Name of components
    vector<string> Cname; 
    /// Critical temperature of components
    Type_A_r<vector<OCP_DBL>> Tc;
    /// Critical pressure of components
    Type_A_r<vector<OCP_DBL>> Pc; 
    /// Critical volume of components
    Type_A_r<vector<OCP_DBL>> Vc;
    /// Critical Z-factor of components
    Type_A_r<vector<OCP_DBL>> Zc; 
    /// Molecular Weight of components
    Type_A_r<vector<OCP_DBL>> MW;
    /// Acentric factor of components
    Type_A_r<vector<OCP_DBL>> Acf; 
    /// OMEGA_A of components
    Type_A_r<vector<OCP_DBL>> OmegaA; 
    /// OMEGA_B of components
    Type_A_r<vector<OCP_DBL>> OmegaB; 
    /// Volume shift of components
    Type_A_r<vector<OCP_DBL>> Vshift;  
    /// PARACHOR of components
    Type_A_r<vector<OCP_DBL>> parachor; 

    // for viscosity calculation
    ///< Critical volume used for viscosity calculations only.
    Type_A_r<vector<OCP_DBL>> Vcvis; 
    /// Critical Z-factor used for viscosity calculations only.
    Type_A_r<vector<OCP_DBL>> Zcvis; 
    /// LBC coefficients for viscosity calculation
    vector<OCP_DBL>           LBCcoef;     
    /// Binary interaction
    vector<vector<OCP_DBL>>   BIC;

    // Thermal only
    /// component molar density at reference temperature and reference pressure, lb/ft3
    Type_A_r<vector<OCP_DBL>> molden; 
    /// component compressibility, 1/psi
    Type_A_r<vector<OCP_DBL>> cp;     
    /// the first ifThermal expansion coefficient, 1/F
    Type_A_r<vector<OCP_DBL>> ct1;  
    /// the second ifThermal expansion coefficient, 1/F
    Type_A_r<vector<OCP_DBL>> ct2;  
    /// the coefficient of density dependence on temperature and pressure, 1/psi-F
    Type_A_r<vector<OCP_DBL>> cpt;  
    /// coefficients in the component liquid enthalpy calculations, Btu/lbmol/F
    Type_A_r<vector<OCP_DBL>> cpl1; 
    /// coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^2
    Type_A_r<vector<OCP_DBL>> cpl2; 
    /// coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^3
    Type_A_r<vector<OCP_DBL>> cpl3; 
    /// coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^4
    Type_A_r<vector<OCP_DBL>> cpl4; 
    /// coefficients in the component liquid enthalpy calculations, Btu/lbmol/F
    Type_A_r<vector<OCP_DBL>> cpg1; 
    /// coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^2
    Type_A_r<vector<OCP_DBL>> cpg2; 
    /// coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^3
    Type_A_r<vector<OCP_DBL>> cpg3; 
    /// coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^4
    Type_A_r<vector<OCP_DBL>> cpg4; 
    /// coefficients in the component gas enthalpy calculations, Btu/lbmol
    Type_A_r<vector<OCP_DBL>> hvapr; 
    /// coefficients in the vaporization enthalpy calculations
    Type_A_r<vector<OCP_DBL>> hvr; 
    /// coefficients in the vaporization enthalpy calculations
    Type_A_r<vector<OCP_DBL>> ev; 
    /// coefficients in water and oil viscosity correlation formulae
    Type_A_r<vector<OCP_DBL>> avisc;
    /// coefficients in water and oil viscosity correlation formulae
    Type_A_r<vector<OCP_DBL>> bvisc;
    /// coefficients Ak in gas viscosity correlation formulae
    Type_A_r<vector<OCP_DBL>> avg;
    /// coefficients Bk in gas viscosity correlation formulae
    Type_A_r<vector<OCP_DBL>> bvg; 

    /// viscosity-versus-temperature dependence, This table can specify the
    /// viscosity-versus-temperature-pressure dependence.
    Table2Set viscTab;

    /// reference pressure
    vector<OCP_DBL> Pref; 
    /// reference temperature
    vector<OCP_DBL> Tref;
    /// Params for Solving Phase Spliting with SSM
    vector<string> SSMparamSTA; 
    /// Params for Solving Phase Spliting with NR
    vector<string> NRparamSTA;  
    /// Params for Solving Phase Spliting with SSM
    vector<string> SSMparamSP;  
    /// Params for Solving Phase Spliting with NR
    vector<string> NRparamSP; 
    /// Params for Solving Rachford-Rice equations
    vector<string> RRparam;     
};


class Miscstr
{
public:
    // vector<OCP_DBL> surTenRef;
    OCP_BOOL        ifMiscible{ OCP_FALSE };
    OCP_DBL         surTenRef{ -1 };  ///< reference surface tension - flow is immiscible when the surface tension is greater than or equal to this value.
    OCP_DBL         surTenEpt{ -1 };  ///< maximum surface tension expected, it should be greater than surTenRef.
    OCP_DBL         surTenPc{ -1 };   ///< maximum surface tension used to scale the input capillary pressure curves.
    OCP_DBL         surTenExp{ 0.25 };  ///< exponent of the surface tension ratio
};


/// ParamReservoir is an internal structure used to stores the information of
/// reservoir(except wells) from input files. It is an intermediate interface and
/// independent of the main simulator. After all file inputting finishs, the params in
/// it will pass to corresponding modules.
class ParamReservoir
{

public:

    string                   unitType;

    /// Temperature for reservoir.
    OCP_DBL                  rsTemp; 
    /// a set of rock params
    vector<RockParam>        rockSet;
    /// Heat loss property
    HLoss                    hLoss;  
    /// reference Miscibility surface tension
    Miscstr                  miscstr; 
    /// params for Brooks-Corey model
    vector<BrooksCoreyParam> BCparam;
    /// params for boundary condition
    vector<BoundaryParam>    BDparam;

    // phase property
    Type_A_r<OCP_DBL> density; ///< Density of oil, water, gas in standard conditions.
    Type_A_r<OCP_DBL> gravity; ///< Gravity of oil, water, gas in standard conditions.
    OCP_BOOL          ifThcon{ OCP_FALSE };
    OCP_DBL           thcono{24}; ///< oil ifThermal conductivity
    OCP_DBL           thcong{24}; ///< gas ifThermal conductivity
    OCP_DBL           thconw{24}; ///< water ifThermal conductivity
    OCP_DBL           thconr{24}; ///< Rock ifThermal conductivity.

    // mixture Models
    OCP_BOOL blackOil{OCP_FALSE}; ///< If ture, blackoil model will be used.
    OCP_BOOL comps{OCP_FALSE};    ///< If OCP_TRUE, compositional model will be used.
    OCP_BOOL thermal{OCP_FALSE};  ///< If OCP_TRUE, ifThermal model will be used.
    OCP_BOOL oil{OCP_FALSE};      ///< If OCP_TRUE, oil phase could exist.
    OCP_BOOL gas{OCP_FALSE};      ///< If OCP_TRUE, gas phase could exist.
    OCP_BOOL water{OCP_FALSE};    ///< If OCP_TRUE, water phase could exist.
    OCP_BOOL disGas{OCP_FALSE};   ///< If OCP_TRUE, dissolve gas could exist in oil phase.

    // flow model
    /// Use gravity drainage for dual porosity runs
    OCP_BOOL GRAVDR{ OCP_FALSE };

    ComponentParam comsParam; ///< information for components

    // SAT Region & PVT Region
    USI               NTSFUN{1}; ///< Num of SAT regions.
    USI               NTPVT{1};  ///< Num of PVT regions.
    USI               NTROOC{1}; ///< Num of Rock regions.

    // Saturation tables & bubble point pressure
    TableSet SWFN_T; ///< Table set of SWFN.
    TableSet SWOF_T; ///< Table set of SWOF.
    TableSet SGFN_T; ///< Table set of SGFN.
    TableSet SGOF_T; ///< Table set of SGOF.
    TableSet SOF3_T; ///< Table set of SOF3.
    TableSet PBVD_T; ///< Table set of PBVD.
    // initial zi vs depth
    TableSet           ZMFVD_T;  ///< Table set of ZMFVD
    TableSet           TEMPVD_T; ///< Table set of TEMPVD
    vector<EQUILParam> EQUIL;    ///< See ParamEQUIL.

    // PVT properties
    USI numPhase; ///< Number of phases
    USI numCom; ///< Number of components(hydrocarbon components), used in Compositional
                ///< Model when input
    TableSet PVCO_T; ///< Table set of PVCO.
    TableSet PVDO_T; ///< Table set of PVDO.
    TableSet PVCDO_T; ///< Table set of PVCDO.
    TableSet PVDG_T; ///< Table set of PVDG.
    TableSet PVTW_T; ///< Table set of PVTW.

    /// Use Garcia's method to calculate the water density
    OCP_BOOL  GARCIAW{ OCP_FALSE };
    /// PVT property for H2O
    Table2Set PVTH2O;
    /// PVT property for CO2
    Table2Set PVTCO2;
    /// Pressure in surface condition.
    OCP_DBL   Psurf;
    /// Temperature in surface condition.
    OCP_DBL   Tsurf;

public:

    /// Find corresponding variable according to the name of variable.
    /// It is used for the scope of the table.
    TableSet* FindPtrTable(const string& varName);
    Table2Set* FindPtrTable2(const string& varName);

    /// Initialize the default value in reservoir, such as temperature, density, table.
    void Init();

    /// Initialize the tables' name and num of colum.
    void InitTable();

    /// Input the keyword: COMPS. COMPS is used in compositional model, which gives the
    /// num of components.
    void InputCOMPS(ifstream& ifs);

    /// Input the keyword: RTEMP. RTEMP gives the temperature of reservoir.
    void InputRTEMP(ifstream& ifs);

    /// Input table, fro example, PVTtable and SATtable such as SWOF, PVCO.
    void InputTABLE(ifstream& ifs, const string& tabName);
    /// Input Table2, for example, VISCTAB
    void InputTABLE2(ifstream& ifs, const string& tabName);

    /// Input the keyword: ROCK. ROCK contains the compressibility factor and reference
    /// pressure at initial porosity.
    void InputROCK(ifstream& ifs);
    /// Input Rock information for ifThermal model
    void InputROCKT(ifstream& ifs);
    /// Input heat loss property for overburden rock and underburden rock
    void InputHLOSS(ifstream& ifs);
    /// Input params for Brooks-Corey model
    void InputBrooksCorey(ifstream& ifs);

    /// Input the Miscibility information
    void InputMISCSTR(ifstream& ifs);

    /// Input the reference gravity of oil, water, and air in standard condition.
    void InputGRAVITY(ifstream& ifs);

    /// Input the reference density of oil, water, and air in standard condition.
    void InputDENSITY(ifstream& ifs);

    /// Input the phase ifThermal conductivity
    void InputTHCON(ifstream& ifs, const string& keyword);

    /// EQUIL contains initial information of reservoir; see ParamEQUIL.
    void InputEQUIL(ifstream& ifs);

    // SATNUM & PVTNUM  -- Region
    /// TABDIMS contains the num of saturation region and PVT region.
    void InputTABDIMS(ifstream& ifs);

    // Input ComponentParam
    // Basic params
    void InputNCOMPS(ifstream& ifs) { 
        vector<string> vbuf; 
        ReadLine(ifs, vbuf);
        comsParam.numCom = stoi(vbuf[0]);
        numCom           = comsParam.numCom;
    }
    void InputCNAMES(ifstream& ifs) { comsParam.InputCNAMES(ifs); };
    void InputCOMPONENTS(ifstream& ifs, const string& keyword)
    {
        comsParam.InputCOMPONENTS(ifs, keyword);
    }
    void InputLBCCOEF(ifstream& ifs) { comsParam.InputLBCCOEF(ifs); }
    void InputBIC(ifstream& ifs) { comsParam.InputBIC(ifs); };
    void InputRefPR(ifstream& ifs, const string& keyword)
    {
        comsParam.InputRefPR(ifs, keyword);
    };

    // PEC Method params
    void InputSSMSTA(ifstream& ifs) { comsParam.InputSSMSTA(ifs); };
    void InputNRSTA(ifstream& ifs) { comsParam.InputNRSTA(ifs); };
    void InputSSMSP(ifstream& ifs) { comsParam.InputSSMSP(ifs); };
    void InputNRSP(ifstream& ifs) { comsParam.InputNRSP(ifs); };
    void InputRR(ifstream& ifs) { comsParam.InputRR(ifs); };

    /// Input boundary conditons
    void InputBoundary(ifstream& ifs);

    // check
    /// Check the reservoir param from input file.
    void CheckParam();

    /// Check Rock
    void CheckRock();

    /// Check cpl1, cpl2, cpl3, cpl4
    void CheckCPL();

    /// Check cpg1, cpg2, cpg3, cpg4
    void CheckCPG();
};

#endif /* end if __PARAMRESERVOIR_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/09/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/