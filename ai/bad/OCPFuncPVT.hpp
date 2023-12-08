 /*! \file    OCPFuncPVT.hpp
*   \brief   Functions for PVT in OCP
*   \author  Shizhe Li
*   \date    Jun/18/2023
*
*-----------------------------------------------------------------------------------
*  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
*  Released under the terms of the GNU Lesser General Public License 3.0 or later.  
*-----------------------------------------------------------------------------------
*/

#ifndef __OCPFUNCPVT_HEADER__ 
#define __OCPFUNCPVT_HEADER__  

// OpenCAEPoroX header files
#include "OCPFuncTable.hpp"  
#include "UtilMath.hpp"
using namespace std;

/** @defgroup PVTFunctions PVT Functions
* @brief PVT functions for fluid properties
* @details This file contains various PVT functions to calculate fluid properties like density, viscosity etc.
* @{
*/

/** @brief Water PVT functions 
* @details Calculates water properties like density, viscosity etc. as a function of pressure
*/
class OCP_PVTW : public OCPFuncTable
{
    public:
        //! Constructor
        OCP_PVTW() = default; 

        //! Set up the table
        void Setup(const vector<vector<OCP_DBL>>& src, const OCP_DBL& stdRhoWin, const OCP_DBL& stdVwin);
        
        //! Calculate water compressibility
        OCP_DBL CalXiW(const OCP_DBL& P) const;  

        //! Calculate water density 
        OCP_DBL CalRhoW(const OCP_DBL& P) const;

        //! Calculate water density and derivatives
        void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu, 
                           OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const;

    protected:
        //! Calculate water formation volume factor
        OCP_DBL CalBw(const OCP_DBL& P) const;

        //! Calculate water formation volume factor and viscosity derivatives
        void CalBwMuwDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu, OCP_DBL& bP, OCP_DBL& muP) const;

    protected:
        //! Mass density of water phase in standard condition 
        OCP_DBL stdRhoW;
        
        //! Molar volume of water phase in standard conditions
        OCP_DBL stdVw;  
};

/** @brief Live oil PVT properties  
* @details Calculates live oil properties like density, viscosity etc. as a function of pressure
*/
class OCP_PVCO : public OCPFuncTable  
{
    public:
        //! Constructor 
        OCP_PVCO() = default;

        //! Set up the table
        void Setup(const vector<vector<OCP_DBL>>& src, const OCP_DBL& stdRhoOin, 
                   const OCP_DBL& stdRhoGin, const OCP_DBL& stdVoin, const OCP_DBL& stdVgin);
        
        //! Calculate oil density
        OCP_DBL CalRhoO(const OCP_DBL& P, const OCP_DBL& Pb) const;

        //! Calculate oil compressibility
        OCP_DBL CalXiO(const OCP_DBL& P, const OCP_DBL& Pb) const;

        //! Calculate density, viscosity etc. derivatives for saturated oil 
        void CalRhoXiMuRsDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu, 
                             OCP_DBL& rs, OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP, OCP_DBL& rsP) const;
                             
        //! Calculate density, viscosity etc. derivatives for undersaturated oil
        void CalRhoXiMuDer(const OCP_DBL& rs, const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, 
                           OCP_DBL& mu, OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP, 
                           OCP_DBL& rhoRs, OCP_DBL& xiRs, OCP_DBL& muRs) const;
                           
        //! Calculate dissolved gas-oil ratio for saturated oil
        OCP_DBL CalRs(const OCP_DBL& P) const;

    protected:
        //! Calculate derivatives for saturated oil  
        void CalRsBoMuoDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& rs, OCP_DBL& mu,
                           OCP_DBL& bP, OCP_DBL& rsP, OCP_DBL& muP) const;
                           
        //! Calculate derivatives for undersaturated oil
        void CalBoMuoDer(const OCP_DBL& rs, const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu,  
                         OCP_DBL& bP, OCP_DBL& muP, OCP_DBL& bRs, OCP_DBL& muRs) const;
                         
    protected:
        //! Mass density of oil phase in standard condition
        OCP_DBL stdRhoO;   
        
        //! Mass density of gas phase in standard condition 
        OCP_DBL stdRhoG;
        
        //! Molar volume of oil in standard condition 
        OCP_DBL stdVo;
        
        //! Molar volume of gas in standard condition
        OCP_DBL stdVg;              
};

/** @brief Dry gas PVT properties
* @details Calculates dry gas properties like density, viscosity etc. as a function of pressure  
*/
class OCP_PVDG : public OCPFuncTable
{
    public:
        //! Constructor
        OCP_PVDG() = default;

        //! Set up the table  
        void Setup(const vector<vector<OCP_DBL>>& src, 
                   const OCP_DBL& stdRhoGin, const OCP_DBL& stdVgin);
                   
        //! Calculate gas compressibility
        OCP_DBL CalXiG(const OCP_DBL& P) const;

        //! Calculate gas density
        OCP_DBL CalRhoG(const OCP_DBL& P) const;

        //! Calculate gas density and derivatives 
        void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu,
                           OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const;
                           
    protected:
        //! Calculate gas formation volume factor  
        OCP_DBL CalBg(const OCP_DBL& P) const;

        //! Calculate gas formation volume factor and viscosity derivatives 
        void CalBgMugDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu, 
                         OCP_DBL& bP, OCP_DBL& muP) const;
                         
    protected:
        //! Mass density of gas phase in standard condition
        OCP_DBL stdRhoG;
        
        //! Molar volume of gas phase in standard conditions 
        OCP_DBL stdVg;
};

/** @brief Dead oil PVT properties
* @details Calculates dead oil properties like density, viscosity etc. as a function of pressure
*/  
class OCP_PVDO : public OCPFuncTable
{
    public:
        //! Constructor
        OCP_PVDO() = default;

        //! Set up the table 
        virtual void Setup(const vector<vector<OCP_DBL>>& src, 
                           const OCP_DBL& stdRhoOin, const OCP_DBL& stdVoin);
                           
        //! Calculate oil compressibility
        virtual OCP_DBL CalXiO(const OCP_DBL& P) const;

        //! Calculate oil density
        virtual OCP_DBL CalRhoO(const OCP_DBL& P) const;

        //! Calculate oil density and derivatives
        virtual void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, 
                                   OCP_DBL& mu, OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const;
                                   
    protected:
        //! Calculate oil formation volume factor  
        virtual OCP_DBL CalBo(const OCP_DBL& P) const;

        //! Calculate oil formation volume factor and viscosity derivatives
        virtual void CalBoMuoDer(const OCP_DBL& P, OCP_DBL& bo, OCP_DBL& muo, 
                                 OCP_DBL& dBodP, OCP_DBL& dMuodP) const;
                                 
    protected:
        //! Mass density of oil phase in standard conditions
        OCP_DBL stdRhoO;    
        
        //! Molar volume of oil phase in standard conditions
        OCP_DBL stdVo;     
};

/** @brief Dead oil PVT properties with constant compressibility 
* @details Calculates dead oil properties like density, viscosity etc. with constant compressibility
*/
class OCP_PVCDO : public OCP_PVDO
{
    public:
        //! Constructor
        OCP_PVCDO() = default;

        //! Set up the table
        void Setup(const vector<vector<OCP_DBL>>& src, 
                   const OCP_DBL& stdRhoOin, const OCP_DBL& stdVoin) override;
                   
        //! Calculate oil compressibility  
        OCP_DBL CalXiO(const OCP_DBL& P) const override;

        //! Calculate oil density
        OCP_DBL CalRhoO(const OCP_DBL& P) const override;

        //! Calculate oil density and derivatives
        void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, 
                           OCP_DBL& mu, OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const override;
                           
    protected:
        //! Calculate oil formation volume factor
        OCP_DBL CalBo(const OCP_DBL& P) const;

        //! Calculate oil formation volume factor and viscosity derivatives  
        void CalBoMuoDer(const OCP_DBL& P, OCP_DBL& bo, OCP_DBL& muo, 
                         OCP_DBL& dBodP, OCP_DBL& dMuodP) const;
                         
    protected:
        //! Reference pressure
        OCP_DBL Pref;
        
        //! Oil formation volume factor at the reference pressure
        OCP_DBL Bref;
        
        //! The oil compressibility 
        OCP_DBL Cb;
        
        //! The oil viscosity at the reference pressure
        OCP_DBL muref;
        
        //! The oil "viscosibility"  
        OCP_DBL Cmu;
};

/** @brief General PVT table with temperature dependence  
* @details Calculates fluid properties like density, viscosity etc. as a function of pressure and temperature
*/
class OCP_PVT2 : public OCPFuncTable2
{
    public:
        //! Calculate density
        OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& T) const;

        //! Calculate density, viscosity and corresponding solubility 
        void CalRhoMuSol(const OCP_DBL& P, const OCP_DBL& T, 
                         OCP_DBL& rho, OCP_DBL& mu, OCP_DBL& sol) const;
                         
        //! Calculate density, viscosity, solubility and derivatives
        void CalRhoMuSolDer(const OCP_DBL& P, const OCP_DBL& T, 
                            OCP_DBL& rho, OCP_DBL& mu, OCP_DBL& sol,
                            OCP_DBL& rhoP, OCP_DBL& muP, OCP_DBL& solP) const;
};

// Typedefs for specialized PVT tables
typedef OCP_PVT2 OCP_PVTCO2; 
typedef OCP_PVT2 OCP_PVTH2O;

/** @brief Garciaw model for water density with CO2 dissolution  
* @details Calculates water density considering CO2 dissolution
*/
class Garciaw
{
    public:
        //! Set up 
        void Setup(const OCP_BOOL& flag);

        //! Check if model is used
        auto IfUse() const;

        //! Calculate water density
        void CalRho(const OCP_DBL& T, const OCP_DBL& xGw, OCP_DBL& rhow) const;

        //! Calculate water density and derivatives  
        void CalRhoDer(const OCP_DBL& T, const OCP_DBL& xGw, const OCP_DBL& xGwP, 
                       OCP_DBL& rhow, OCP_DBL& rhowP, OCP_DBL& drhow_dxGw) const;
                       
        //! Calculate water density and derivatives
        void CalRhoDer(const OCP_DBL& T, const OCP_DBL& xGw, const OCP_DBL& xGwP,  
                       const OCP_DBL& xGwT, OCP_DBL& rhow, OCP_DBL& rhowP, 
                       OCP_DBL& rhowT, OCP_DBL& drhow_dxGw) const;
                       
    protected:
        //! If use this model
        OCP_BOOL ifUse;
        
        //! Molecular weight of CO2
        const OCP_DBL MWCO2; 
};

/** @brief Viscosity calculation parameters  
*/
class ViscosityParams 
{
    public:
        //! Constructor
        ViscosityParams(const OCP_DBL* Pin, const OCP_DBL* Tin, const OCP_DBL* xin);

        //! Constructor  
        ViscosityParams(const OCP_DBL* Pin, const OCP_DBL* Tin, 
                        const OCP_DBL* xin, const OCP_DBL* xiin);
                        
        //! Constructor
        ViscosityParams(const OCP_DBL* Pin, const OCP_DBL* Tin,
                        const OCP_DBL* xin, const OCP_DBL* xiin,  
                        const OCP_DBL* xiPin, const OCP_DBL* xiTin,
                        const OCP_DBL* xixin);
                        
    public:
        //! Pressure
        const OCP_DBL* P;
        
        //! Temperature 
        const OCP_DBL* T;
        
        //! Molar fraction
        const OCP_DBL* x;
        
        //! Molar density 
        const OCP_DBL* xi;
        
        //! Derivative of molar density w.r.t. pressure
        const OCP_DBL* xiP;
        
        //! Derivative of molar density w.r.t. temperature 
        const OCP_DBL* xiT;
        
        //! Derivative of molar density w.r.t. composition
        const OCP_DBL* xix;
};

/** @brief Base viscosity calculation class  
*/
class ViscosityMethod 
{
    public:
        //! Constructor
        ViscosityMethod() = default;

        //! Calculate viscosity  
        virtual OCP_DBL CalViscosity(const ViscosityParams& vp) = 0;

        //! Calculate viscosity and derivatives 
        virtual OCP_DBL CalViscosity(const ViscosityParams& vp, 
                                     OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) = 0;
                                     
    protected:
        //! Number of components
        USI nc;
        
        //! Viscosity of components  
        vector<OCP_DBL> muc;
        
        //! Derivative of viscosity w.r.t. pressure
        vector<OCP_DBL> mucP;
        
        //! Derivative of viscosity w.r.t. temperature
        vector<OCP_DBL> mucT;
};

/** @brief Viscosity calculation using table and linear mixing  
*/
class ViscosityMethod01 : public ViscosityMethod
{
    public:
        //! Constructor 
        ViscosityMethod01(const Table2& tab);

        //! Calculate viscosity
        OCP_DBL CalViscosity(const ViscosityParams& vp) override;

        //! Calculate viscosity and derivatives
        OCP_DBL CalViscosity(const ViscosityParams& vp,  
                             OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) override;
                             
    protected:
        //! Viscosity table (function of temperature and pressure) 
        OCPTable2 viscTab;
};

/** @brief Viscosity calculation using correlation parameters and linear mixing
*/ 
class ViscosityMethod02 : public ViscosityMethod
{
    public:
        //! Constructor
        ViscosityMethod02(const vector<OCP_DBL>& av, const vector<OCP_DBL>& bv);

        //! Calculate viscosity 
        OCP_DBL CalViscosity(const ViscosityParams& vp) override;

        //! Calculate viscosity and derivatives
        OCP_DBL CalViscosity(const ViscosityParams& vp,  
                             OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) override;
                             
    protected:
        //! Coefficients for viscosity correlation  
        vector<OCP_DBL> avisc;
        
        //! Coefficients for viscosity correlation
        vector<OCP_DBL> bvisc;
};

/** @brief Lohrenz-Bray-Clark viscosity calculation  
*/
class ViscosityMethod03 : public ViscosityMethod 
{
    public:
        //! Constructor
        ViscosityMethod03(const ComponentParam& param, const USI& tarId);

        //! Calculate viscosity
        OCP_DBL CalViscosity(const ViscosityParams& vp) override;

        //! Calculate viscosity and derivatives 
        OCP_DBL CalViscosity(const ViscosityParams& vp,  
                             OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) override;
                             
    protected:
        //! Number of components
        USI nc;
        
        //! LBC coefficients

