/*! \file    MixtureComp.hpp
 *  \brief   MixtureComp class declaration
 *  \author  Shizhe Li
 *  \date    Nov/30/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __MIXTURECOMP_HEADER__
#define __MIXTURECOMP_HEADER__

// Standard header files
#include <algorithm>
#include <math.h>
#include <vector>

// OpenCAEPoroX header files
#include "DenseMat.hpp"
#include "UtilMath.hpp"
#include "MixtureUnit.hpp"
#include "OCPFuncPVT.hpp"

using namespace std;

/// Params for SSM in Phase Stability Analysis
class SSMparamSTA
{
public:
    USI      maxIt;      ///< Max Iteration
    OCP_DBL  tol;        ///< Tolerance
    OCP_DBL  Ktol{1E-4}; ///< tolerace^2 for K
    OCP_DBL  dYtol{1E-6};
    OCP_DBL  eYt{1E-8}; ///< if Yt > 1 + eYt, then single phase is unstable
    OCP_DBL  tol2;      ///< tol*tol
    OCP_DBL  realTol;   ///< Real tol
    OCP_BOOL conflag;   ///< convergence flag, if converges, conflag = OCP_TRUE
                        // test
    USI     curIt;      ///< current Iterations
    OCP_DBL curSk;
};

/// Params for NR in Phase Stability Analysis
class NRparamSTA
{
public:
    USI      maxIt;   ///< Max Iteration
    OCP_DBL  tol;     ///< Tolerance
    OCP_DBL  tol2;    ///< tol*tol
    OCP_DBL  realTol; ///< Real tol
    OCP_BOOL conflag; ///< convergence flag, if converges, conflag = OCP_TRUE
    // test
    USI curIt; ///< current Iters
};

/// Params for SSM in Phase Split
class SSMparamSP
{
public:
    USI      maxIt;   ///< Max Iteration
    OCP_DBL  tol;     ///< Tolerance
    OCP_DBL  tol2;    ///< tol*tol
    OCP_DBL  realTol; ///< Real tol
    OCP_BOOL conflag; ///< convergence flag, if converges, conflag = OCP_TRUE
    // test
    USI curIt; ///< current Iters
};

/// Params for NR in Phase Split
class NRparamSP
{
public:
    USI      maxIt;   ///< Max Iteration
    OCP_DBL  tol;     ///< Tolerance
    OCP_DBL  tol2;    ///< tol*tol
    OCP_DBL  realTol; ///< Real tol
    OCP_BOOL conflag; ///< convergence flag, if converges, conflag = OCP_TRUE
    // test
    USI curIt; ///< current Iters
};

/// Param for Solving Rachford-Rice Equations
class RRparam
{
public:
    USI     maxIt; ///< Max Iteration
    OCP_DBL tol;   ///< Tolerance
    OCP_DBL tol2;  ///< tol*tol
    // test
    USI curIt; ///< current Iters
};

class EoScontrol
{
    friend class MixtureComp;

private:
    SSMparamSTA SSMsta;
    NRparamSTA  NRsta;
    SSMparamSP  SSMsp;
    NRparamSP   NRsp;
    RRparam     RR;
};


class MixtureComp : public MixtureUnit
{

public:
    OCP_DBL GetErrorPEC() override { return ePEC; }
    void    OutMixtureIters() const override;

private:
    // total iters
    OCP_ULL itersSSMSTA{0};
    OCP_ULL itersNRSTA{0};
    OCP_ULL itersSSMSP{0};
    OCP_ULL itersNRSP{0};
    OCP_ULL itersRR{0};
    // total counts, one count may contain many iters
    OCP_ULL countsSSMSTA{0};
    OCP_ULL countsNRSTA{0};
    OCP_ULL countsSSMSP{0};
    OCP_ULL countsNRSP{0};

    OCP_ULL countsFailed{0};
    // phase equilibrium calculation error
    // if NP = 1, it's from phase stable analysis, if skiped, it's 0
    // if NP > 1, it's from phase spliting calculation
    OCP_DBL ePEC;

public:
    MixtureComp() = default;

    MixtureComp(const ParamReservoir& rs_param, const USI& i)
        : MixtureComp(rs_param.comsParam, i)
    {
        // water property
        if (rs_param.PVTW_T.data.size() != 0) {
            if (rs_param.gravity.activity)
                std_RhoW = RHOW_STD * rs_param.gravity.data[1];
            if (rs_param.density.activity) std_RhoW = RHOW_STD;
            PVTW.Setup(rs_param.PVTW_T.data[i], std_RhoW);
        }

        // Input Miscible Params
        InputMiscibleParam(rs_param, i);
    };
    /// Input and Setup basic component params
    MixtureComp(const ComponentParam& param, const USI& i);

    void InitPTZ(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Ziin)
    {
        P = Pin;
        T = Tin;
        Dcopy(NC, &zi[0], Ziin);
    }
    void InitPTN(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
    {
        P = Pin;
        T = Tin;
        Dcopy(numCom, &Ni[0], Niin);
        Nh = Dnorm1(NC, &Ni[0]);
        for (USI i = 0; i < NC; i++) zi[i] = Ni[i] / Nh;
    }

    void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) override;

    void InitFlashIMPEC(const OCP_DBL& Pin,
                        const OCP_DBL& Pbbin,
                        const OCP_DBL& Tin,
                        const OCP_DBL* Sjin,
                        const OCP_DBL& Vpore,
                        const OCP_DBL* Ziin,
                        const OCP_USI& bId) override;

    void InitFlashFIM(const OCP_DBL& Pin,
                      const OCP_DBL& Pbbin,
                      const OCP_DBL& Tin,
                      const OCP_DBL* Sjin,
                      const OCP_DBL& Vpore,
                      const OCP_DBL* Ziin,
                      const OCP_USI& bId) override;

    // ftype = 0, flash from single phase
    // ftype = 1, skip phase stability analysis and num of phase = 1
    // ftype = 1, skip phase stability analysis and num of phase = 2
    void FlashIMPEC(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const OCP_DBL* Niin,
                    const USI&     lastNP,
                    const OCP_DBL* xijin,
                    const OCP_USI& bId) override;

    void CalFlash();

    void FlashFIM(const OCP_DBL& Pin,
                  const OCP_DBL& Tin,
                  const OCP_DBL* Niin,
                  const OCP_DBL* Sjin,
                  const USI&     lastNP,
                  const OCP_DBL* xijin,
                  const OCP_USI& bId) override;


    OCP_DBL
    XiPhase(const OCP_DBL& Pin,
            const OCP_DBL& Tin,
            const vector<OCP_DBL>& Ziin,
            const USI&     tarPhase) override;

    OCP_DBL
    RhoPhase(const OCP_DBL& Pin,
             const OCP_DBL& Pbb,
             const OCP_DBL& Tin,
             const vector<OCP_DBL>& Ziin,
             const USI&     tarPhase) override;

    // For Well
    void SetupWellOpt(WellOpt&                  opt,
                      const vector<SolventINJ>& sols,
                      const OCP_DBL&            Psurf,
                      const OCP_DBL&            Tsurf) override;
    void CalProdWeight(const OCP_DBL&         Pin,
                       const OCP_DBL&         Tin,
                       const OCP_DBL*         Niin,
                       const vector<OCP_DBL>& prodPhase,
                       vector<OCP_DBL>&       prodWeight) override;

    void CalProdRate(const OCP_DBL&   Pin,
                     const OCP_DBL&   Tin,
                     const OCP_DBL*   Niin,
                     vector<OCP_DBL>& prodRate) override;

    OCP_DBL CalInjWellEnthalpy(const OCP_DBL& Tin, const OCP_DBL* Ziin) override
    {
        OCP_ABORT("Can not be used in Compositional Model!");
    }


protected:
    // Basic Components Informations
    /// num of hydrocarbon components
    USI             NC;   
    /// name of hydrocarbon components
    vector<string>  Cname;  
    /// critical temperature of hydrocarbon components
    vector<OCP_DBL> Tc; 
    /// critical pressure of hydrocarbon components
    vector<OCP_DBL> Pc;
    /// critical volume of hydrocarbon components
    vector<OCP_DBL> Vc; 
    /// molecular Weight of hydrocarbon components
    vector<OCP_DBL> MWC; 
    /// acentric factor of hydrocarbon components
    vector<OCP_DBL> Acf;
    // for viscosity calculation
    /// critical volume used for viscosity calculations only
    vector<OCP_DBL> Vcvis;
    /// critical Z-factor used for viscosity calculations only
    vector<OCP_DBL> Zcvis;
    /// LBC coefficients for viscosity calculation
    vector<OCP_DBL> LBCcoef;

protected:
    // water property
    /// PVTW
    OCP_PVTW        PVTW;
    /// mass density of water phase in standard condition.
    OCP_DBL         std_RhoW;

protected:
    // EoS Calculations
    EoSCalculation  eos;

protected:
    // Flash method Control
    /// allowable maximum num of hydrocarbon phases
    USI             NPmax;
    /// method params for solving phase equilibrium
    EoScontrol      EoSctrl;

protected:
    // Initial properties for flash   
    /// current Pressure
    OCP_DBL         P;
    /// current Temperature
    OCP_DBL         T;
    /// mole fraction of hydrocarbon components
    vector<OCP_DBL> zi;

protected:
    // The following variables are all for hydrocarbon phases and components
    /// last num of phase in last NR Step from external iterations
    USI                      lNP; 
    /// current num of phase
    USI                      NP;    
    /// total moles of components
    OCP_DBL                  Nh;  
    /// volume of phase
    vector<OCP_DBL>          vC; 
    /// molar fraction of phase
    vector<OCP_DBL>          nu; 
    /// molar fraction of i-th component in jth phase
    vector<vector<OCP_DBL>>  x;
    /// fugacity coefficient of i-th component in j-th phase
    vector<vector<OCP_DBL>>  phi; 
    /// fugacity of i-th component in j-th phase
    vector<vector<OCP_DBL>>  fug;
    /// moles of i-th component in jth phase
    vector<vector<OCP_DBL>>  n;
    /// last n in NR iterations in phase spliting calculations
    vector<vector<OCP_DBL>>  ln;  
    /// Gibbs energy, before flash (not true value)
    OCP_DBL                  GibbsEnergyB; 
    /// Gibbs energy, after flash (not true value)
    OCP_DBL                  GibbsEnergyE; 
    /// molar density of phase
    vector<OCP_DBL>          xiC;
    /// mass density of phase
    vector<OCP_DBL>          rhoC;  
    /// molecular weight of phase
    vector<OCP_DBL>          MW;
    /// label of phase
    vector<USI>              phaseLabel;
    /// molar volume of phase
    vector<OCP_DBL>          vmj;
    /// d vm / dP
    vector<OCP_DBL>          vmP;
    /// d vm / dx
    vector<vector<OCP_DBL>>  vmx;

protected:
    // Phase Function
    /// Allocate memoery for phase variables
    void AllocatePhase();
    /// Calculate molecular weight of phase
    void CalMW();
    void CalVfXiRho();
    void CalSaturation();
    USI  FindMWmax();
    /// x[j][i] -> n[j][i]
    void x2n();

public:
    // Method Function
    // Allocate memoery for Method variables
    void     AllocateMethod();
    void     PhaseEquilibrium();
    void     CalKwilson();
    OCP_BOOL PhaseStable();
    OCP_BOOL StableSSM(const USI& Id);   ///< strict SSM
    OCP_BOOL StableSSM01(const USI& Id); ///< relaxed SSM
    OCP_BOOL StableNR(const USI& Id);
    void     AssembleJmatSTA();
    OCP_BOOL CheckSplit();
    void     PhaseSplit();
    void     SplitSSM(const OCP_BOOL& flag);
    void     SplitSSM2(const OCP_BOOL& flag);
    void     SplitSSM3(const OCP_BOOL& flag);
    void     RachfordRice2();  ///< Used when NP = 2
    void     RachfordRice2P(); ///< Used when NP = 2, improved Rachford-Rice2
    void     RachfordRice3();  ///< Used when NP > 2
    void     UpdateXRR();      ///< Update X according to RR
    void     SplitBFGS();      ///< Use BFGS to calculate phase splitting
    void     SplitNR();        ///< Use NR to calculate phase splitting
    void     CalResSP();
    void     PrintFugN();
    void     AssembleJmatSP();
    OCP_DBL  CalStepNRsp();

private:
    // Method Variables
    USI testPId;                 ///< Index of the testing phase in stability analysis
    vector<vector<OCP_DBL>> Kw;  ///< Equilibrium Constant of Whilson
    vector<vector<OCP_DBL>> Ks;  ///< Approximation of Equilibrium Constant in SSM
    vector<OCP_DBL>         lKs; ///< last Ks
    vector<OCP_DBL> phiSta; ///< Fugacity coefficient used in phase stability analysis
    vector<OCP_DBL> fugSta; ///< Fugacity used in phase stability analysis
    // SSM in Stability Analysis
    vector<OCP_DBL> Y;  ///< x[i] / Yt
    OCP_DBL         Yt; ///< Sum Y
    vector<OCP_DBL> di; ///< phi(id) * x(id), id is the index of testing phase
    // NR in Stability Analysis
    vector<OCP_DBL>         resSTA;
    vector<OCP_DBL>         JmatSTA; ///< d g / d Y
    vector<vector<OCP_DBL>> fugX;    ///< d ln f / d X

    // SSM in Phase Split
    vector<OCP_DBL> resRR; ///< Error in Rachford-Rice equations.
    // NR in Phase Split
    vector<OCP_DBL> lresSP; ///< last resSP, used in BFGS
    vector<OCP_DBL> resSP;  ///< d G / d nij, G is Gibbs free energy: ln fij - ln fi,np
    vector<OCP_DBL> JmatSP; ///< Jacobian Matrix of (ln fij - ln fi,np) wrt. nij
    vector<vector<OCP_DBL>>
                    fugN;       ///< d ln fij / d nkj, in each subvector, ordered by k.
    // for linearsolve with lapack
    vector<OCP_INT> pivot;     ///< used in dgesv_ in lapack
    vector<OCP_DBL> JmatWork;  ///< work space for Jmat in STA and SP
    OCP_INT         lJmatWork; ///< length of JmatWork
    char            uplo{'U'};

protected:
    // After Phase Equilibrium Calculation finishs, properties and some auxiliary
    // variables will be calculated.
    void AllocateOthers();
    void IdentifyPhase();
    /// Copy the basic properties from MixtureComp to Mixture
    void CopyPhase();
    void CalViscosity();
    void CalViscoLBC();
    void CalViscoHZYT();

    void CalXiPNX_partial();
    void CalRhoPX_partial();
    void CalMuPX_partial();
    void CalMuPXLBC_partial();
    void CalXiRhoMuPN_full();

    // For Phase num : any
    void CalVfiVfp_full01();
    void AssembleMatVfiVfp_full01();
    void AssembleRhsVfiVfp_full01();
    void CaldXsdXpAPI01();

    // For Phase num : <=2
    void CalVfiVfp_full02();
    void AssembleMatVfiVfp_full02();
    void AssembleRhsVfiVfp_full02();
    void CaldXsdXpAPI02();

private:
    // Phase properties and auxiliary variables
    vector<OCP_DBL> muC; ///< Viscosity of phase
    vector<vector<OCP_DBL>>
        muAux; ///< Auxiliary variables for Viscosity, used to calculate Derivative
    vector<OCP_DBL>
        muAux1I; ///< Auxiliary variables for Viscosity, used to calculate Derivative
    vector<OCP_DBL>         sqrtMWi;
    vector<vector<OCP_DBL>> fugP; ///< d ln fij / d P
    vector<OCP_DBL>         Zp;   ///< d Z / d P

    vector<OCP_DBL> JmatTmp; ///< Temp Mat for transpose of a matrix
    vector<OCP_DBL> JmatDer; ///< Used to store Jacobian Mat for calculating derivates
    /// rhs or d nij / d Nk, d nij / dP in calVtiVtp
    /// rhs or dXs / dXp in Cal dXsdXp
    vector<OCP_DBL> rhsDer;

    vector<OCP_DBL> vjp; ///< dvj / dp, used in 2 hydrocarbon phase in EOS
    vector<vector<OCP_DBL>> vji; ///< dvj / dNi, used in 2 hydrocarbon phase in EOS; or dvj / dnij
    vector<OCP_DBL> xixC; ///< d xi / d xij
    vector<OCP_DBL> xiPC; ///< d xi / d P
    vector<OCP_DBL> xiNC; ///< d xi / d Nk
    vector<OCP_DBL> muN;  ///< d mu[j] / d N[i]: numphase * numCom
    vector<OCP_DBL> xiN;  ///< d xi[j] / d N[i]: numphase * numCom
    vector<OCP_DBL> rhoN; ///< d rho[j] / d N[i]: numphase * numCom

    /////////////////////////////////////////////////////////////////////
    // Optional Features
    /////////////////////////////////////////////////////////////////////

public:
    void SetupOptionalFeatures(OptionalFeatures& optFeatures) override;

    /////////////////////////////////////////////////////////////////////
    // Accelerate PVT
    /////////////////////////////////////////////////////////////////////

protected:
    /// Allocate memory for variables used in skipping stability analysis
    void AllocateSkip();
    /// Calculate d ln phi[i][j] / d n[k][j]
    void CalPhiNSkip();
    /// Assemble matrix to Calculated eigen value used for skipping
    void AssembleSkipMatSTA();
    /// Calculate skip info for next step
    void CalSkipForNextStep();
    /// Calculate Flash type for IMPEC
    void CalFtypeIMPEC() { ftype = skipSta->CalFtypeIMPEC(P, T, Nh, Ni, bulkId); }
    /// Calculate Flash type for FIM
    void CalFtypeFIM(const OCP_DBL* Sjin)
    {
        ftype = skipSta->CalFtypeFIM(P, T, Nh, Ni, Sjin, lNP, bulkId);
    }

protected:
    SkipStaAnaly* skipSta; ///< Skip analysis Term pointing to OptionalFeature

    USI ftype{0}; ///< Decide the start point of flash
    /// If ture, then skipping could be try,
    //  if ftype = 0 also, then new range should be calculated
    OCP_BOOL        flagSkip;
    vector<OCP_DBL> phiN;       ///< d ln phi[i][j] / d n[k][j]
    vector<OCP_SIN> skipMatSTA; ///< matrix for skipping Stability Analysis
    /// eigen values of matrix for skipping Skip Stability Analysis.
    //  Only the minimum eigen value will be used
    vector<OCP_SIN> eigenSkip;
    vector<OCP_SIN> eigenWork; ///< work space for computing eigenvalues with ssyevd_

    /////////////////////////////////////////////////////////////////////
    // Miscible
    /////////////////////////////////////////////////////////////////////

protected:
    void InputMiscibleParam(const ParamReservoir& rs_param, const USI& tarId);
    void CalSurfaceTension();

protected:
    Miscible*        miscible;    ///< Miscible term pointing to OptionalFeature
    USI              mIndex;      ///< method index of miscible factor calculation
    SurTenMethod01Params stm01Params; ///< params used in surface tension method01
    MisFacMethod01Params mfm01Params; ///< params used in miscible factor method01
};


#endif //__MIXTURECOMP_HEADER__
