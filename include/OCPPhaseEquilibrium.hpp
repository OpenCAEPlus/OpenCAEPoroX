/*! \file    OCPPhaseEquilibrium.hpp
 *  \brief   OCPPhaseEquilibrium class declaration
 *  \author  Shizhe Li
 *  \date    Jul/28/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPPHASEEQUILIBRIUM_HEADER__
#define __OCPPHASEEQUILIBRIUM_HEADER__

#include "OCPConst.hpp"
#include "OCPEoS.hpp"
#include "OptionalFeatures.hpp"
#include <vector>

using namespace std;


class PEIterTol
{
public:
    /// Maximum Iterations
    USI      maxIt;
    /// Tolerance
    OCP_DBL  tol;
    /// tol*tol
    OCP_DBL  tol2;
    /// residual
    OCP_DBL  res;
    /// convergence flag, if converges, conflag = OCP_TRUE
    OCP_BOOL conflag;
    /// current Iters
    USI curIt;
};

/// Params for SSM in Phase Stability Analysis
class SSMparamSTA : public PEIterTol
{
public:
    OCP_DBL  Ktol{ 1E-4 }; ///< tolerace^2 for K
    OCP_DBL  dYtol{ 1E-6 };
    OCP_DBL  eYt{ 1E-8 }; ///< if Yt > 1 + eYt, then single phase is unstable
    OCP_DBL curSk;
};

/// Params for NR in Phase Stability Analysis
class NRparamSTA : public PEIterTol { };

/// Params for SSM in Phase Split
class SSMparamSP : public PEIterTol { };

/// Params for NR in Phase Split
class NRparamSP : public PEIterTol { };

/// Param for Solving Rachford-Rice Equations
class RRparam : public PEIterTol { };


class FlashCtrl
{
public:
    SSMparamSTA SSMsta;
    NRparamSTA  NRsta;
    SSMparamSP  SSMsp;
    NRparamSP   NRsp;
    RRparam     RR;
};


class OCPPhaseEquilibrium
{
public:
	OCPPhaseEquilibrium() = default;
    /// Setup PhaseEquilibrium
    void Setup(const ComponentParam& param, const USI& tarId, EoSCalculation* eosin);
    /// PhaseEquilibrium API
    void PhaseEquilibrium(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* ziin,
                          const USI& ftypein, const OCP_DBL& lNPin, const OCP_DBL* lKsin);
    /// Get num of present phases
    const auto& GetNP() const { return NP; }
    /// return molar fraction of components in j-phase
    const auto& GetX(const USI& j) const { return x[j]; }
    /// return molar fraction of j-phase
    const auto& GetNu(const USI& j) const { return nu[j]; }
    /// Return ftype
    const auto& GetFtype() const { return ftype; }

protected:
    void SetInitalValue(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin, 
                        const USI& ftypein, const OCP_DBL& lNPin, const OCP_DBL* lKsin);

///////////////////////////////////////////////
// Basic Component properties
///////////////////////////////////////////////

protected:
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


///////////////////////////////////////////////
// Euqations of State
///////////////////////////////////////////////

protected:
    // EoS Calculations
    EoSCalculation*  eos;


///////////////////////////////////////////////
// Initial Variables
///////////////////////////////////////////////

protected: 
    /// current Pressure
    OCP_DBL         P;
    /// current Temperature
    OCP_DBL         T;
    /// mole fraction of hydrocarbon components
    vector<OCP_DBL> zi;


///////////////////////////////////////////////
// Basic variables
///////////////////////////////////////////////

protected:
    /// Allocate memoery for basic variables
    void AllocateBasicVars();
    /// x -> n
    void x2n();
    /// Calculate Molecular Weight
    void CalMW();
    /// Find the weightest phase
    USI  FindMWmax();

protected:
    /// last num of phase in last NR Step from external iterations
    USI                      lNP;
    /// current num of phase
    USI                      NP;
    /// molar fraction of phase
    vector<OCP_DBL>          nu;
    /// molar fraction of i-th component in jth phase
    vector<vector<OCP_DBL>>  x;
    /// moles of i-th component in jth phase
    vector<vector<OCP_DBL>>  n;
    /// last n in NR iterations in phase spliting calculations
    vector<vector<OCP_DBL>>  ln;
    /// Gibbs energy, before flash (not true value)
    OCP_DBL                  GibbsEnergyB;
    /// Gibbs energy, after flash (not true value)
    OCP_DBL                  GibbsEnergyE;
    /// Molecular Weight of phase
    vector<OCP_DBL>          MW;
    /// fugacity coefficient of i-th component in j-th phase
    vector<vector<OCP_DBL>>  phi;
    /// fugacity of i-th component in j-th phase
    vector<vector<OCP_DBL>>  fug;


///////////////////////////////////////////////
// Method
///////////////////////////////////////////////

protected:
    /// Allocate memory for PhaseEquilibrium
    void     AllocateMethodVars();

protected:
    // Flash method Control
    /// allowable maximum num of hydrocarbon phases
    USI             NPmax;
    /// method params for solving phase equilibrium
    FlashCtrl       flashCtrl;

protected:
    void     CalKwilson();


    ///////////////////////////////////////////////
    // Phase-Stability Analysis
    ///////////////////////////////////////////////
    
protected:
    /// hase-Stability Analysis API 
    OCP_BOOL PhaseStable();
    /// Successive Substitution Methods
    OCP_BOOL StableSSM(const USI& Id);
    /// Relaxed Successive Substitution Methods
    OCP_BOOL StableSSM01(const USI& Id);
    /// NR Methis
    OCP_BOOL StableNR(const USI& Id);
    /// Assemble Jacobian matrix for StableNR
    void     AssembleJmatSTA();

protected:
    // Method Variables
    /// Index of the testing phase in stability analysis
    USI                     testPId;     
    /// Equilibrium Constant of Whilson
    vector<vector<OCP_DBL>> Kw; 
    /// Approximation of Equilibrium Constant in SSM
    vector<vector<OCP_DBL>> Ks;  
    /// last Ks
    vector<OCP_DBL>         lKs;
    /// Fugacity coefficient used in phase stability analysis
    vector<OCP_DBL>         phiSta; 
    /// Fugacity used in phase stability analysis
    vector<OCP_DBL>         fugSta; 

    // SSM in Stability Analysis
    /// x[i] / Yt
    vector<OCP_DBL>         Y;  
    /// Sum of Y
    OCP_DBL                 Yt;
    /// phi[id] * x[id], id is the index of testing phase
    vector<OCP_DBL>         di; 

    // NR in Stability Analysis
    /// resiual of TPD equations 
    vector<OCP_DBL>         resSTA;
    ///< Jacobian matrix of TPD equations 
    vector<OCP_DBL>         JmatSTA;
    /// d ln fij / d xkj, in each subvector, ordered by k.
    vector<vector<OCP_DBL>> lnfugX;

    ///////////////////////////////////////////////
    // Phase-Splitting Calculations
    ///////////////////////////////////////////////

    /// Phase-Splitting Calculations API
    void     PhaseSplit();
    /// Successive Substitution Methods API
    void     SplitSSM(const OCP_BOOL& flag);
    /// Successive Substitution Methods for 2 phases
    void     SplitSSM2(const OCP_BOOL& flag);
    /// Successive Substitution Methods for >=3 phases
    void     SplitSSM3(const OCP_BOOL& flag);
    /// Solve RR equations with NR when NP = 2
    void     RachfordRice2();
    /// Solve RR equations with NR when NP = 2 (improved but seemd less robust)
    void     RachfordRice2P();
    /// Solve RR equations with NR when NP >= 3
    void     RachfordRice3();
    /// Update x with updated nu
    void     UpdateXRR();
    /// Solve fugacity equilibrium equations with BFGS
    void     SplitBFGS();
    /// Solve fugacity equilibrium equations with NR
    void     SplitNR();
    /// Calculate resiual of fugacity equilibrium equations 
    void     CalResSP();
    /// Assemble Jacobian matrix of fugacity equilibrium equations 
    void     AssembleJmatSP();
    /// Calculate NR-step for NR
    OCP_DBL  CalStepNRsp();
    /// Check the correctness of Phase-Splitting Calculations
    OCP_BOOL CheckSplit();

    // SSM in Phase Split
    /// resiual of Rachford-Rice equations.
    vector<OCP_DBL> resRR; 

    // NR in Phase Split
    /// resiual of fugacity equilibrium equations.
    vector<OCP_DBL>         resSP;  
    /// Jacobian matrix of fugacity equilibrium equations.
    vector<OCP_DBL>         JmatSP;
    /// last resSP, used in BFGS(to do)
    vector<OCP_DBL>         lresSP;
    /// d ln fij / d nkj, in each subvector, ordered by k.
    vector<vector<OCP_DBL>> lnfugN;

    // for linearsolve with lapack
    /// used in dgesv_ in lapack
    vector<OCP_INT>         pivot;
    /// work space for Jmat in STA and SP
    vector<OCP_DBL>         JmatWork;
    /// length of JmatWork
    OCP_INT                 lenJmatWork;
    /// stored in upper triangle when using dsysv_
    char                    uplo{ 'U' };


 ///////////////////////////////////////////////
 // Statistcs
 ///////////////////////////////////////////////

public:
    /// Print the number of iterations used in PE in total
    void OutMixtureIters() const;

protected:
    /// total iterations for SSMSTA
    OCP_ULL itersSSMSTA{ 0 };
    /// total iterations for NRSTA
    OCP_ULL itersNRSTA{ 0 };
    /// total iterations for SSMSP
    OCP_ULL itersSSMSP{ 0 };
    /// total iterations for NRSP
    OCP_ULL itersNRSP{ 0 };
    /// total iterations for RR
    OCP_ULL itersRR{ 0 };

    /// total counts of implementing SSMSTA
    OCP_ULL countsSSMSTA{ 0 };
    /// total counts of implementing NRSTA
    OCP_ULL countsNRSTA{ 0 };
    /// total counts of implementing SSMSP
    OCP_ULL countsSSMSP{ 0 };
    /// total counts of implementing NRSP
    OCP_ULL countsNRSP{ 0 };


///////////////////////////////////////////////
// Other Variables
///////////////////////////////////////////////

protected:
    /// Decide the start point of flash
    USI         ftype{ 0 };
};



#endif /* end if __OCPPHASEEQUILIBRIUM_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/28/2023      Create file                          */
/*----------------------------------------------------------------------------*/