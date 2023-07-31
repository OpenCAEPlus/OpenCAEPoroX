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
#include "OCPPhaseEquilibrium.hpp"

using namespace std;


class MixtureComp : public MixtureUnit
{

public:
    void    OutMixtureIters() const override;


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
            wIdP = numPhase - 1;
            wIdC = numCom - 1;
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
        copy(Ziin, Ziin + NC, zi.begin());
        Nh = 1;
    }
    void InitPTN(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
    {
        P = Pin;
        T = Tin;
        Nh = 0;
        for (USI i = 0; i < NC; i++) {
            Ni[i] = Niin[i];
            Nh    += Ni[i];
        }
        Ni[wIdC] = Niin[wIdC];
        Nt       = Nh + Ni[wIdC];
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
    /// maximum number of hydrocarbon phase
    USI             NPmax;
    /// number of hydrocarbon components
    USI             NC;
    /// number of existing hydrocarbon phase
    USI             NP;
    /// molar fraction for hydrocarbon components
    vector<OCP_DBL> zi;
    OCP_DBL         Nh;
    USI             lNP;
    vector<OCP_DBL> MW;
    vector<OCP_DBL> vm;
    vector<OCP_DBL> vmP;
    vector<vector<OCP_DBL>> vmx;
    // for linearsolve with lapack
    /// used in dgesv_ in lapack
    vector<OCP_INT>         pivot;
    vector<USI>  epIndex;
    vector<USI>  epIndexH;
    /// molecular Weight of hydrocarbon components
    vector<OCP_DBL> MWC;
    /// critical temperature of hydrocarbon components
    vector<OCP_DBL> Tc;
    /// critical pressure of hydrocarbon components
    vector<OCP_DBL> Pc;
    /// critical volume of hydrocarbon components
    vector<OCP_DBL> Vc;
    /// phase labels, used to identify phase
    vector<USI> phaseLabel;

protected:
    void CalPropertyW(const OCP_DBL& vw);

protected:
    // water property
    /// PVTW
    OCP_PVTW        PVTW;
    /// mass density of water phase in standard condition.
    OCP_DBL         std_RhoW;
    /// index of water phase
    USI             wIdP;
    /// index of water component
    USI             wIdC;

protected:
    /// Phase equilibrium calculations
    OCPPhaseEquilibrium PE;
    // EoS Calculations
    EoSCalculation      eos;

protected:
    // viscosity calculations
    ViscosityCalculation visCal;


protected:
    /// Calculate molecular weight of phase   
    void CalVfS();

protected:

    /// Copy phase properties from PhaseEquilibrium
    void CopyPhaseFromPE();
    /// Calculate Molecular Weight
    void CalMW();
    /// Calculate moalr volume and volume of phase
    void CalVmVj();
    /// Calculate phase property
    void CalPropertyH();
    /// Calculate phases' molar density, mass density and viscosity
    void CalXiRhoMu();
    /// Calculate phase property and derivatives
    void CalPropertyHDer();
    /// Calculate phases' molar density, mass density and viscosity and derivatives
    void CalXiRhoMuDer();


protected:
    /// Identify Phase
    void IdentifyPhase();
    /// ReOrder Phase
    void ReOrderPhase();
    /// ReOrder Phase and ders
    void ReOrderPhaseDer();

protected:
    /// work space for reorder
    vector<OCP_DBL> rowork;


protected:
    // After Phase Equilibrium Calculation finishs, properties and some auxiliary
    // variables will be calculated.
    void AllocateOthers();
    
    

    // For Phase num : any
    void CalVfiVfp_full01();
    void AssembleMatVfiVfp_full01();
    void AssembleRhsVfiVfp_full01();
    void CaldXsdXp01();

    // For Phase num : <=2
    void CalVfiVfp_full02();
    void AssembleMatVfiVfp_full02();
    void AssembleRhsVfiVfp_full02();
    void CaldXsdXp02();

private:
    vector<vector<OCP_DBL>> lnfugP; ///< d ln fij / d P
    vector<vector<OCP_DBL>> lnfugN; ///< d ln fij / d nkj

    vector<OCP_DBL> JmatDer; ///< Used to store Jacobian Mat for calculating derivates
    /// rhs or d nij / d Nk, d nij / dP in calVtiVtp
    /// rhs or dXs / dXp in Cal dXsdXp
    vector<OCP_DBL> rhsDer;

    vector<OCP_DBL> vjp; ///< dvj / dp, used in 2 hydrocarbon phase in EOS
    vector<vector<OCP_DBL>> vji; ///< dvj / dNi, used in 2 hydrocarbon phase in EOS; or dvj / dnij

    /////////////////////////////////////////////////////////////////////
    // Optional Features
    /////////////////////////////////////////////////////////////////////

public:
    void SetupOptionalFeatures(OptionalFeatures& optFeatures) override;


    ///////////////////////////////////////////////
    // Skip Phase-Stability Analysis
    ///////////////////////////////////////////////

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
    /// Skip analysis Term pointing to OptionalFeature
    SkipStaAnaly* skipSta;
    /// Decide the start point of flash
    USI ftype{ 0 };
    ///  if try to skip
    OCP_BOOL        flagSkip;
    /// d ln phi[i][j] / d n[k][j]
    vector<OCP_DBL> lnphiN;
    /// matrix for skipping Stability Analysis,    
    vector<OCP_SIN> skipMatSTA;
    /// eigen values of matrix for skipping Skip Stability Analysis.
    /// Only the minimum eigen value will be used
    vector<OCP_SIN> eigenSkip;
    /// work space for computing eigenvalues with ssyevd_
    vector<OCP_SIN> eigenWork;



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
