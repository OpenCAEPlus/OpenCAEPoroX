/*! \file    AcceleratePEC.hpp
 *  \brief   AcceleratePEC class declaration
 *  \author  Shizhe Li
 *  \date    Dec/25/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __ACCELERATEPEC_HEADER__
#define __ACCELERATEPEC_HEADER__

#include "OCPMixtureComp.hpp"

#include <vector>

using namespace std;


class SkipPSAVarset
{

public:  
    void SetNb(const USI& nbin) { nb = nbin; }
    /// Reset SkipPSA vars to last time step
    void ResetToLastTimeStep();
    /// Update SkipPSA vars at last time step
    void UpdateLastTimeStep();

public:
    OCP_BOOL         ifSetup{ OCP_FALSE };  ///< Only one setup is needed.
    OCP_USI          nb;                    ///< Num of bulk num
    USI              np;                    ///< Num of phase used in phase equilibrium calculation
    USI              nc;                    ///< Num of components used in phase equilibrium calculation

    vector<OCP_BOOL> flag;     ///< If true, skip will be test
    vector<OCP_DBL>  minEigen; ///< minimum eigenvalue used for testing skipping
    vector<OCP_DBL>  P;        ///< Pressure at last step
    vector<OCP_DBL>  T;        ///< Temperature at last step
    vector<OCP_DBL>  zi;       ///< Mole fraction of components(for test) at last step

    vector<OCP_BOOL> lflag;     ///< Last flag
    vector<OCP_DBL>  lminEigen; ///< Last min eigenvalue
    vector<OCP_DBL>  lP;        ///< Last P
    vector<OCP_DBL>  lT;        ///< Last T
    vector<OCP_DBL>  lzi;       ///< Last zi
};


class SkipPSAMethod
{
public:
    SkipPSAMethod() = default;
    /// Calculate the ftype without predicted saturations
    virtual USI CalFtype(const OCP_DBL& Pin,
                         const OCP_DBL& Tin,
                         const OCP_DBL* Niin,
                         const OCP_USI& n) = 0;
    /// Calculate the ftype with predicted saturations
    virtual USI CalFtype(const OCP_DBL& Pin,
                         const OCP_DBL& Tin,
                         const OCP_DBL* Niin,
                         const OCP_DBL* S,
                         const USI&     np,
                         const OCP_USI& n) = 0;
    /// Calculate skip info for next step
    virtual void CalSkipForNextStep(const OCP_USI& bId) = 0;

};


////////////////////////////////////////////////////////////////
// SkipPSAMethod01
////////////////////////////////////////////////////////////////


class SkipPSAMethod01 : public SkipPSAMethod
{
public:
    SkipPSAMethod01(SkipPSAVarset* vsin, const OCPMixtureComp* compsin);
    /// Calculate the ftype without predicted saturations
    USI CalFtype(const OCP_DBL& Pin,
                 const OCP_DBL& Tin,
                 const OCP_DBL* Niin,
                 const OCP_USI& bId) override;
    /// Calculate the ftype with predicted saturations
    USI CalFtype(const OCP_DBL& Pin,
                 const OCP_DBL& Tin,
                 const OCP_DBL* Niin,
                 const OCP_DBL* S,
                 const USI&     np,
                 const OCP_USI& bId) override;
    void CalSkipForNextStep(const OCP_USI& bId) override;

protected:
    OCP_BOOL IfSkip(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const OCP_DBL* Niin,
                    const OCP_USI& bId) const;

protected:
    /// pointer of variables set
    SkipPSAVarset*        vs;
    /// support modules
    const OCPMixtureComp* comps;
    
    /// d ln phi[i][j] / d n[k][j]
    vector<OCP_DBL>       lnphiN;
    /// matrix for skipping Stability Analysis,    
    vector<OCP_SIN>       skipMatSTA;
    /// eigen values of matrix for skipping Skip Stability Analysis.
    /// Only the minimum eigen value will be used
    vector<OCP_SIN>       eigenSkip;
    /// work space for computing eigenvalues with ssyevd_
    vector<OCP_SIN>       eigenWork;
};


class SkipPSA
{

public:
    /// Setup SkipPSA
    USI Setup(const OCP_USI& nb, const OCPMixtureComp* compsin);
    /// Set ifUse to true or false
    void SetUseSkip(const OCP_BOOL& flag) { ifUse = flag; }
    /// Return ifUse
    OCP_BOOL IfUseSkip() const { return ifUse; }
    /// Calculate the ftype without predicted saturations
    USI CalFtype(const OCP_DBL& Pin,
                 const OCP_DBL& Tin,
                 const OCP_DBL* Niin,
                 const OCP_USI& bId,
                 const USI&     mIndex)
    {
        if (ifUse) return sm[mIndex]->CalFtype(Pin, Tin, Niin, bId);
        else           return 0;
    }
    /// Calculate the ftype with predicted saturations
    USI CalFtype(const OCP_DBL& Pin,
                 const OCP_DBL& Tin,
                 const OCP_DBL* Niin,
                 const OCP_DBL* S,
                 const USI&     np,
                 const OCP_USI& bId,
                 const USI&     mIndex)
    {
        if (ifUse) return sm[mIndex]->CalFtype(Pin, Tin, Niin, S, np, bId);
        else           return 0;
    }
    void CalSkipForNextStep(const OCP_USI& bId, const USI& mIndex)
    {
        if (ifUse) 
            sm[mIndex]->CalSkipForNextStep(bId);        
    }
    /// Reset SkipPSA vars to last time step
    void ResetToLastTimeStep() { if(ifUse) vs.ResetToLastTimeStep(); }
    /// Update SkipPSA vars at last time step
    void UpdateLastTimeStep() { if (ifUse) vs.UpdateLastTimeStep(); }

protected:
    /// If use skipping PSA option
    OCP_BOOL                 ifUse{ OCP_TRUE };
    /// variables set
    SkipPSAVarset            vs;
    /// Skipping Method
    vector<SkipPSAMethod*>   sm;
};

#endif /* end if __ACCELERATEPEC_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/25/2022      Create file                          */
/*----------------------------------------------------------------------------*/