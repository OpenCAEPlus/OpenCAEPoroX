/*! \file    OCPNRsuite.hpp
 *  \brief   data structure used in non-linear iterations
 *  \author  Shizhe Li
 *  \date    Oct/30/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPNRSUITE_HEADER__
#define __OCPNRSUITE_HEADER__

// Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "Domain.hpp"
#include "Reservoir.hpp"

using namespace std;


/// NR dataset for nonlinear solution
class OCPNRsuite
{
public:
    /// Setup for method which use NR
    void Setup(const OCP_BOOL& ifthermal, const BulkVarSet& bvs, const OCP_USI& nw, const Domain& domain);
    /// Setup for method which does not use NR
    void Setup(const BulkVarSet& bvs, const Domain& domain);
    /// Reset 
    void InitStep(const BulkVarSet& bvs);
    /// Calculate max change for themral model
    void CalMaxChangeNR(const Reservoir& rs);
    /// Get dP
    OCP_DBL DP(const OCP_USI& n) const { return dP[n]; }
    /// Get dNi
    OCP_DBL DN(const OCP_USI& n, const USI& i) const { return dN[n * nc + i]; }
    /// Get current max dP
    OCP_DBL DPBmaxNRc() const { return dPBmaxNR.back(); };
    /// Get current max dS
    OCP_DBL DSmaxNRc() const { return dSmaxNR.back(); };
    /// Get all bulk max dP
    const auto& DPBmaxNR() const { return dPBmaxNR; };
    /// Get all max dT
    const auto& DTmaxNR() const { return dTmaxNR; };
    /// Get all max dN
    const auto& DNmaxNR() const { return dNmaxNR; };
    /// Get all max dS
    const auto& DSmaxNR() const { return dSmaxNR; };
    /// Get all well max dP
    const auto& DPWmaxNR() const { return dPWmaxNR; };

public:
    /// residual
    OCPNRresidual   res;

protected:
    /// Communicator
    MPI_Comm        myComm;
    OCP_INT         numproc, myrank;
    /// if use NR
    OCP_BOOL        ifUseNR{ OCP_FALSE };
    /// model
    OCP_BOOL        ifThermal{ OCP_FALSE };
    /// numBulk
    OCP_USI         nb;
    /// numPhase, numCom
    USI             np, nc;

    // between NR-step

    /// P at last NR steps
    vector<OCP_DBL> lP;
    /// T at last NR steps
    vector<OCP_DBL> lT;
    /// Ni at last NR steps
    vector<OCP_DBL> lN;
    /// S at last NR steps
    vector<OCP_DBL> lS;
    /// P change between NR steps
    vector<OCP_DBL> dP;
    /// T change between NR steps
    vector<OCP_DBL> dT;
    /// Ni change between NR steps
    vector<OCP_DBL> dN;
    /// saturation change between NR steps
    vector<OCP_DBL> dS;

    /// Max pressure difference of all NR steps within a time step (bulk)
    vector<OCP_DBL> dPBmaxNR;
    /// Max temperature difference of all NR steps within a time step
    vector<OCP_DBL> dTmaxNR;
    /// Max Ni difference of all NR steps within a time step
    vector<OCP_DBL> dNmaxNR;
    /// Max saturation difference of all NR steps within a time step
    vector<OCP_DBL> dSmaxNR;
    /// Max pressure difference of all NR steps within a time step (well)
    vector<OCP_DBL> dPWmaxNR;

    // CFL
public:
    /// Calculate CFL number
    void CalCFL(const Reservoir& rs, const OCP_DBL& dt, const OCP_BOOL& ifComm);
    /// Get maxCFL
    OCP_DBL GetMaxCFL() const { return maxCFL; }
    /// Get CFL
    const OCP_DBL& GetCFL(const OCP_USI& n, const USI& j) const { return cfl[n * np + j]; }
    /// Check CFL
    OCP_BOOL CheckCFL(const OCP_DBL& cflLim) const;

protected:
    /// CFL number for each bulk
    vector<OCP_DBL> cfl;  
    /// max CFL number for global
    OCP_DBL         maxCFL{ 0 };
    /// local maxCFL
    OCP_DBL         maxCFL_loc{ 0 }; 


    // between time step
public:
    /// calculate maximum change between time step
    void CalMaxChangeTime(const Reservoir& rs);
    /// Return dPmax.
    OCP_DBL DPmaxT() const { return dPmaxT; }
    /// Return dPmax.
    OCP_DBL DPBmaxT() const { return dPBmaxT; }
    /// Return dPmax.
    OCP_DBL DPWmaxT() const { return dPWmaxT; }
    /// Return dTmax
    OCP_DBL DTmaxT() const { return dTmaxT; }
    /// Return dNmax.
    OCP_DBL DNmaxT() const { return dNmaxT; }
    /// Return dSmax.
    OCP_DBL DSmaxT() const { return dSmaxT; }
    /// Return eVmax.
    OCP_DBL EVmaxT() const { return eVmaxT; }


protected:
    /// Max change in pressure during the current time step.(bulk and well)
    OCP_DBL dPmaxT; 
    /// Max change in pressure during the current time step.(bulk and well)
    OCP_DBL dPBmaxT;
    /// Max change in pressure during the current time step.(bulk and well)
    OCP_DBL dPWmaxT;
    /// Max change in temperature during the current time step.
    OCP_DBL dTmaxT; 
    /// Max relative change in moles of component during the current time step.
    OCP_DBL dNmaxT;
    /// Max change in saturation during the current time step.
    OCP_DBL dSmaxT; 
    /// Max relative diff between fluid and pore volume during the current time step.
    OCP_DBL eVmaxT; 


    // Iterations
public:
    /// initialize iters when begining a new time step
    void InitIter();
    /// update iters when finishing a LS step
    void UpdateIter(const USI& lsIter);
    /// Reset iters when resetting time step
    void ResetIter();
    /// Get iterNR
    auto GetIterNR() const { return iterNR; }
    /// Get iterLS
    auto GetIterLS() const { return iterLS; }
    /// Get iterNRw
    auto GetIterNRw() const { return iterNRw; }
    /// Get iterLSw
    auto GetIterLSw() const { return iterLSw; }
    /// Get iterNRLS
    const auto& GetIterNRLS() const { return iterNRLS; };

protected:
    /// number of Newton-Raphson iterations within a time step
    USI         iterNR{ 0 };
    /// number of linear solver iterations within a time step
    USI         iterLS{ 0 };
    /// wasted number of Newton-Raphson iterations within a time step
    USI         iterNRw{ 0 };
    /// wasted number of linear solver iterations within a time step
    USI         iterLSw{ 0 };
    /// LS iterations of each NR within a time step
    vector<USI> iterNRLS;
};


#endif  /* end if __OCPNRSUITE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/30/2021      Create file                          */
/*----------------------------------------------------------------------------*/