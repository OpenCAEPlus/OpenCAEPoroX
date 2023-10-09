/*! \file    OCPNRsuite.hpp
 *  \brief   data structure used in NR iterations
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
#include "BulkVarSet.hpp"
#include "Domain.hpp"

using namespace std;

class OCPRes
{
public:
    void SetupIsoT(const OCP_USI& nb, const OCP_USI& nw, const OCP_USI& nc);
    void SetupT(const OCP_USI& nb, const OCP_USI& nw, const OCP_USI& nc);
    void SetZero();

    /// residual for all equations for each bulk
    vector<OCP_DBL> resAbs; 
    /// 2-norm of relative residual wrt. pore volume for all equations of each bulk
    vector<OCP_DBL> resRelV; 
    /// 2-norm of relative residual wrt. total moles for mass conserve equations of each bulk
    vector<OCP_DBL> resRelN; 
    /// 2-norm of relative residual wrt. total energy for energy conserve equations of each bulk
    vector<OCP_DBL> resRelE; 
    /// (initial) maximum relative residual wrt. pore volume for each bulk,
    OCP_DBL maxRelRes0_V; 
    /// (iterations) maximum relative residual wrt. pore volume for each bulk
    OCP_DBL maxRelRes_V; 
    /// maximum relative residual wrt. total moles for each bulk
    OCP_DBL maxRelRes_N;  
    /// maximum relative residual wrt. total energy for each bulk
    OCP_DBL maxRelRes_E; 
    /// maximum relative residual wrt. total moles for each well
    OCP_DBL maxWellRelRes_mol; 

    // use negative number to represent well number (ToDo)
    /// index of bulk which has maxRelRes_V
    OCP_INT maxId_V; 
    /// index of bulk which has maxRelRes_N
    OCP_INT maxId_N; 
    /// index of bulk which has maxRelRes_E
    OCP_INT maxId_E; 
};


/// NR dataset for nonlinear solution
class OCPNRsuite
{
public:
    /// Setup for themral model
    void Setup(const OCP_BOOL& ifthermal, const BulkVarSet& bvs, const OCP_USI& nw, const Domain& domain);
    /// Reset 
    void InitStep(const BulkVarSet& bvs);
    /// Calculate max change for themral model
    void CaldMax(const BulkVarSet& bvs);
    /// Get dP
    OCP_DBL DP(const OCP_USI& n) const { return dP[n]; }
    /// Get dNi
    OCP_DBL DN(const OCP_USI& n, const USI& i) const { return dN[n * nc + i]; }
    /// Get max dP
    OCP_DBL DPmax() const { return dPmax.back(); };
    /// Get max dS
    OCP_DBL DSmax() const { return dSmax.back(); };

public:
    /// residual
    OCPRes          res;

protected:
    /// Communicator
    MPI_Comm        myComm;
    OCP_INT         numproc, myrank;
    /// model
    OCP_BOOL        ifThermal{ OCP_FALSE };
    /// numBulk
    OCP_USI         nb;
    /// numPhase, numCom
    USI             np, nc;
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

    /// Max pressure difference of all NR steps within a time step
    vector<OCP_DBL> dPmax;
    /// Max temperature difference of all NR steps within a time step
    vector<OCP_DBL> dTmax;
    /// Max Ni difference of all NR steps within a time step
    vector<OCP_DBL> dNmax;
    /// Max saturation difference of all NR steps within a time step
    vector<OCP_DBL> dSmax;
};


#endif  /* end if __ISOTHERMALMETHOD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/30/2021      Create file                          */
/*----------------------------------------------------------------------------*/