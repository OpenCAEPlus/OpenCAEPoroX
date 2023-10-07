/*! \file    OCPStructure.hpp
 *  \brief   Some Structure in OpenCAEPoroX
 *  \author  Shizhe Li
 *  \date    Oct/30/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPSTRUCTURE_HEADER__
#define __OCPSTRUCTURE_HEADER__

// Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "BulkVarSet.hpp"

using namespace std;

class OCPRes
{
public:
    void SetupIsoT(const OCP_USI& nb, const OCP_USI& nw, const OCP_USI& nc)
    {
        OCP_USI reslen = (nb + nw) * (nc + 1);
        resAbs.resize(reslen);
        resRelV.resize(nb);
        resRelN.resize(nb);
    }
    void SetupT(const OCP_USI& nb, const OCP_USI& nw, const OCP_USI& nc)
    {
        OCP_USI reslen = (nb + nw) * (nc + 2);
        resAbs.resize(reslen);
        resRelV.resize(nb);
        resRelN.resize(nb);
        resRelE.resize(nb);
    }
    void SetZero()
    {
        fill(resAbs.begin(), resAbs.end(), 0);
        fill(resRelV.begin(), resRelV.end(), 0);
        fill(resRelN.begin(), resRelN.end(), 0);
        fill(resRelE.begin(), resRelE.end(), 0);
        maxRelRes_V       = 0;
        maxRelRes_N       = 0;
        maxRelRes_E       = 0;
        maxWellRelRes_mol = 0;
        maxId_V           = 0;
        maxId_N           = 0;
        maxId_E           = 0;
    }
    void SetInitRes() { maxRelRes0_V = maxRelRes_V; }

    /// residual for all equations for each bulk
    vector<OCP_DBL> resAbs; 
    /// 2-norm of relative residual wrt. pore volume for all equations of each bulk
    vector<OCP_DBL> resRelV; 
    /// 2-norm of relative residual wrt. total moles for mass conserve equations of each bulk
    vector<OCP_DBL> resRelN; 
    /// 2-norm of relative residual wrt. total energy for energy conserve equations of each bulk
    vector<OCP_DBL> resRelE; 
    /// (initial) maximum relative residual wrt. pore volume for each bulk
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
class OCPNRsuit
{
public:
    void SetupIsoT(const OCP_USI& nb, const OCP_USI& nw, const USI& np, const USI& nc) {
        res.SetupIsoT(nb, nw, nc);
        lP.resize(nb);
    }
    void SetupT(const OCP_USI& nb, const OCP_USI& nw, const USI& np, const USI& nc) {
        res.SetupT(nb, nw, nc);
    }
public:
    /// residual of nonlinear equations
    OCPRes          res;
    /// P at last NR steps
    vector<OCP_DBL> lP;
    /// T at last NR steps
    vector<OCP_DBL> lT;
    /// Ni at last NR steps
    vector<OCP_DBL> lN;
    /// S at last NR steps
    vector<OCP_DBL> lS;
    /// P change between NR steps
    vector<OCP_DBL> dPNR;
    /// T change between NR steps
    vector<OCP_DBL> dTNR;
    /// Ni change between NR steps
    vector<OCP_DBL> dNNR;
    /// saturation change between NR steps
    vector<OCP_DBL> dSNR;

    /// Max pressure difference in an NR step
    OCP_DBL         NRdPmax;
    /// Max temperature difference in an NR step
    OCP_DBL         NRdTmax;
    /// Max Ni difference in an NR step
    OCP_DBL         NRdNmax;
    /// Max saturation difference in an NR step(Real)
    OCP_DBL         NRdSmax;
};


#endif