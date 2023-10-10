/*! \file    OCPNRresidual.hpp
 *  \brief   data structure used in non-linear iterations
 *  \author  Shizhe Li
 *  \date    Oct/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPNRRESIDUAL_HEADER__
#define __OCPNRRESIDUAL_HEADER__

 // Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"

using namespace std;

class OCPNRresidual
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



#endif  /* end if __OCPNRRESIDUAL_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/