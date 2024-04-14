/*! \file    OCPControlNR.hpp
 *  \brief   OCPControlNR class declaration
 *  \author  Shizhe Li
 *  \date    Dec/05/2023 
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPCONTROLNR_HEADER__
#define __OCPCONTROLNR_HEADER__

 // Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "ParamControl.hpp"
#include "OCPNRsuite.hpp"

using namespace std;


/// Params for Newton iterations controls
/// Note: Important for convergence of solution methods
class ControlNRParam
{
    friend class ControlNR;

public:
    ControlNRParam() = default;
    ControlNRParam(const vector<OCP_DBL>& src);

protected:
    /// Maximum number of Newton iterations in a time step
    USI     maxIter;
    /// Maximum non-linear convergence error
    OCP_DBL tol;
    /// Maximum Pressure change in a Newton iteration
    OCP_DBL dPmax;
    /// Maximum Saturation change in a Newton iteration
    OCP_DBL dSmax;
    /// Minimum Pressure change in a Newton iteration
    OCP_DBL dPmin;
    /// Minimum Saturation change in a Newton iteration
    OCP_DBL dSmin;
    /// Maximum Verr (vol error b/w fluid and pore) in a Newton step
    OCP_DBL Verrmax;
};


class ControlNR
{
public:
    /// Setup control param
    void SetCtrlParam(const vector<OCP_DBL>& src) {
        ps.push_back(ControlNRParam(src));
    }
    /// Setup communicator
    void SetupComm(const Domain& domain) {
        myComm = domain.global_comm;
        numproc = domain.numproc;
        myrank = domain.myrank;
    }
    /// Set param for next TSTEP
    void SetNextTSTEP(const USI& i) { wp = &ps[i]; }
    /// Get dSmax
    auto DSmax() const { return wp->dSmax; }
    /// Get dPmax
    auto DPmax() const { return wp->dPmax; }
    /// If NR iterations converge
    OCPNRStateC CheckConverge(const OCPNRsuite& NRs, const initializer_list<string>& il) const;

protected:
    MPI_Comm         myComm;
    OCP_INT          numproc, myrank;

protected:
    /// control param set
    vector<ControlNRParam> ps;
    /// current param
    const ControlNRParam* wp;
};

#endif /* end if __OCPControlNR_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/05/2023      Create file                          */
/*----------------------------------------------------------------------------*/