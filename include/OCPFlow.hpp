/*! \file    OCPFlow.hpp
 *  \brief   OCPFlow class declaration
 *  \author  Shizhe Li
 *  \date    Jul/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPFLOW_HEADER__
#define __OCPFLOW_HEADER__

#include "OCPFlowVarSet.hpp"


using namespace std;


/////////////////////////////////////////////////////
// OCPFlow
/////////////////////////////////////////////////////

class OCPFlow
{
public:
    OCPFlow() = default;
    auto FlowType() const { return vs.flowType; }
    OCPFlowVarSet& GetVarSet() { return vs; }
    
    virtual OCP_DBL GetSwco() const = 0;
    virtual OCP_DBL GetMaxPcow() const = 0;
    virtual OCP_DBL GetMinPcow() const = 0;
    virtual OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const = 0;
    virtual OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const = 0;
protected:
    
    OCPFlowVarSet vs;
};


#endif /* end if __OCPFLOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/