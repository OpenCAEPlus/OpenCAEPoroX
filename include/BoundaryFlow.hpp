/*! \file    BoundaryFlow.hpp
 *  \brief   BoundaryFlow class declaration
 *  \author  Shizhe Li
 *  \date    Sep/27/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BOUNDARYFLOW_HEADER__
#define __BOUNDARYFLOW_HEADER__


 // OpenCAEPoroX header files
#include "ParamReservoir.hpp"
#include "BulkVarSet.hpp"

#include <vector>

using namespace std;


class BoundaryFlowVarSet
{
    friend class BoundaryFlow;
    friend class BoundaryFlowMethod01;

public:

    void ResetToLastTimeStep()
    {

    }
    void UpdateLastTimeStep()
    {
    }

protected:

};


class BoundaryFlowMethod
{
public:
    BoundaryFlowMethod() = default;
};


/// Const pressure flow
class BoundaryFlowMethod01 : public BoundaryFlowMethod
{
public:
    BoundaryFlowMethod01(const BoundaryParam& bP, BoundaryFlowVarSet& bfvs) {
        name = bP.name;
        P    = bP.P;
    }

protected:
    string  name;
    OCP_DBL P;
};


class BoundaryFlow
{

    friend class BulkAccumuTerm01;
public:
    BoundaryFlow() = default;
    auto IfUse(const OCP_USI& n) const {
        if (!ifUse)              return OCP_FALSE;
        else if (mIndex[n] < 0)  return OCP_FALSE;
        else                     return OCP_TRUE;
    }
    void Setup(const ParamReservoir& rs_param, const BulkVarSet& bvs, const vector<string>& boundName, const vector<USI>& boundIndex);
    void ResetToLastTimeStep() { if (ifUse)  vs.ResetToLastTimeStep(); }
    void UpdateLastTimeStep() { if (ifUse)  vs.UpdateLastTimeStep(); }

protected:
    /// If use heat loss
    OCP_BOOL                    ifUse{ OCP_FALSE };
    /// Heat loss varsets
    BoundaryFlowVarSet          vs;
    /// Method Index
    vector<INT>                 mIndex;
    /// method for heat loss calculation
    vector<BoundaryFlowMethod*> bfM;
};


#endif /* end if __BOUNDARYFLOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/27/2023      Create file                          */
/*----------------------------------------------------------------------------*/
