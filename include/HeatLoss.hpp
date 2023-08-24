/*! \file    HeatLoss.hpp
 *  \brief   HeatLoss class declaration
 *  \author  Shizhe Li
 *  \date    Aug/24/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __HEATLOSS_HEADER__
#define __HEATLOSS_HEADER__


 // OpenCAEPoroX header files
#include "BulkVarSet.hpp"


#include <vector>

using namespace std;


class HeatLossVarSet
{
public:
    /// Auxiliary variable
    vector<OCP_DBL> I;     
    /// Auxiliary variable
    vector<OCP_DBL> p;
    /// dP / dT
    vector<OCP_DBL> pT;
    /// last I
    vector<OCP_DBL> lI;
    /// last p
    vector<OCP_DBL> lp;
    /// last pT
    vector<OCP_DBL> lpT;
};


class HeatLossMethod
{
public:
    HeatLossMethod() = default;
    virtual void CalHeatLoss(const BulkVarSet& bvs) = 0;
};


class HeatLossMethod01 : public HeatLossMethod
{
public:
    HeatLossMethod01() = default;
    void CalHeatLoss(const BulkVarSet& bvs) override;
};


class HeatLoss01
{
public:
    HeatLoss01() = default;
    void Setup();
    void CalHeatLoss(const BulkVarSet& bvs) { hlM->CalHeatLoss(bvs); }

protected:
    /// Heat loss varsets
    HeatLossVarSet  vs;
    /// method for heat loss calculation
    HeatLossMethod* hlM;
};


#endif /* end if __HeatLoss_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/24/2023      Create file                          */
/*----------------------------------------------------------------------------*/
