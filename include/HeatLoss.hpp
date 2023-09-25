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
#include "ParamReservoir.hpp"
#include "BulkVarSet.hpp"

#include <vector>

using namespace std;


class HeatLossVarSet
{
    friend class HeatLoss;
    friend class HeatLossMethod01;

public:
    void SetNb(const OCP_USI& nbin) { nb = nbin; }
    void ResetToLastTimeStep()
    {
        I   = lI;
        hl  = lhl;
        hlT = lhlT;
    }
    void UpdateLastTimeStep()
    {
        lI   = I;
        lhl  = hl;
        lhlT = hlT;
    }

protected:   
    /// number of bulks
    OCP_USI         nb;
    /// Auxiliary variable
    vector<OCP_DBL> I;     
    /// heat loss rate
    vector<OCP_DBL> hl;
    /// dhl / dT
    vector<OCP_DBL> hlT;
    /// last I
    vector<OCP_DBL> lI;
    /// last hl
    vector<OCP_DBL> lhl;
    /// last hlT
    vector<OCP_DBL> lhlT;
};


class HeatLossMethod
{
public:
    HeatLossMethod() = default;
    virtual void CalHeatLoss(const OCP_USI& bId, HeatLossVarSet& hlvs, const BulkVarSet& bvs, const OCP_DBL& t, const OCP_DBL& dt) = 0;
};


class HeatLossMethod01 : public HeatLossMethod
{
public:
    HeatLossMethod01(const HLoss& param, HeatLossVarSet& hlvs);
    void CalHeatLoss(const OCP_USI& bId, HeatLossVarSet& hlvs, const BulkVarSet& bvs, const OCP_DBL& t, const OCP_DBL& dt) override;

protected:
    /// Thermal conductivity of overburden rock
    OCP_DBL         obK;
    /// Thermal conductivity of underburden rock
    OCP_DBL         ubK;
    /// Thermal diffusivity of overburden rock
    OCP_DBL         obD;
    /// Thermal diffusivity of underburden rock
    OCP_DBL         ubD;
};


class HeatLoss
{
public:
    HeatLoss() = default;
    auto IfUse() const { return ifUse; }
    void Setup(const ParamReservoir& rs_param, const OCP_USI& nb);
    void CalHeatLoss(const BulkVarSet& bvs, const OCP_DBL& t, const OCP_DBL& dt);
    const OCP_DBL& GetHl(const OCP_USI& bId) const { OCP_ASSERT(ifUse, "Inavailable!");  return vs.hl[bId]; };
    const OCP_DBL& GetHlT(const OCP_USI& bId) const { OCP_ASSERT(ifUse, "Inavailable!"); return vs.hlT[bId]; };
    void ResetToLastTimeStep() { if (ifUse)  vs.ResetToLastTimeStep(); }
    void UpdateLastTimeStep() { if (ifUse)  vs.UpdateLastTimeStep(); }

protected:
    /// If use heat loss
    OCP_BOOL                ifUse{ OCP_FALSE };
    /// Heat loss varsets
    HeatLossVarSet          vs;
    /// method for heat loss calculation
    vector<HeatLossMethod*> hlM;
};


#endif /* end if __HeatLoss_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/24/2023      Create file                          */
/*----------------------------------------------------------------------------*/
