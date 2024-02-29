/*! \file    OCPDiffusion.hpp
 *  \brief   OCPDiffusion class declaration
 *  \author  Shizhe Li
 *  \date    Feb/28/2024
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPDIFFUSION_HEADER__
#define __OCPDIFFUSION_HEADER__


 // OpenCAEPoroX header files
#include "ParamReservoir.hpp"
#include "BulkVarSet.hpp"
#include "BulkConnVarSet.hpp"
#include "BulkConnFunc.hpp"

#include <vector>

using namespace std;


class OCPDiffusionVarSet
{
    friend class OCPDiffusion;
    friend class OCPDiffusionMethod01;

public:
    void SetUp(const USI& _np, const USI& _nc, const INT& oIndex, const INT& gIndex, const INT& wIndex) {
        np      = _np;
        nc      = _nc;
        o       = oIndex;
        g       = gIndex;
        w       = wIndex;
    }
    void ResetToLastTimeStep()
    {

    }
    void UpdateLastTimeStep()
    {

    }

public:
    /// number of bulk connections
    OCP_USI         numConn;
    /// number of phases
    USI             np;
    /// number of components
    USI             nc;
    /// oil, gas, water index
    INT             o, g, w;
};


class OCPDiffusionMethod
{
public:
    OCPDiffusionMethod() = default;

    virtual void CalFlux(const BulkConnPair& bp, const OCPDiffusionVarSet& dvs, const BulkVarSet& bvs, FluxVarSet& fvs) const = 0;
    virtual void AssembleFIM(const BulkConnPair& bp, const OCPDiffusionVarSet& dvs, const BulkVarSet& bvs, FluxVarSet& fvs) const = 0;

public:
    /// Calculate diffusion
    void CalDiffu(BulkConnPair& bp, const Bulk& bk) { bcd.CalDiffu(bp, bk); }
protected:
    BulkConnDiffu   bcd;
};


class OCPDiffusionMethod01 : public OCPDiffusionMethod
{
public:
    OCPDiffusionMethod01(const ParamReservoir& rs_param, OCPDiffusionVarSet& dvs);
    void CalFlux(const BulkConnPair& bp, const OCPDiffusionVarSet& dvs, const BulkVarSet& bvs, FluxVarSet& fvs) const override;
    void AssembleFIM(const BulkConnPair& bp, const OCPDiffusionVarSet& dvs, const BulkVarSet& bvs, FluxVarSet& fvs) const override;

protected:
    /// diffusion coefficient of components in phases
    vector<vector<OCP_DBL>>   diffuCP;
};


class OCPDiffusion
{
public:
    OCPDiffusion() = default;
    void Setup(const ParamReservoir& rs_param, const BulkVarSet& bvs);
    /// Calculate diffusion
    void CalDiffu(BulkConnPair& bp, const Bulk& bk);
    /// Calculate the component flux caused by diffusion
    void CalFlux(const BulkConnPair& bp, const BulkVarSet& bvs, FluxVarSet& fvs);
    const auto& GetVarSet() const { return vs; }
    void ResetToLastTimeStep() { if (ifUse)  vs.ResetToLastTimeStep(); }
    void UpdateLastTimeStep() { if (ifUse)  vs.UpdateLastTimeStep(); }


protected:
    /// If use mass diffusion
    OCP_BOOL                    ifUse{ OCP_FALSE };
    /// mass diffusion varsets
    OCPDiffusionVarSet          vs;
    /// method for mass diffusion calculation
    vector<OCPDiffusionMethod*> dM;
};


#endif /* end if __OCPDiffusion_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Feb/28/2024      Create file                          */
/*----------------------------------------------------------------------------*/
