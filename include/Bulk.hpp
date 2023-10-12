/*! \file    Bulk.hpp
 *  \brief   Bulk class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BULK_HEADER__
#define __BULK_HEADER__

// Standard header files
#include <iostream>
#include <vector>

// OpenCAEPoroX header files
#include "DenseMat.hpp"
#include "FlowUnit.hpp"
#include "LinearSystem.hpp"
#include "ParamReservoir.hpp"
#include "Domain.hpp"
#include "PreParamGridWell.hpp"
#include "BulkVarSet.hpp"
#include "PVTModule.hpp"
#include "SATModule.hpp"
#include "ROCKModule.hpp"
#include "BulkInitializer.hpp"
#include "BulkAccumuModule.hpp"
#include "OptionalModules.hpp"


using namespace std;


class BulkTypeAIM
{
public:
    void Setup(const OCP_USI& nb) { indicator.resize(nb, -1); }
    void Init(const OCP_INT& flag) { fill(indicator.begin(), indicator.end(), flag); }
    void SetBulkType(const OCP_USI& n, const OCP_INT& flag) { indicator[n] = flag; }
    OCP_BOOL IfFIMbulk(const OCP_USI& n) const { return indicator[n] >= 0; }
    OCP_BOOL IfIMPECbulk(const OCP_USI& n) const { return indicator[n] < 0; }
    OCP_INT GetBulkType(const OCP_USI& n) const { return indicator[n]; }
    OCP_USI GetNumFIMBulk() const {
        numFIMBulk = 0;
        for (const auto& v : indicator) {
            if (v >= 0) numFIMBulk++;
        }
        return numFIMBulk;
    }

protected:
    mutable OCP_USI numFIMBulk; ///< num of FIM bulk
    vector<OCP_INT> indicator; ///< Stores the index of FIM bulk in equations, 
                               ///< FIM bulk: >=0; IMPEC bulk: <0;
};

/// Physical information of each active reservoir bulk.
//  Note: Bulk contains main physical infomation of active grids. It describes the
//  actual geometric domain for simulating. Variables are stored bulk by bulk, and then
//  phase by phase, then component by component. The bulks are ordered in the alphabetic
//  order, i.e. the X-axis indices first, followed by the Y-axis and Z-axis indices.
//  Operations on each bulk are also defined here.
class Bulk
{
    friend class BulkConn;
    friend class Well;
    friend class PeacemanWell;
    friend class PeacemanWellIsoT;
    friend class PeacemanWellT;
 
    friend class Out4RPT;
    friend class Out4VTK;
    friend class OCPFlux01;
    friend class OCPFlux02;
    friend class OCPFluxT01;
    

    // temp
    friend class Reservoir;
    friend class IsothermalMethod;
    friend class IsoT_IMPEC;
    friend class IsoT_FIM;
    friend class IsoT_AIMc;
    friend class T_FIM;

    /////////////////////////////////////////////////////////////////////
    // Input Param and Setup
    /////////////////////////////////////////////////////////////////////

public:
    /// Input param from internal data structure ParamReservoir.
    void InputParam(const ParamReservoir& rs_param);

    /// Allocate memory for fluid grid for isothermal model.
    void Setup();

    /////////////////////////////////////////////////////////////////////
    // General Variables
    /////////////////////////////////////////////////////////////////////

public:
    auto& GetVarSet() const { return vs; }
    void Initialize(const Domain& domain) { INITm.Initialize(vs, PVTm, SATm, optMs, domain); }

protected:
    /// Bsaic variable set
    BulkVarSet        vs;
    /// PVT Module
    PVTModule         PVTm;
    /// SAT Module
    SATModule         SATm;
    /// Rock Module
    ROCKModule        ROCKm;
    /// Initializer
    BulkInitializer   INITm;
    /// Accumulation term
    BulkAccumuModule  ACCm;

    /// optional modules
    OptionalModules   optMs; 

public:
    /// Return the number of bulks.
    OCP_USI GetBulkNum() const { return vs.nb; }
    /// Return the number of bulks.
    OCP_USI GetInteriorBulkNum() const { return vs.nbI; }
    /// Return the number of phases.
    USI GetPhaseNum() const { return vs.np; }
    /// Return the number of components.
    USI GetComNum() const { return vs.nc; }



protected:
    OCP_DBL          rsTemp;    ///< Reservoir temperature.

    /////////////////////////////////////////////////////////////////////
    // Region
    /////////////////////////////////////////////////////////////////////

public:

    /// Output iterations in MixtureUnit
    void OutMixtureIters() const { PVTm.OutputIters(0); }


    /////////////////////////////////////////////////////////////////////
    // Basic Fluid Information
    /////////////////////////////////////////////////////////////////////

public:
    /// Calculate average pressure in reservoir.
    OCP_DBL CalFPR(OCP_DBL& vtmp) const;
    /// Calculate average Temperature in reservoir.
    OCP_DBL CalFTR(OCP_DBL& vtmp) const;
    /// Return pressure of the n-th bulk.
    OCP_DBL GetP(const OCP_USI& n) const { return vs.P[n]; }
    /// Return oil saturation of the n-th bulk.
    OCP_DBL GetSOIL(const OCP_USI& n) const
    {
        return vs.S[n * vs.np + vs.o];
    }
    /// Return gas saturation of the n-th bulk.
    OCP_DBL GetSGAS(const OCP_USI& n) const
    {
        return vs.S[n * vs.np + vs.g];
    }
    /// Return water saturation of the n-th bulk.
    OCP_DBL GetSWAT(const OCP_USI& n) const
    {
        return vs.S[n * vs.np + vs.w];
    }


    /////////////////////////////////////////////////////////////////////
    // Important Indicator Variable and Check
    /////////////////////////////////////////////////////////////////////

public:
    /// Check if negative P occurs
    ReservoirState CheckP() const;
    /// Check if negative T occurs
    ReservoirState CheckT() const;
    /// Check if negative Ni occurs
    ReservoirState CheckNi();
    /// Check if relative volume error is outranged.
    ReservoirState CheckVe(const OCP_DBL& Vlim) const;

public:
    /// push back an element for wellBulkId
    void AddWellBulkId(const OCP_USI& n) { wellBulkId.push_back(n); }

protected:
    vector<OCP_USI> wellBulkId; ///< Index of bulks which are penetrated by wells and
                                ///< their K-neighbor
    BulkTypeAIM bulkTypeAIM;
};

#endif /* end if __BULK_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/09/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/