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
#include "OCPStructure.hpp"
#include "ParamReservoir.hpp"
#include "Domain.hpp"
#include "PreParamGridWell.hpp"
#include "BulkVarSet.hpp"
#include "PVTModule.hpp"
#include "SATModule.hpp"
#include "ROCKModule.hpp"
#include "BoundaryCondition.hpp"
#include "BulkInitializer.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////
// General Error Yype
/////////////////////////////////////////////////////////////////////

const int BULK_SUCCESS                   = 0;
const int BULK_NEGATIVE_PRESSURE         = -1;
const int BULK_NEGATIVE_TEMPERATURE      = -2;
const int BULK_NEGATIVE_COMPONENTS_MOLES = -3;
const int BULK_OUTRANGED_VOLUME_ERROR    = -4;
const int BULK_OUTRANGED_CFL             = -5;




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
    friend class Out4RPT;
    friend class Out4VTK;
    friend class OCPFlux_IsoT;
    friend class OCPFlux_T;
    friend class OCPFlux_T_NF;

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
    void InputParam(const ParamReservoir& rs_param, OptionalFeatures& opts);
    void InputParamBLKOIL(const ParamReservoir& rs_param, OptionalFeatures& opts);
    void InputParamCOMPS(const ParamReservoir& rs_param, OptionalFeatures& opts);
    void InputParamTHERMAL(const ParamReservoir& rs_param, OptionalFeatures& opts);

    /// Allocate memory for fluid grid for isothermal model.
    void SetupIsoT(const Domain& domain);
    /// Allocate memory for fluid grid for ifThermal model.
    void SetupT(const Domain& domain);

    /////////////////////////////////////////////////////////////////////
    // General Variables
    /////////////////////////////////////////////////////////////////////

public:
    auto& GetVarSet() const { return vs; }
    void Initialize(const Domain& domain) { INITm.Initialize(vs, PVTm, SATm, domain); }

protected:
    /// Bsaic variable set
    BulkVarSet        vs;
    /// PVT Module
    PVTModule         PVTm;
    /// SAT Module
    SATModule         SATm;
    /// Rock Module
    ROCKModule        ROCKm;
    /// Boundary condition Module
    BoundaryCondition BCm;

    /// Initialize
    BulkInitializer   INITm;

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
    vector<OCP_DBL>  thconp;    ///< Phase thermal conductivity: numPhase

    /////////////////////////////////////////////////////////////////////
    // Region
    /////////////////////////////////////////////////////////////////////

public:

    /// Output iterations in MixtureUnit
    void OutMixtureIters() const { PVTm.GetPVT(0)->OutMixtureIters(); }


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
        return vs.S[n * vs.np + vs.oIndex];
    }
    /// Return gas saturation of the n-th bulk.
    OCP_DBL GetSGAS(const OCP_USI& n) const
    {
        return vs.S[n * vs.np + vs.gIndex];
    }
    /// Return water saturation of the n-th bulk.
    OCP_DBL GetSWAT(const OCP_USI& n) const
    {
        return vs.S[n * vs.np + vs.wIndex];
    }


    /////////////////////////////////////////////////////////////////////
    // Important Indicator Variable and Check
    /////////////////////////////////////////////////////////////////////

public:
    /// Check if negative P occurs
    OCP_INT CheckP() const;
    /// Check if negative T occurs
    OCP_INT CheckT() const;
    /// Check if negative Ni occurs
    OCP_INT CheckNi();
    /// Check if relative volume error is outranged.
    OCP_INT CheckVe(const OCP_DBL& Vlim) const;
    /// Check if Cfl is outranged.
    OCP_INT CheckCFL(const OCP_DBL& cflLim) const;

    /// Calculate max change of indicator variables.
    void CalMaxChange();

    /// Return dPmax.
    OCP_DBL GetdPmax() const { return dPmax; }
    /// Return dTmax
    OCP_DBL GetdTmax() const { return dTmax; }
    /// Return dNmax.
    OCP_DBL GetdNmax() const { return dNmax; }
    /// Return dSmax.
    OCP_DBL GetdSmax() const { return dSmax; }
    /// Return eVmax.
    OCP_DBL GeteVmax() const { return eVmax; }
    /// Return maxCFL
    OCP_DBL GetMaxCFL() const { return maxCFL; }

protected:
    OCP_DBL dPmax; ///< Max change in pressure during the current time step.
    OCP_DBL dTmax; ///< Max change in temperature during the current time step.
    OCP_DBL dSmax; ///< Max change in saturation during the current time step.
    OCP_DBL dNmax; ///< Max change in moles of component during the current time step.
    OCP_DBL eVmax; ///< Max relative diff between fluid and pore volume during the
                   ///< current time step.

    mutable vector<OCP_DBL> cfl;             ///< CFL number for each bulk
    mutable OCP_DBL         maxCFL{ 0 };     ///< max CFL number
    mutable OCP_DBL         maxCFL_loc{ 0 }; ///< local maxCFL


public:
    /// push back an element for wellBulkId
    void AddWellBulkId(const OCP_USI& n) { wellBulkId.push_back(n); }

protected:
    vector<OCP_USI> wellBulkId; ///< Index of bulks which are penetrated by wells and
                                ///< their K-neighbor
    BulkTypeAIM bulkTypeAIM;
    vector<OCP_DBL> xijNR; ///< store the current NR step's xij in AIM
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