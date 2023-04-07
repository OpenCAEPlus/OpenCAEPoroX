/*! \file    Reservoir.hpp
 *  \brief   Reservoir class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __RESERVOIR_HEADER__
#define __RESERVOIR_HEADER__


// OpenCAEPoro header files
#include "AllWells.hpp"
#include "Bulk.hpp"
#include "BulkConn.hpp"
#include "OptionalFeatures.hpp"
#include "ParamRead.hpp"
#include "PreProcess.hpp"

template <typename T>
class VarInfo {
public:
    VarInfo(const string n, const T* src, T* dst) :name(n), src_ptr(src), dst_ptr(dst) {};
    VarInfo(const string n, const T* src) :name(n), src_ptr(src) { dst_ptr = const_cast<T*>(src_ptr); };
    const string   name;      // varname
    const T*       src_ptr;   // ptr of src data
    T*             dst_ptr;   // ptr of dst data
};


class VarInfoGrid
{
public:
    VarInfoGrid() = default;
    VarInfoGrid(const vector<VarInfo<vector<OCP_DBL>>> dbl,
                const vector<VarInfo<vector<OCP_USI>>> usi):
        var_dbl(dbl), var_usi(usi) { CalNumByte(); };
    void CalNumByte() { 
        numvar_dbl    = var_dbl.size();
        numvar_usi    = var_usi.size();
        numByte_dbl   = sizeof(OCP_DBL);
        numByte_usi   = sizeof(OCP_USI);
        numByte_total = numvar_dbl * numByte_dbl + numvar_usi * numByte_usi; }
    vector<VarInfo<vector<OCP_DBL>>>  var_dbl;
    vector<VarInfo<vector<OCP_USI>>>  var_usi;

    OCP_INT  numvar_dbl;
    OCP_INT  numvar_usi;
    OCP_INT  numByte_dbl;
    OCP_INT  numByte_usi;
    OCP_INT  numByte_total; ///< num of byte 
};


/// Reservoir is the core component in our simulator, it contains the all reservoir
/// information, and all operations on it.
///
/// Reservoir has four Core components.
/// Grids contains the basic informations of all grids as a database of reservoir.
/// Bulk only stores active grids, which defines the area used for calculation.
/// AllWells contains the well information, it's used to manage operations related to
/// wells. BulkConn contains connections between bulks(active grids).
class Reservoir
{
    friend class OCPControl;
    friend class Summary;
    friend class CriticalInfo;
    friend class Out4RPT;
    friend class Out4VTK;

    // temp
    friend class IsoT_IMPEC;
    friend class IsoT_FIM;
    friend class IsoT_FIMn;
    friend class IsoT_AIMc;
    friend class T_FIM;
    friend class Solver;


    /////////////////////////////////////////////////////////////////////
    // Param Distribute
    /////////////////////////////////////////////////////////////////////

public:
    void InputParam(PreProcess& prepro, ParamRead& rsparam);

    /// Setup Domain
    void SetupDomain(Domain& myDomain) { swap(domain, myDomain); }
    /// Grid-based, Conn-based
    void InputDistParamGrid(ParamReservoir& param, PreParamGridWell& prepro);
    /// Well-based
    void InputDistParamOthers(const ParamRead& param);

    const Domain& GetDomain() const { return domain; }

protected:
    // Get name and address of send_var
    OCP_INT GetSendVarInfo(const VarInfoGrid& varInfo, VarInfoGrid& send_var);


    /////////////////////////////////////////////////////////////////////
    // General
    /////////////////////////////////////////////////////////////////////

public:

    /// Setup static information for reservoir with input params for Isothermal model
    void SetupIsoT();
    /// Setup static information for reservoir with input params for Thermal model
    void SetupT();
    /// Apply the control of ith critical time point.
    void ApplyControl(const USI& i);
    /// Calculate num of Injection, Production
    void CalIPRT(const OCP_DBL& dt);
    /// Calculate Maximum Change of some reference variables for IMPEC
    void CalMaxChange();
    /// Return the num of Bulk
    OCP_USI GetBulkNum() const { return bulk.GetBulkNum(); }
    /// Return the num of Bulk
    OCP_USI GetInteriorBulkNum() const { return bulk.GetInteriorBulkNum(); }
    /// Return the num of Well
    USI GetWellNum() const { return allWells.GetWellNum(); }
    /// Return the num of Components
    USI GetComNum() const { return bulk.GetComNum(); }
    /// Return num of open well
    USI GetNumOpenWell() const { return allWells.GetNumOpenWell(); }

protected:
    // Grid             grid;        ///< Init Grid info.
    Bulk             bulk;        ///< Active grid info.
    AllWells         allWells;    ///< Wells class info.
    BulkConn         conn;        ///< Bulk's connection info.
    OptionalFeatures optFeatures; ///< optional features.
    Domain           domain;      ///< domain decomposition

public:
    /// Calculate the CFL number, including bulks and wells for IMPEC
    OCP_DBL CalCFL(const OCP_DBL& dt) const;
    /// Return NRdPmax
    OCP_DBL GetNRdPmax() { return bulk.GetNRdPmax(); }
    /// Return NRdSmax
    OCP_DBL GetNRdSmax(OCP_USI& index) { return bulk.CalNRdSmax(index); }
    /// Return NRdNmax
    OCP_DBL GetNRdNmax() { return bulk.GetNRdNmax(); }
    void    PrintSolFIM(const string& outfile) const;
    void    OutInfoFinal() const { bulk.OutMixtureIters(); }
};

#endif /* end if __RESERVOIR_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/