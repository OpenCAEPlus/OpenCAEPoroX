/*! \file    BulkConn.hpp
 *  \brief   BulkConn class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BULKCONN_HEADER__
#define __BULKCONN_HEADER__

// OpenCAEPoroX header files
#include "FLUXModule.hpp"

using namespace std;


/// Properties and operations on connections between bulks (active grids).
//  Note: BulkConn is a core component of reservoir, it contains all properties and
//  operations on connections between bulks (active grids). You can traverse all the
//  connections through an effective iterator. Flow calculations between active
// bulks and matrix assembling with contributions from bulks only are included.
class BulkConn
{
    friend class Reservoir;
    friend class NRsuite;
    // temp
    friend class IsoT_FIM;
    friend class IsoT_IMPEC;
    friend class IsoT_AIMc;
    friend class IsoT_FIMddm;
    friend class T_FIM;

public:

    /////////////////////////////////////////////////////////////////////
    // General Variables
    /////////////////////////////////////////////////////////////////////

public:
    /// Input params
    void InputParam(const ParamReservoir& rs_param, const Bulk& bk) {
        vs.numConn = numConn;
        optMs.Setup(rs_param, bk);
        FLUXm.InputParam(rs_param, iteratorConn, bk, optMs);      
    }
    /// Get variable set
    auto& GetVarSet() const { return vs; }
    /// Get num of connection
    auto GetNumConn() const { return numConn; }
    /// Calculate transmissibility for all connections
    void CalTrans(const Bulk& bk) {
        for (OCP_USI c = 0; c < numConn; c++) {
            FLUXm.GetFlux(c)->CalTrans(iteratorConn[c], bk);
        }
    }
    /// Calculate diffusity for all connections
    void CalDiffu(const Bulk& bk) {
        for (OCP_USI c = 0; c < numConn; c++) {
            FLUXm.GetFlux(c)->CalDiffu(iteratorConn[c], bk);
        }
    }

protected:

    /// Number of connections between bulks.
    OCP_USI numConn; 

    /// All connections (pair of indices) between bulks: numConn.
    //  Note: In each pair, the index of first bulk is less than the second. The data
    //  in iteratorConn is generated from neighbor.
    vector<BulkConnPair>    iteratorConn;

    /// Neighboring information of each bulk: activeGridNum.
    //  Note: The i-th entry stores the i-th bulk's neighbors, which is sorted in an
    //  increasing order.
    vector<vector<OCP_USI>> neighbor;

    /////////////////////////////////////////////////////////////////////
    // Physical Variables
    /////////////////////////////////////////////////////////////////////

    /// basic varset for bulk connections
    BulkConnVarSet          vs;   

    /////////////////////////////////////////////////////////////////////
    // Flux
    /////////////////////////////////////////////////////////////////////

    /// flux term
    FLUXModule              FLUXm; 
    /// optional modules
    BulkConnOptionalModules optMs;  
};

#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/17/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/