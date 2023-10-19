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
#include "Bulk.hpp"
#include "OCPFlux.hpp"

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
    friend class T_FIM;

public:

    /////////////////////////////////////////////////////////////////////
    // General Variables
    /////////////////////////////////////////////////////////////////////

public:
    /// Input params
    void InputParam(const ParamReservoir& rs_param, const Bulk& bk);
    /// Get variable set
    auto& GetVarSet() const { return vs; }
    /// Get num of connection
    auto GetNumConn() const { return numConn; }

protected:
    OCP_USI numConn; ///< Number of connections between bulks.

    /// All connections (pair of indices) between bulks: numConn.
    //  Note: In each pair, the index of first bulk is less than the second. The data
    //  in iteratorConn is generated from neighbor.
    vector<BulkConnPair> iteratorConn;

    /// Neighboring information of each bulk: activeGridNum.
    //  Note: The i-th entry stores the i-th bulk's neighbors, which is sorted in an
    //  increasing order.
    vector<vector<OCP_USI>> neighbor;

    /////////////////////////////////////////////////////////////////////
    // Physical Variables
    /////////////////////////////////////////////////////////////////////

    BulkConnVarSet          vs;   ///< values between connections

    /////////////////////////////////////////////////////////////////////
    // Flux
    /////////////////////////////////////////////////////////////////////


    vector<OCPFlux*>        flux;    ///< flux term
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