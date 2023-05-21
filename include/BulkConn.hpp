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

// Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "Bulk.hpp"
#include "OCPFlux.hpp"
#include "DenseMat.hpp"
#include "LinearSystem.hpp"
#include "OCPStructure.hpp"

using namespace std;


/// Properties and operations on connections between bulks (active grids).
//  Note: BulkConn is a core component of reservoir, it contains all properties and
//  operations on connections between bulks (active grids). You can traverse all the
//  connections through an effective iterator. Flow calculations between active
// bulks and matrix assembling with contributions from bulks only are included.
class BulkConn
{
    friend class Reservoir;
    // temp
    friend class IsoT_FIM;
    friend class IsoT_IMPEC;
    friend class IsoT_AIMc;
    friend class IsoT_FIMn;
    friend class T_FIM;

public:

    /////////////////////////////////////////////////////////////////////
    // General Variables
    /////////////////////////////////////////////////////////////////////

public:
    /// Setup active connections
    void SetupIsoT(const Bulk& bk);
    void SetupT(const Bulk& bk);
    /// Calculate the effective area used for flow
    void CalAkd(const Bulk& bk);
    void SetConnTypeIsoT(const Bulk& bk);
    void SetConnTypeT(const Bulk& bk);

protected:
    OCP_USI numConn; ///< Number of connections between bulks.

    /// All connections (pair of indices) between bulks: numConn.
    //  Note: In each pair, the index of first bulk is greater than the second. The data
    //  in iteratorConn is generated from neighbor.
    vector<BulkPair> iteratorConn;

    /// Neighboring information of each bulk: activeGridNum.
    //  Note: The i-th entry stores the i-th bulk's neighbors, which is sorted in an
    //  increasing order.
    vector<vector<OCP_USI>> neighbor;

    /////////////////////////////////////////////////////////////////////
    // Physical Variables
    /////////////////////////////////////////////////////////////////////

    BulkConnVal               bcval;   ///< values between connections

    mutable vector<OCPFlux*>  flux;    ///< flux term

protected:
    //  Note: Upblock is identified by difference of pressure between phases.
    vector<OCP_USI> upblock;          ///< Index of upwinding bulk of connections : numConn * numPhase.
    vector<OCP_DBL> upblock_Rho;      ///< Mass density of phase from upblock: numConn * numPhase.
    vector<OCP_DBL> upblock_Trans;    ///< Transmissibility of phase from upblock: numConn * numPhase.
    vector<OCP_DBL> upblock_Velocity; ///< Volume flow rate of phase from upblock: numConn * numsPhase.
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