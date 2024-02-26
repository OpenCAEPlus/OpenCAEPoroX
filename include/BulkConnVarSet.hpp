/*! \file    BulkConnVarSet.hpp
 *  \brief   BulkConnVarSet class declaration
 *  \author  Shizhe Li
 *  \date    Aug/23/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BULKCONNVARSET_HEADER__
#define __BULKCONNVARSET_HEADER__


 // OpenCAEPoroX header files
#include "OCPConst.hpp"


#include <vector>

using namespace std;


/// Connection between two bulks (bId, eId); usually, indices bId > eId.
//  Note: Bulks are the active grid cells.
class BulkConnPair
{
    friend class BulkConnTransMethod01;
    friend class BulkConnDiffuMethod01;

public:
    /// Default constructor.
    BulkConnPair() = default;

    /// Setup BulkPair with bId and eId.
    BulkConnPair(const OCP_USI& BId,
        const OCP_USI& EId,
        const ConnDirect& direct,
        const OCP_DBL& AreaB,
        const OCP_DBL& AreaE, 
        const OCP_DBL& transM)
        : bId(BId)
        , eId(EId)
        , direction(direct)
        , areaB(AreaB)
        , areaE(AreaE)
        , transMult(transM) {};
    const auto& BId() const { return bId; }
    const auto& EId() const { return eId; } 
    const auto& Trans() const { return trans; }
    const auto& Diffu() const { return diffu; }
    const auto& AreaB() const { return areaB; }
    const auto& AreaE() const { return areaE; }
    const auto& Direction() const { return direction; }

protected:
    /// first index of a pair
    OCP_USI    bId; 
    /// second index of a pair.
    OCP_USI    eId;
    /// Effective area for permeability
    OCP_DBL    trans;
    /// Effective area for diffusity
    OCP_DBL    diffu;
    /// Connection direction
    ConnDirect direction;
    /// Area of intersecting face from first bulk
    OCP_DBL    areaB;
    /// Area of intersecting face from second bulk
    OCP_DBL    areaE;
    /// Transmissibility multipliers
    OCP_DBL    transMult{ 1.0 };
};


class BulkConnVarSet
{

    /////////////////////////////////////////////////////////////////////
    // General Information
    /////////////////////////////////////////////////////////////////////

public:


public:
    /// Index of upwinding bulk of connections for each phase
    vector<OCP_USI> upblock;
    /// Mass density of phase for connections
    vector<OCP_DBL> rho;
    /// Volume flow rate of phase from upblock
    vector<OCP_DBL> flux_vj;
    /// mole flow rate of components 
    vector<OCP_DBL> flux_ni;  

    /// last upblock
    vector<OCP_USI> lupblock;
    /// last rho
    vector<OCP_DBL> lrho;      

};



#endif /* end if __BulkConnVarSet_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/23/2023      Create file                          */
/*----------------------------------------------------------------------------*/
