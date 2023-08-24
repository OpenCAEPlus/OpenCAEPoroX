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
    friend class BulkConn;
    friend class BulkConnAreaMethod01;

public:
    /// Default constructor.
    BulkConnPair() = default;

    /// Setup BulkPair with bId and eId.
    BulkConnPair(const OCP_USI& BId,
        const OCP_USI& EId,
        const USI& direct,
        const OCP_DBL& AreaB,
        const OCP_DBL& AreaE)
        : bId(BId)
        , eId(EId)
        , direction(direct)
        , areaB(AreaB)
        , areaE(AreaE) {};

    OCP_USI BId() const { return bId; }
    OCP_USI EId() const { return eId; } 
    OCP_DBL Area() const { return area; }
    USI Type() const { return type; }
    OCP_DBL AreaB() const { return areaB; }
    OCP_DBL AreaE() const { return areaE; }
    USI     Direction() const { return direction; }

protected:
    /// first index of a pair
    OCP_USI bId; 
    /// second index of a pair.
    OCP_USI eId;
    /// Effective area
    OCP_DBL area;
    /// Connection type
    USI     type; 
    /// Connection direction
    USI     direction; 
    /// Area of intersecting face from first bulk
    OCP_DBL areaB;
    /// Area of intersecting face from second bulk
    OCP_DBL areaE;       
};


class BulkConnVarSet
{

    /////////////////////////////////////////////////////////////////////
    // General Information
    /////////////////////////////////////////////////////////////////////

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
