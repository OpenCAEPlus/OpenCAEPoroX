/*! \file    OCPFlux.hpp
 *  \brief   OCPFlux class declaration
 *  \author  Shizhe Li
 *  \date    May/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPFLUX_HEADER__
#define __OCPFLUX_HEADER__

 // Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "Bulk.hpp"


/// Connection between two bulks (bId, eId); usually, indices bId > eId.
//  Note: Bulks are the active grid cells.
class BulkPair
{
    friend class BulkConn;

public:
    /// Default constructor.
    BulkPair() = default;

    /// Setup BulkPair with bId and eId.
    BulkPair(const OCP_USI& BId,
        const OCP_USI& EId,
        const USI& direct,
        const OCP_DBL& AreaB,
        const OCP_DBL& AreaE)
        : bId(BId)
        , eId(EId)
        , direction(direct)
        , areaB(AreaB)
        , areaE(AreaE) {};

    OCP_USI BId() const { return bId; }   ///< Return beginning index.
    OCP_USI EId() const { return eId; }   ///< Return ending index.
    OCP_DBL Area() const { return area; } ///< Return effective area
    OCP_DBL AreaB() const { return areaB; }
    OCP_DBL AreaE() const { return areaE; }
    USI     Direction() const { return direction; }

protected:
    OCP_USI bId;       ///< Beginning index of a pair.
    OCP_USI eId;       ///< Ending index of a pair.
    OCP_DBL area;      ///< Effective area
    USI     direction; ///< 1-x, 2-y, 3-z
    OCP_DBL areaB;     ///< Area of intersecting faces from Begin grid
    OCP_DBL areaE;     ///< Area of intersecting faces from End grid
};


// Physical Variables of BulkConn
class BulkConnVal
{
protected:
    //  Note: Upblock is identified by difference of pressure between phases.
    vector<OCP_USI>
        upblock; ///< Index of upwinding bulk of connections : numConn * numPhase.
    vector<OCP_DBL>
        upblock_Rho; ///< Mass density of phase from upblock: numConn * numPhase.
    vector<OCP_DBL>
        upblock_Trans; ///< Transmissibility of phase from upblock: numConn * numPhase.
    vector<OCP_DBL> upblock_Velocity; ///< Volume flow rate of phase from upblock:
                                      ///< numConn * numsPhase.
    vector<OCP_DBL> Adkt;             ///< Thermal conductivity between neighbors

    // Last time step
    vector<OCP_USI> lupblock;          ///< last upblock
    vector<OCP_DBL> lupblock_Rho;      ///< last upblock_Rho
    vector<OCP_DBL> lupblock_Trans;    ///< last upblock_Trans
    vector<OCP_DBL> lupblock_Velocity; ///< last upblock_Velocity
    vector<OCP_DBL> lAdkt;             ///< last Adkt

    // Derivatives
    vector<OCP_DBL> AdktP; ///< d Adkt / d P, order: connections -> bId.P -> eId.P
    vector<OCP_DBL> AdktT; ///< d Adkt / d T, order: connections -> bId.T -> eId.T
    vector<OCP_DBL>
        AdktS; ///< d Adkt / d S, order: connections -> bId.phase -> eId.phase

    // Last time step
    vector<OCP_DBL> lAdktP; ///< last AdktP
    vector<OCP_DBL> lAdktT; ///< last AdktT
    vector<OCP_DBL> lAdktS; ///< last AdktS

};

class OCPFlux
{
public:
    OCPFlux() = default;
    // Calculate flux of components and phases
    virtual void CalFlux(const BulkPair& bp) = 0;


protected:
    vector<OCP_DBL> flux_Ni;  ///< flux of componnets
    vector<OCP_DBL> flux_nj;  ///< flux of phases

};


using namespace std;



#endif // __OCPFLUX_HEADER__

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           May/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/
