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
    OCP_DBL    trans{ 0.0 };
    /// Effective area for diffusity
    OCP_DBL    diffu{ 0.0 };
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
    /// Number of connections between bulks.
    OCP_USI         numConn;

public:
    /// For Convection
    /// Index of upwinding bulk of connections for each phase
    vector<OCP_USI> upblock;
    /// Pressure difference between connection bulks for each phase
    vector<OCP_DBL> dP;
    /// Volume flow rate of phase from upblock
    vector<OCP_DBL> flux_vj;
    /// mole flow rate of components 
    vector<OCP_DBL> flux_ni;  

    /// last upblock
    vector<OCP_USI> lupblock;
    /// last dP
    vector<OCP_DBL> ldP;     

};


class FluxVarSet
{
public:
    void Allocate(const USI& _np, const USI& _nc, const USI& _ncol1, const USI& _ncol2) {
        np    = _np;
        nc    = _nc;
        ncol1 = _ncol1;
        ncol2 = _ncol2;

        flux_ni.resize(nc);
        dFdXpB.resize(ncol1 * ncol1);
        dFdXpE.resize(ncol1 * ncol1);
        dFdXsB.resize(ncol1 * ncol2);
        dFdXsE.resize(ncol1 * ncol2);
    }
    void SetZeroFluxNi() {
        fill(flux_ni.begin(), flux_ni.end(), 0.0);
    }
    void SetZeroFIM() {
        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
        fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
        fill(dFdXsE.begin(), dFdXsE.end(), 0.0);
    }
    void SetZeroIMPEC() {
        valbb = 0;
        valee = 0;
        rhsb = 0;
        rhse = 0;
    }

public:
    /// number of phase
    USI              np;
    /// number of components
    USI              nc;

    /// mole flow rate of components 
    vector<OCP_DBL>  flux_ni;

    // for FIM
    USI              ncol1;
    USI              ncol2;
    /// dF / dXp for bId bulk
    vector<OCP_DBL>  dFdXpB;
    /// dF / dXp for eId bulk
    vector<OCP_DBL>  dFdXpE;
    /// dF / dXs for bId bulk
    vector<OCP_DBL>  dFdXsB;
    /// dF / dXs for eId bulk
    vector<OCP_DBL>  dFdXsE;

    // for IMPEC
    /// val in b-b, -val in b-e
    OCP_DBL          valbb;
    /// val in e-e, -val in e-b
    OCP_DBL          valee;
    /// rhs in b
    OCP_DBL          rhsb;
    /// rhs in e
    OCP_DBL          rhse;
};



#endif /* end if __BulkConnVarSet_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/23/2023      Create file                          */
/*----------------------------------------------------------------------------*/
