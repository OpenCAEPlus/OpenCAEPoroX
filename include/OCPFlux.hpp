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

using namespace std;


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
    USI Type() const { return type; }     ///< Get connection type
    OCP_DBL AreaB() const { return areaB; }
    OCP_DBL AreaE() const { return areaE; }
    USI     Direction() const { return direction; }

protected:
    OCP_USI bId;       ///< Beginning index of a pair.
    OCP_USI eId;       ///< Ending index of a pair.
    OCP_DBL area;      ///< Effective area
    USI     type;      ///< connection type
    USI     direction; ///< 1-x, 2-y, 3-z
    OCP_DBL areaB;     ///< Area of intersecting faces from Begin grid
    OCP_DBL areaE;     ///< Area of intersecting faces from End grid   
};


class BulkConnVal
{
    friend class Reservoir;
    friend class OCPFlux_IsoT;
    friend class OCPFlux_T;
    friend class IsoT_FIM;
    friend class IsoT_IMPEC;
    friend class IsoT_AIM;
    friend class T_FIM;

protected:

    //  Note: Upblock is identified by difference of pressure between phases.
    vector<OCP_USI> upblock;  ///< Index of upwinding bulk of connections : numConn * numPhase.
    vector<OCP_DBL> rho;      ///< Mass density of phase from upblock: numConn * numPhase.
    vector<OCP_DBL> velocity; ///< Volume flow rate of phase from upblock: numConn * numsPhase.
    vector<OCP_DBL> flux_ni;  ///< mole velocity of components 
    vector<OCP_DBL> Adkt;     ///< Thermal conductivity between neighbors

    // Last time step
    vector<OCP_USI> lupblock;  ///< last upblock
    vector<OCP_DBL> lrho;      ///< last upblock_Rho
    vector<OCP_DBL> lvelocity; ///< last upblock_Velocity
    vector<OCP_DBL> lflux_ni;  ///< mole velocity of components
    vector<OCP_DBL> lAdkt;     ///< last Adkt

    // Derivatives
    vector<OCP_DBL> AdktP; ///< d Adkt / d P, order: connections -> bId.P -> eId.P
    vector<OCP_DBL> AdktT; ///< d Adkt / d T, order: connections -> bId.T -> eId.T
    vector<OCP_DBL> AdktS; ///< d Adkt / d S, order: connections -> bId.phase -> eId.phase

    // Last time step
    vector<OCP_DBL> lAdktP; ///< last AdktP
    vector<OCP_DBL> lAdktT; ///< last AdktT
    vector<OCP_DBL> lAdktS; ///< last AdktS
};


class OCPFlux
{
    friend class IsoT_FIM;

public:
    OCPFlux() = default;
    void Allocate(const USI& np, const USI& nc) {
        numPhase = np;
        numCom   = nc;
        upblock.resize(np, 0.0);
        rho.resize(np, 0.0);
        flux_vj.resize(np, 0.0);
        flux_ni.resize(nc, 0.0);
    }
    // Calculate flux of components and phases
    virtual void CalFlux(const BulkPair& bp, const Bulk& bk) = 0;
    virtual void AssembleMatFIM(const BulkPair& bp, const OCP_USI& c, const BulkConnVal& bcv, const Bulk& bk) = 0;
    virtual void AssembleMatAIM(const BulkPair& bp, const OCP_USI& c, const BulkConnVal& bcv, const Bulk& bk) = 0;
    virtual void AssembleMatIMPEC(const BulkPair& bp, const OCP_USI& c, const BulkConnVal& bcv, const Bulk& bk) = 0;

    
    const vector<OCP_USI>& GetUpblock() const { return upblock; }
    const vector<OCP_DBL>& GetRho() const { return rho; }
    const vector<OCP_DBL>& GetFluxVj() const { return flux_vj; }
    const vector<OCP_DBL>& GetFluxNi() const { return flux_ni; }
    OCP_DBL GetAdkt() const { return Adkt; }

    const vector<OCP_DBL>& GetdFdXpB() const { return dFdXpB; }
    const vector<OCP_DBL>& GetdFdXpE() const { return dFdXpE; }
    const vector<OCP_DBL>& GetdFdXsB() const { return dFdXsB; }
    const vector<OCP_DBL>& GetdFdXsE() const { return dFdXsE; }
    OCP_DBL GetValbb() const { return valbb; }
    OCP_DBL GetValee() const { return valee; }
    OCP_DBL GetRhsb() const { return  rhsb; }
    OCP_DBL GetRhse() const { return  rhse; }

protected:

    // Common variables 
    USI              numPhase; ///< num of phase
    USI              numCom;   ///< num of components
    vector<OCP_USI>  upblock;  ///< upblock of connections
    vector<OCP_DBL>  rho;      ///< weighted density of phase between bulks
    vector<OCP_DBL>  flux_vj;  ///< volume velocity of phase from upblock
    vector<OCP_DBL>  flux_ni;  ///< mole velocity of components

    // for thermal
    OCP_DBL          Adkt;     ///< efficient thermal conduction coefficient

    // for FIM
    vector<OCP_DBL>  dFdXpB;   ///< dF / dXp for bId bulk
    vector<OCP_DBL>  dFdXpE;   ///< dF / dXp for eId bulk
    vector<OCP_DBL>  dFdXsB;   ///< dF / dXs for bId bulk
    vector<OCP_DBL>  dFdXsE;   ///< dF / dXs for eId bulk

    // for IMPEC
    OCP_DBL          valbb;    ///< val in b-b, -val in b-e
    OCP_DBL          valee;    ///< val in e-e, -val in e-b
    OCP_DBL          rhsb;     ///< rhs in b
    OCP_DBL          rhse;     ///< rhs in e
};


class OCPFlux_IsoT : public OCPFlux
{
public:
    OCPFlux_IsoT() = default;
    OCPFlux_IsoT(const Bulk& bk) {
        Allocate(bk.numPhase, bk.numCom);
        dFdXpB.resize((numCom + 1) * (numCom + 1));
        dFdXpE.resize((numCom + 1) * (numCom + 1));
        dFdXsB.resize((numCom + 1) * (numCom + 1) * numPhase);
        dFdXsE.resize((numCom + 1) * (numCom + 1) * numPhase);
    }
    void CalFlux(const BulkPair& bp, const Bulk& bk) override;
    void AssembleMatFIM(const BulkPair& bp, const OCP_USI& c, const BulkConnVal& bcv, const Bulk& bk) override;
    void AssembleMatAIM(const BulkPair& bp, const OCP_USI& c, const BulkConnVal& bcv, const Bulk& bk) override;
    void AssembleMatIMPEC(const BulkPair& bp, const OCP_USI& c, const BulkConnVal& bcv, const Bulk& bk) override;
};


class OCPFlux_T : public OCPFlux
{
    // flux for thermal model, flux between connetions with fluid flow
    // fluid flow or heat conduction
public:
    OCPFlux_T() = default;
    OCPFlux_T(const Bulk& bk) {
        Allocate(bk.numPhase, bk.numCom);
        dFdXpB.resize((numCom + 2) * (numCom + 2));
        dFdXpE.resize((numCom + 2) * (numCom + 2));
        dFdXsB.resize((numCom + 2) * (numCom + 1) * numPhase);
        dFdXsE.resize((numCom + 2) * (numCom + 1) * numPhase);
    }
    void CalFlux(const BulkPair & bp, const Bulk & bk) override;
    void AssembleMatFIM(const BulkPair & bp, const OCP_USI& c, const BulkConnVal& bcv, const Bulk & bk) override;
    void AssembleMatAIM(const BulkPair & bp, const OCP_USI& c, const BulkConnVal& bcv, const Bulk & bk) override{}
    void AssembleMatIMPEC(const BulkPair & bp, const OCP_USI& c, const BulkConnVal& bcv, const Bulk & bk) override{}
};


#endif // __OCPFLUX_HEADER__

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           May/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/
