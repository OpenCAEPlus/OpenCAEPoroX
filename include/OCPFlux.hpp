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


class BulkConnValSet
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

    // Last time step
    vector<OCP_USI> lupblock;  ///< last upblock
    vector<OCP_DBL> lrho;      ///< last upblock_Rho
};


class OCPFlux
{
    friend class IsoT_FIM;

public:
    OCPFlux() = default;
    void Allocate(const USI& npin, const USI& ncin) {
        np = npin;
        nc = ncin;
        upblock.resize(np, 0.0);
        rho.resize(np, 0.0);
        flux_vj.resize(np, 0.0);
        flux_ni.resize(nc, 0.0);
    }
    // Calculate flux of components and phases
    virtual void CalFlux(const BulkPair& bp, const Bulk& bk) = 0;
    virtual void AssembleMatFIM(const BulkPair& bp, const OCP_USI& c, const BulkConnValSet& bcv, const Bulk& bk) = 0;
    virtual void AssembleMatAIM(const BulkPair& bp, const OCP_USI& c, const BulkConnValSet& bcv, const Bulk& bk) = 0;
    virtual void AssembleMatIMPEC(const BulkPair& bp, const OCP_USI& c, const BulkConnValSet& bcv, const Bulk& bk) = 0;

    
    const vector<OCP_USI>& GetUpblock() const { return upblock; }
    const vector<OCP_DBL>& GetRho() const { return rho; }
    const vector<OCP_DBL>& GetFluxVj() const { return flux_vj; }
    const vector<OCP_DBL>& GetFluxNi() const { return flux_ni; }
    const vector<OCP_DBL>& GetFluxHj() const { return flux_Hj; }
    OCP_DBL GetConductH() const { return conduct_H; }

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
    USI              np; ///< num of phase
    USI              nc;   ///< num of components
    vector<OCP_USI>  upblock;  ///< upblock of connections
    vector<OCP_DBL>  rho;      ///< weighted density of phase between bulks
    vector<OCP_DBL>  flux_vj;  ///< volume velocity of phase from upblock
    vector<OCP_DBL>  flux_ni;  ///< mole velocity of components
    vector<OCP_DBL>  flux_Hj;  ///< enthalpy velocity of phase from upblock

    // for thermal
    OCP_DBL          conduct_H;     ///< efficient thermal conduction coefficient

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
        Allocate(bk.vs.np, bk.vs.nc);
        dFdXpB.resize((nc + 1) * (nc + 1));
        dFdXpE.resize((nc + 1) * (nc + 1));
        dFdXsB.resize((nc + 1) * (nc + 1) * np);
        dFdXsE.resize((nc + 1) * (nc + 1) * np);
    }
    void CalFlux(const BulkPair& bp, const Bulk& bk) override;
    void AssembleMatFIM(const BulkPair& bp, const OCP_USI& c, const BulkConnValSet& bcv, const Bulk& bk) override;
    void AssembleMatAIM(const BulkPair& bp, const OCP_USI& c, const BulkConnValSet& bcv, const Bulk& bk) override;
    void AssembleMatIMPEC(const BulkPair& bp, const OCP_USI& c, const BulkConnValSet& bcv, const Bulk& bk) override;
};


class OCPFlux_T : public OCPFlux
{
    // flux for thermal model,
    // fluid flow or heat conduction
public:
    OCPFlux_T() = default;
    OCPFlux_T(const Bulk& bk) {
        Allocate(bk.vs.np, bk.vs.nc);
        flux_Hj.resize(np);
        dFdXpB.resize((nc + 2) * (nc + 2));
        dFdXpE.resize((nc + 2) * (nc + 2));
        dFdXsB.resize((nc + 2) * (nc + 1) * np);
        dFdXsE.resize((nc + 2) * (nc + 1) * np);
    }
    void CalFlux(const BulkPair & bp, const Bulk& bk) override;
    void AssembleMatFIM(const BulkPair& bp, const OCP_USI& c, const BulkConnValSet& bcv, const Bulk& bk) override;
    void AssembleMatAIM(const BulkPair& bp, const OCP_USI& c, const BulkConnValSet& bcv, const Bulk& bk) override{}
    void AssembleMatIMPEC(const BulkPair& bp, const OCP_USI& c, const BulkConnValSet& bcv, const Bulk& bk) override{}
};


#endif // __OCPFLUX_HEADER__

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           May/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/
