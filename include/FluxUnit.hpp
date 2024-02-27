/*! \file    FluxUnit.hpp
 *  \brief   FluxUnit class declaration
 *  \author  Shizhe Li
 *  \date    Feb/26/2024
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __FLUXUNIT_HEADER__
#define __FLUXUNIT_HEADER__

 // Standard header files
#include <cassert>

// OpenCAEPoroX header files
#include "OCPFlux.hpp"

using namespace std;


class FluxUnit
{
public:
    FluxUnit(const USI& index, const USI& np, const USI& nc) {
        switch (index)
        {
        case 0:
            flux = new OCPFlux01(np, nc);
            break;
        case 1:
            flux = new OCPFlux02(np, nc);
            break;
        case 2:
            flux = new OCPFluxT01(np, nc);
            break;
        default:
            OCP_ABORT("Wrong flux type!");
            break;
        }
    }
    /// Calculate transmissibility
    void CalTrans(BulkConnPair& bp, const Bulk& bk) const {
        flux->CalTrans(bp, bk);
    }
    /// Calculate diffusity for all connections
    void CalDiffu(BulkConnPair& bp, const Bulk& bk) const {
        flux->CalDiffu(bp, bk);
    }
    /// Calculate flux of components and phases
    void CalFlux(const BulkConnPair& bp, const Bulk& bk) const {
        flux->CalFlux(bp, bk);
    }
    /// Assemble matrix for FIM
    void AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) const {
        flux->AssembleMatFIM(bp, c, bcvs, bk);
    }
    /// Assemble matrix for AIM
    void AssembleMatAIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) const {
        flux->AssembleMatAIM(bp, c, bcvs, bk);
    }
    /// Assemble matrix for IMPEC
    void AssembleMatIMPEC(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) const {
        flux->AssembleMatIMPEC(bp, c, bcvs, bk);
    }


    const vector<OCP_USI>& GetUpblock() const { return flux->GetUpblock(); }
    const vector<OCP_DBL>& GetRho() const { return flux->GetRho(); }
    const vector<OCP_DBL>& GetFluxVj() const { return flux->GetFluxVj(); }
    const vector<OCP_DBL>& GetFluxNi() const { return flux->GetFluxNi(); }
    const vector<OCP_DBL>& GetFluxHj() const { return flux->GetFluxHj(); }
    OCP_DBL GetConductH() const { return flux->GetConductH(); }

    const vector<OCP_DBL>& GetdFdXpB() const { return flux->GetdFdXpB(); }
    const vector<OCP_DBL>& GetdFdXpE() const { return flux->GetdFdXpE(); }
    const vector<OCP_DBL>& GetdFdXsB() const { return flux->GetdFdXsB(); }
    const vector<OCP_DBL>& GetdFdXsE() const { return flux->GetdFdXsE(); }
    OCP_DBL GetValbb() const { return flux->GetValbb(); }            
    OCP_DBL GetValee() const { return flux->GetValee(); }
    OCP_DBL GetRhsb() const { return  flux->GetRhsb(); }
    OCP_DBL GetRhse() const { return  flux->GetRhse(); }

protected:

	OCPFlux* flux;

};



#endif /* end if __PVTMODULE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Feb/26/2024      Create file                          */
/*----------------------------------------------------------------------------*/