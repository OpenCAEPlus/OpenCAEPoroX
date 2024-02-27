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
#include "OCPConvection.hpp"

using namespace std;


class FluxUnit
{
public:
    FluxUnit(const USI& index, const USI& np, const USI& nc) {
        switch (index)
        {
        case 0:
            convect = new OCPConvection01(np, nc);
            break;
        case 1:
            convect = new OCPConvection02(np, nc);
            break;
        case 2:
            convect = new OCPConvectionT01(np, nc);
            break;
        default:
            OCP_ABORT("Wrong flux type!");
            break;
        }
    }
    /// Calculate transmissibility
    void CalTrans(BulkConnPair& bp, const Bulk& bk) const {
        convect->CalTrans(bp, bk);
    }
    /// Calculate diffusity for all connections
    void CalDiffu(BulkConnPair& bp, const Bulk& bk) const {
        convect->CalDiffu(bp, bk);
    }
    /// Calculate flux of components and phases
    void CalFlux(const BulkConnPair& bp, const Bulk& bk) const {
        convect->CalFlux(bp, bk);
    }
    /// Assemble matrix for FIM
    void AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) const {
        convect->AssembleMatFIM(bp, c, bcvs, bk);
    }
    /// Assemble matrix for AIM
    void AssembleMatAIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) const {
        convect->AssembleMatAIM(bp, c, bcvs, bk);
    }
    /// Assemble matrix for IMPEC
    void AssembleMatIMPEC(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) const {
        convect->AssembleMatIMPEC(bp, c, bcvs, bk);
    }


    const vector<OCP_USI>& GetUpblock() const { return convect->GetUpblock(); }
    const vector<OCP_DBL>& GetRho() const { return convect->GetRho(); }
    const vector<OCP_DBL>& GetFluxVj() const { return convect->GetFluxVj(); }
    const vector<OCP_DBL>& GetFluxNi() const { return convect->GetFluxNi(); }
    const vector<OCP_DBL>& GetFluxHj() const { return convect->GetFluxHj(); }
    OCP_DBL GetConductH() const { return convect->GetConductH(); }

    const vector<OCP_DBL>& GetdFdXpB() const { return convect->GetdFdXpB(); }
    const vector<OCP_DBL>& GetdFdXpE() const { return convect->GetdFdXpE(); }
    const vector<OCP_DBL>& GetdFdXsB() const { return convect->GetdFdXsB(); }
    const vector<OCP_DBL>& GetdFdXsE() const { return convect->GetdFdXsE(); }
    OCP_DBL GetValbb() const { return convect->GetValbb(); }            
    OCP_DBL GetValee() const { return convect->GetValee(); }
    OCP_DBL GetRhsb() const { return  convect->GetRhsb(); }
    OCP_DBL GetRhse() const { return  convect->GetRhse(); }

protected:

	OCPConvection* convect;

};



#endif /* end if __PVTMODULE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Feb/26/2024      Create file                          */
/*----------------------------------------------------------------------------*/