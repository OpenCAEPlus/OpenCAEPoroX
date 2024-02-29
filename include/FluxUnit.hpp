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
#include "BulkConnOptionalModules.hpp"

using namespace std;


class FluxUnit
{
public:
    FluxUnit(const USI& index, const USI& np, const USI& nc, BulkConnOptionalModules& optMs) {
        switch (index)
        {
        case 0:
            fluxvs.Allocate(np, nc, nc + 1, (nc + 1) * np);
            convect = new OCPConvection01(np, nc);
            break;
        case 1:
            fluxvs.Allocate(np, nc, nc + 1, (nc + 1) * np);
            convect = new OCPConvection02(np, nc);
            break;
        case 2:
            fluxvs.Allocate(np, nc, nc + 2, (nc + 1) * np);
            convect     = new OCPConvectionT01(np, nc);           
            break;
        default:
            OCP_ABORT("Wrong flux type!");
            break;
        }
        diffusion   = &optMs.diffusion;
        heatConduct = &optMs.heatConduct;
    }
    /// Calculate flux coefficients
    void CalFluxCoeff(BulkConnPair& bp, const Bulk& bk) const {
        convect->CalTrans(bp, bk);
        diffusion->CalDiffu(bp, bk);
    }
    /// Calculate transmissibility
    void CalTrans(BulkConnPair& bp, const Bulk& bk) const {
        convect->CalTrans(bp, bk);
    }
    /// Calculate diffusity for all connections
    void CalDiffu(BulkConnPair& bp, const Bulk& bk) const {
        diffusion->CalDiffu(bp, bk);
    }
    /// Calculate flux of components and phases
    void CalFlux(const BulkConnPair& bp, const Bulk& bk) const {
        fluxvs.SetZeroFluxNi();
        convect->CalFlux(bp, bk, fluxvs);
        diffusion->CalFlux(bp, bk.GetVarSet(), fluxvs);
        heatConduct->CalFlux(bp, bk.GetVarSet());
    }
    /// Assemble matrix for FIM
    void AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) const {
        fluxvs.SetZeroFIM();
        convect->AssembleMatFIM(bp, c, bcvs, bk, fluxvs);
        diffusion->AssembleMatFIM(bp, bk.GetVarSet(), fluxvs);
        heatConduct->AssembleMatFIM(bp, bk.GetVarSet(), fluxvs);
    }
    /// Assemble matrix for AIM
    void AssembleMatAIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) const {
        fluxvs.SetZeroFIM();
        convect->AssembleMatAIM(bp, c, bcvs, bk, fluxvs);
    }
    /// Assemble matrix for IMPEC
    void AssembleMatIMPEC(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) const {
        fluxvs.SetZeroIMPEC();
        convect->AssembleMatIMPEC(bp, c, bcvs, bk, fluxvs);
    }


    const vector<OCP_USI>& GetConvectUpblock() const { return convect->GetUpblock(); }
    const vector<OCP_DBL>& GetConvectDP() const { return convect->GetDP(); }
    const vector<OCP_DBL>& GetConvectVj() const { return convect->GetVj(); }
    const vector<OCP_DBL>& GetConvectHj() const { return convect->GetHj(); }


    OCP_DBL GetConductH() const { return heatConduct->GetConductH(); }

    const vector<OCP_DBL>& GetFluxNi() const { return fluxvs.flux_ni; }
    const vector<OCP_DBL>& GetdFdXpB() const { return fluxvs.dFdXpB; }
    const vector<OCP_DBL>& GetdFdXpE() const { return fluxvs.dFdXpE; }
    const vector<OCP_DBL>& GetdFdXsB() const { return fluxvs.dFdXsB; }
    const vector<OCP_DBL>& GetdFdXsE() const { return fluxvs.dFdXsE; }
    OCP_DBL GetValbb() const { return fluxvs.valbb; }            
    OCP_DBL GetValee() const { return fluxvs.valee; }
    OCP_DBL GetRhsb() const { return  fluxvs.rhsb; }
    OCP_DBL GetRhse() const { return  fluxvs.rhse; }

protected:
    mutable FluxVarSet     fluxvs;
	OCPConvection*         convect;
    OCPDiffusion*          diffusion;
    HeatConduct*           heatConduct; 
};



#endif /* end if __FLUXUNIT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Feb/26/2024      Create file                          */
/*----------------------------------------------------------------------------*/