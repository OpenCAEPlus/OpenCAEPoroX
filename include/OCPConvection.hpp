/*! \file    OCPConvection.hpp
 *  \brief   OCPConvection class declaration
 *  \author  Shizhe Li
 *  \date    May/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPCONVECTION_HEADER__
#define __OCPCONVECTION_HEADER__

 // Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "Bulk.hpp"
#include "BulkConnVarSet.hpp"
#include "BulkConnTrans.hpp"

using namespace std;


class OCPConvection
{

public:
    OCPConvection() = default;
    /// Calculate flux of components and phases
    virtual void CalFlux(const BulkConnPair& bp, const Bulk& bk, FluxVarSet& fvs) = 0;
    /// Assemble matrix for FIM
    virtual void AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk, FluxVarSet& fvs) = 0;
    /// Assemble matrix for AIM
    virtual void AssembleMatAIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk, FluxVarSet& fvs) = 0;
    /// Assemble matrix for IMPEC
    virtual void AssembleMatIMPEC(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk, FluxVarSet& fvs) = 0;

    
    const vector<OCP_USI>& GetUpblock() const { return upblock; }
    const vector<OCP_DBL>& GetDP() const { return dP; }
    const vector<OCP_DBL>& GetVj() const { return vj; }
    const vector<OCP_DBL>& GetHj() const { return Hj; }

protected:
    void Setup(const USI& npin, const USI& ncin) {
        Allocate(npin, ncin);
        bcT.Setup();
        bcD.Setup();
    }
    void Allocate(const USI& npin, const USI& ncin) {
        np = npin;
        nc = ncin;
        upblock.resize(np, 0.0);
        dP.resize(np, 0.0);
        vj.resize(np, 0.0);
    }

protected:

    /// number of phase
    USI              np;
    /// number of components
    USI              nc;
    /// Index of upwinding bulk of connections for each phase
    vector<OCP_USI>  upblock;
    /// Pressure difference between connection bulks for each phase
    vector<OCP_DBL>  dP;
    /// Volume flow rate of phase from upblock
    vector<OCP_DBL>  vj;
    /// enthalpy flow rate of phase from upblock
    vector<OCP_DBL>  Hj;   

public:
    /// Calculate transmissibility
    void CalTrans(BulkConnPair& bp, const Bulk& bk) { bcT.CalTrans(bp, bk); }
    /// Calculate diffusity
    void CalDiffu(BulkConnPair& bp, const Bulk& bk) { bcD.CalDiffu(bp, bk); }
protected:
    /// area calculation of bulk connection
    BulkConnTrans     bcT;
    BulkConnDiffu     bcD;
};


/// For Isothermal darcy flux
class OCPConvection01 : public OCPConvection
{
public:
    OCPConvection01() = default;
    OCPConvection01(const USI& npin, const USI& ncin) {
        Setup(npin, ncin);
    }
    void CalFlux(const BulkConnPair& bp, const Bulk& bk, FluxVarSet& fvs) override;
    void AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk, FluxVarSet& fvs) override;
    void AssembleMatAIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk, FluxVarSet& fvs) override;
    void AssembleMatIMPEC(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk, FluxVarSet& fvs) override;
};


/// For Isothermal darcy flux with standard gravity drainage model
class OCPConvection02 : public OCPConvection
{
public:
    OCPConvection02() = default;
    OCPConvection02(const USI& npin, const USI& ncin) { Setup(npin, ncin); }
    void CalFlux(const BulkConnPair& bp, const Bulk& bk, FluxVarSet& fvs) override;
    void AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk, FluxVarSet& fvs) override;
    void AssembleMatAIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk, FluxVarSet& fvs) override { OCP_ABORT("NOT USED!"); }
    void AssembleMatIMPEC(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk, FluxVarSet& fvs) override { OCP_ABORT("NOT USED!"); }
};


/// For thermal darcy flux
class OCPConvectionT01 : public OCPConvection
{
public:
    OCPConvectionT01() = default;
    OCPConvectionT01(const USI& npin, const USI& ncin) {
        Setup(npin, ncin);
        Hj.resize(np);
    }
    void CalFlux(const BulkConnPair & bp, const Bulk& bk, FluxVarSet& fvs) override;
    void AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk, FluxVarSet& fvs) override;
    void AssembleMatAIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk, FluxVarSet& fvs) override{}
    void AssembleMatIMPEC(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk, FluxVarSet& fvs) override{}
};


#endif // __OCPCONVECTION_HEADER__

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           May/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/
