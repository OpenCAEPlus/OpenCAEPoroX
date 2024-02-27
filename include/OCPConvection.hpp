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
    virtual void CalFlux(const BulkConnPair& bp, const Bulk& bk) = 0;
    /// Assemble matrix for FIM
    virtual void AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) = 0;
    /// Assemble matrix for AIM
    virtual void AssembleMatAIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) = 0;
    /// Assemble matrix for IMPEC
    virtual void AssembleMatIMPEC(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) = 0;

    
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
    void Setup(const USI& npin, const USI& ncin) {
        Allocate(npin, ncin);
        bcT.Setup();
        bcD.Setup();
    }
    void Allocate(const USI& npin, const USI& ncin) {
        np = npin;
        nc = ncin;
        upblock.resize(np, 0.0);
        rho.resize(np, 0.0);
        flux_vj.resize(np, 0.0);
        flux_ni.resize(nc, 0.0);
    }

protected:

    /// number of phase
    USI              np;
    /// number of components
    USI              nc;
    /// Index of upwinding bulk of connections for each phase
    vector<OCP_USI>  upblock;
    /// Mass density of phase for connections
    vector<OCP_DBL>  rho;
    /// Volume flow rate of phase from upblock
    vector<OCP_DBL>  flux_vj;
    /// mole flow rate of components 
    vector<OCP_DBL>  flux_ni;
    /// enthalpy flow rate of phase from upblock
    vector<OCP_DBL>  flux_Hj;  

    // thermal conduction term
    OCP_DBL          conduct_H;

    // for FIM
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
        dFdXpB.resize((nc + 1) * (nc + 1));
        dFdXpE.resize((nc + 1) * (nc + 1));
        dFdXsB.resize((nc + 1) * (nc + 1) * np);
        dFdXsE.resize((nc + 1) * (nc + 1) * np);
    }
    void CalFlux(const BulkConnPair& bp, const Bulk& bk) override;
    void AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override;
    void AssembleMatAIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override;
    void AssembleMatIMPEC(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override;
};


/// For Isothermal darcy flux with standard gravity drainage model
class OCPConvection02 : public OCPConvection
{
public:
    OCPConvection02() = default;
    OCPConvection02(const USI& npin, const USI& ncin) {
        Setup(npin, ncin);
        dFdXpB.resize((nc + 1) * (nc + 1));
        dFdXpE.resize((nc + 1) * (nc + 1));
        dFdXsB.resize((nc + 1) * (nc + 1) * np);
        dFdXsE.resize((nc + 1) * (nc + 1) * np);
    }
    void CalFlux(const BulkConnPair& bp, const Bulk& bk) override;
    void AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override;
    void AssembleMatAIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override { OCP_ABORT("NOT USED!"); }
    void AssembleMatIMPEC(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override { OCP_ABORT("NOT USED!"); }
};


/// For thermal darcy flux
class OCPConvectionT01 : public OCPConvection
{
public:
    OCPConvectionT01() = default;
    OCPConvectionT01(const USI& npin, const USI& ncin) {
        Setup(npin, ncin);
        flux_Hj.resize(np);
        dFdXpB.resize((nc + 2) * (nc + 2));
        dFdXpE.resize((nc + 2) * (nc + 2));
        dFdXsB.resize((nc + 2) * (nc + 1) * np);
        dFdXsE.resize((nc + 2) * (nc + 1) * np);
    }
    void CalFlux(const BulkConnPair & bp, const Bulk& bk) override;
    void AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override;
    void AssembleMatAIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override{}
    void AssembleMatIMPEC(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override{}
};


#endif // __OCPFLUX_HEADER__

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           May/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/
