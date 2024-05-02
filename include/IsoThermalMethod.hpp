/*! \file    IsothermalMethod.hpp
 *  \brief   Declaration of solution methods for fluid part in OpenCAEPoroX
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __ISOTHERMALMETHOD_HEADER__
#define __ISOTHERMALMETHOD_HEADER__

// OpenCAEPoroX header files
#include "LinearSystem.hpp"
#include "OCPControl.hpp"
#include "Reservoir.hpp"
#include "UtilOutput.hpp"
#include "UtilTiming.hpp"
#include "OCPTimeRecord.hpp"
#include "OCPNRsuite.hpp"

class IsothermalMethod
{
public:
    /// Calculate rock
    void CalRock(Bulk& bk) const;
    /// Get NRsuite
    const OCPNRsuite& GetNRsuite() const { return NR; }
    virtual void ExchangeSolutionP(Reservoir& rs) const;
    virtual void ExchangeSolutionNi(Reservoir& rs) const;
    void SetPreMethod(const OCP_BOOL& flag) { preM = flag; }
    void SetWorkLS(const USI& w, const USI& i);
    USI  GetWorkLS()const { return wls; }

protected:
    /// If use as a preconditioner for other method
    OCP_BOOL        preM{ OCP_FALSE };
    /// Newton-Raphson iteration suite
    OCPNRsuite      NR;
    /// Index of linear solver method
    USI             wls;
};

/// IsoT_IMPEC is IMPEC (implicit pressure explict saturation) method.
class IsoT_IMPEC : virtual public IsothermalMethod
{
public:
    /// Setup IMPEC
    void Setup(Reservoir& rs, const OCPControl& ctrl);
    /// Init
    void InitReservoir(Reservoir& rs) const;
    /// Prepare for Assembling matrix.
    void Prepare(Reservoir& rs, OCPControl& ctrl);
    /// Assemble Matrix
    void AssembleMat(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Solve the linear system.
    OCP_BOOL SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl);
    /// Update properties of fluids.
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);
    /// Determine if NR iteration finishes.
    OCP_BOOL FinishNR(const Reservoir& rs);
    void     FinishStep(Reservoir& rs, OCPControl& ctrl);

protected:
    /// Perform Flash with Sj and calculate values needed for FIM
    void InitFlash(Bulk& bk) const;
    /// Calculate relative permeability and capillary pressure needed for FIM
    void CalKrPc(Bulk& bk) const;
    /// Pass value needed for FIM from flash to bulk
    void PassFlashValue(Bulk& bk, const OCP_USI& n) const;

private:
    /// Allocate memory for reservoir
    void AllocateReservoir(Reservoir& rs);
    /// Perform Flash with Ni and calculate values needed for FIM
    void CalFlash(Bulk& bk);
    /// Calculate flux between bulks and wells
    void CalFlux(Reservoir& rs) const;
    /// Calculate flux between bulks
    void CalBulkFlux(Reservoir& rs) const;
    /// Update mole composition of each bulk according to mass conservation for IMPEC
    void MassConserve(Reservoir& rs, const OCP_DBL& dt) const;
    /// Assemble linear system for bulks
    void AssembleMatBulks(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Assemble linear system for wells
    void AssembleMatWells(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Update P, Pj, BHP after linear system is solved
    void GetSolution(Reservoir& rs, vector<OCP_DBL>& u);
    /// Reset variables to last time step
    void ResetToLastTimeStep01(Reservoir& rs, OCPControl& ctrl);
    void ResetToLastTimeStep02(Reservoir& rs, OCPControl& ctrl);
    /// Update values of last step for FIM.
    void UpdateLastTimeStep(Reservoir& rs) const;
};

/// IsoT_FIM is FIM (Fully Implicit Method).
class IsoT_FIM : virtual public IsothermalMethod
{
public:
    /// Setup FIM
    void Setup(Reservoir& rs, const OCPControl& ctrl);
    /// Init
    void InitReservoir(Reservoir& rs);
    /// Prepare for Assembling matrix.
    void Prepare(Reservoir& rs, const OCP_DBL& dt);
    /// Assemble Matrix
    void AssembleMat(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Solve the linear system.
    OCP_BOOL SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl);
    /// Update properties of fluids.
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);
    /// Finish a Newton-Raphson iteration.
    OCP_BOOL FinishNR(Reservoir& rs, OCPControl& ctrl);
    /// Finish a time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl);
    /// Transfer to FIM method
    void TransferToFIM(Reservoir& rs, OCPControl& ctrl);


protected:
    /// Calculate initial residual
    virtual void CalInitRes(Reservoir& rs, const OCP_DBL& dt) { CalRes(rs, dt, OCP_TRUE); }
    /// Calculate residual
    virtual void CalRes(Reservoir& rs, const OCP_DBL& dt) { CalRes(rs, dt, OCP_FALSE); }
    /// Allocate memory for reservoir
    void AllocateReservoir(Reservoir& rs);
    /// Pass value needed for FIM from flash to bulk
    void PassFlashValue(Bulk& bk, const OCP_USI& n);
    /// Calculate relative permeability and capillary pressure needed for FIM
    void CalKrPc(Bulk& bk) const;
    /// Calculate residual
    void CalRes(Reservoir& rs, const OCP_DBL& dt, const OCP_BOOL& initRes0);
    /// Assemble linear system for wells
    void AssembleMatWells(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Reset variables to last time step
    void ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl);
    /// Update values of last step for FIM.
    void UpdateLastTimeStep(Reservoir& rs) const;

private:
    /// Perform Flash with Sj and calculate values needed for FIM
    void InitFlash(Bulk& bk);
    /// Perform Flash with Ni and calculate values needed for FIM
    void CalFlash(Bulk& bk);
    /// Assemble linear system for bulks
    void AssembleMatBulks(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Update P, Ni, BHP after linear system is solved
    void GetSolution(Reservoir& rs, vector<OCP_DBL>& u, const ControlNR& ctrlNR);
};


class IsoT_AIMc : public IsoT_IMPEC, public IsoT_FIM
{
public:
    /// Setup AIMc
    void Setup(Reservoir& rs, const OCPControl& ctrl);
    /// Setup neighbor
    void SetupNeighbor(Reservoir& rs);
    /// Init
    void InitReservoir(Reservoir& rs);
    /// Prepare for Assembling matrix.
    void Prepare(Reservoir& rs, const OCP_DBL& dt);
    /// Assemble Matrix
    void AssembleMat(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Solve the linear system.
    OCP_BOOL SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl);
    /// Update properties of fluids.
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);
    /// Finish a Newton-Raphson iteration.
    OCP_BOOL FinishNR(Reservoir& rs, OCPControl& ctrl);
    /// Finish a time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl);
    /// Get NRsuite
    const OCPNRsuite& GetNRsuite() const { return NR; }

protected:
    /// Allocate memory for reservoir
    void AllocateReservoir(Reservoir& rs);
    /// Determine which bulk are treated Implicit
    void SetFIMBulk(Reservoir& rs);
    // Set K-neighbors
    void SetKNeighbor(const vector<vector<OCP_USI>>& neighbor, const OCP_USI& p, BulkTypeAIM& tar, OCP_INT k);
    /// Perform flash calculation with Ni for Explicit bulk -- Update partial properties
    void CalFlashEp(Bulk& bk);
    /// Perform flash calculation with Ni for Explicit bulk -- Update all properties
    void CalFlashEa(Bulk& bk);
    /// Perform flash calculation with Ni for Implicit bulk
    void CalFlashI(Bulk& bk);
    /// Pass flash value needed for Explicit bulk -- Update partial properties
    void PassFlashValueEp(Bulk& bk, const OCP_USI& n);
    /// Calculate relative permeability and capillary pressure for Explicit bulk
    void CalKrPcE(Bulk& bk);
    /// Calculate relative permeability and capillary pressure for Implicit bulk
    void CalKrPcI(Bulk& bk);
    /// Assemble linear system for bulks
    void AssembleMatBulks(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Update P, Ni, BHP after linear system is solved
    void GetSolution(Reservoir& rs, vector<OCP_DBL>& u, const ControlNR& ctrlNR);
    /// Reset variables to last time step
    void ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl);
    /// Update values of last step for AIMc.
    void UpdateLastTimeStep(Reservoir& rs) const;
};


class IsoT_FIMddm : public IsoT_FIM
{
public:
    /// Setup FIMddm
    void Setup(Reservoir& rs, const OCPControl& ctrl);
    /// Init
    void InitReservoir(Reservoir& rs);
    /// Prepare for Assembling matrix.
    void Prepare(Reservoir& rs, const OCP_DBL& dt);
    /// Assemble Matrix
    void AssembleMat(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Solve the linear system.
    OCP_BOOL SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl);
    /// Update properties of fluids.
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);
    /// Finish a Newton-Raphson iteration.
    OCP_BOOL FinishNR(Reservoir& rs, OCPControl& ctrl);
    /// Finish a time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl);

protected:
    /// Allocate memory for reservoir
    void AllocateReservoir(Reservoir& rs);
    /// Calculate initial residual
    void CalInitRes(Reservoir& rs, const OCP_DBL& dt) override { CalRes(rs, dt, OCP_TRUE); }
    /// Calculate residual
    void CalRes(Reservoir& rs, const OCP_DBL& dt) override { CalRes(rs, dt, OCP_FALSE); }
    void CalRankSet(const Domain& domain);
    void CalFlash(Bulk& bk, const set<OCP_INT>& rankSet, const Domain& domain);
    void CalKrPc(Bulk& bk, const set<OCP_INT>& rankSet, const Domain& domain);
    void CalRock(Bulk& bk, const set<OCP_INT>& rankSet, const Domain& domain);
    void CalRes(Reservoir& rs, const OCP_DBL& dt, const OCP_BOOL& initRes0);
    void CalResConstP(Reservoir& rs, const OCP_DBL& dt, const OCP_BOOL& initRes0);
    void CalResConstV(Reservoir& rs, const OCP_DBL& dt, const OCP_BOOL& initRes0);
    

protected:
    /// Assemble linear system for bulks
    void AssembleMatBulks(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    void AssembleMatBulksConstP(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    void AssembleMatBulksConstV(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Update P, Ni, BHP after linear system is solved
    void GetSolution(Reservoir& rs, vector<OCP_DBL>& u, const ControlNR& ctrlNR);
    /// Update property for ghost grid
    void UpdatePropertyBoundary(Reservoir& rs);
    void ExchangePBoundary(Reservoir& rs) const;
    void ExchangeNiBoundary(Reservoir& rs) const;
    OCP_BOOL IfBulkInLS(const USI& bId, const Domain& domain) const;
    void SetStarBulkSet(const Bulk& bulk, const Domain& domain);
    void ResetBoundary(Reservoir& rs);
    /// Reset variables to last time step
    void ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl);
    /// Update values of last step for AIMc.
    void UpdateLastTimeStep(Reservoir& rs) const;

protected:
    set<OCP_INT>    rankSetInLS;
    set<OCP_INT>    rankSetOutLS;
    vector<OCP_USI> starBulkSet{};
    /// global initial residual
    OCP_DBL         global_res0;
    /// constant pressure for boundary
    const USI       constP = 0;
    /// constant velocity for boundary
    const USI       constV = 1;
    USI             boundCondition{ constP };
    OCP_DBL         dSlim = 1E-2;
};



#endif /* end if __ISOTHERMALMETHOD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/