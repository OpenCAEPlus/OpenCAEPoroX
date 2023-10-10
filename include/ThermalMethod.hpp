/*! \file    ThermalMethod.hpp
 *  \brief   ThermalMethod class declaration
 *  \author  Shizhe Li
 *  \date    Nov/10/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __THERMALMETHOD_HEADER__
#define __THERMALMETHOD_HEADER__

#include "LinearSystem.hpp"
#include "OCPControl.hpp"
#include "Reservoir.hpp"
#include "UtilOutput.hpp"
#include "UtilTiming.hpp"
#include "OCPTimeRecord.hpp"
#include "OCPNRsuite.hpp"

class T_FIM
{
public:
    void     Setup(Reservoir& rs, LinearSystem& ls, const OCPControl& ctrl);
    void     InitReservoir(Reservoir& rs);
    void     Prepare(Reservoir& rs, const OCPControl& ctrl);
    void     AssembleMat(LinearSystem&    ls,
                         const Reservoir& rs,
                         const OCP_DBL&   dt);
    void     SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl);
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);
    OCP_BOOL FinishNR(Reservoir& rs, OCPControl& ctrl);
    void     FinishStep(Reservoir& rs, OCPControl& ctrl);
    const OCPNRsuite& GetNRsuite() const { return NR; }

protected:
    void AllocateReservoir(Reservoir& rs);
    void AllocateLinearSystem(LinearSystem& ls, const Reservoir& rs, const OCPControl& ctrl);
    void InitRock(Bulk& bk) const;
    void CalRock(Bulk& bk) const;
    void InitFlash(Bulk& bk);
    void CalFlash(Bulk& bk);
    void PassFlashValue(Bulk& bk, const OCP_USI& n);
    void CalKrPc(Bulk& bk) const;
    void UpdateLastTimeStep(Reservoir& rs) const;
    void CalRes(Reservoir& rs, const OCP_DBL& dt);
    void AssembleMatBulks(LinearSystem&    ls, const Reservoir& rs, const OCP_DBL&   dt) const;
    void AssembleMatWells(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    void GetSolution(Reservoir& rs, vector<OCP_DBL>& u, const ControlNR& ctrlNR);
    void ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl);

protected:
    /// Newton-Raphson iteration suite
    OCPNRsuite      NR;
};

#endif /* end if __THERMALMETHOD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/