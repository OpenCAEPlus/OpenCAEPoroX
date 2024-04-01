/*! \file    IsothermalSolver.hpp
 *  \brief   IsothermalSolver class declaration
 *  \author  Shizhe Li
 *  \date    Oct/21/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __ISOTHERMALSOLVER_HEADER__
#define __ISOTHERMALSOLVER_HEADER__

// OpenCAEPoroX header files
#include "IsoThermalMethod.hpp"

/// IsothermalSolver class for fluid solution method.
class IsothermalSolver
{

public:
    /// Setup the isothermal solver.
    void SetupMethod(Reservoir& rs, const OCPControl& ctrl);
    /// Initialize the Reservoir.
    void InitReservoir(Reservoir& rs);
    /// calculate one time step.
    const OCPNRsuite& GoOneStep(Reservoir& rs, OCPControl& ctrl);

protected:
    /// Prepare for assembling Mat.
    void Prepare(Reservoir& rs, OCPControl& ctrl);
    /// Assemble jacobian matrix.
    void AssembleMat(const Reservoir& rs, OCPControl& ctrl);
    /// Solve the linear system in single problem.
    OCP_BOOL SolveLinearSystem(Reservoir& rs, OCPControl& ctrl);
    /// Update properties of fluid.
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);
    /// Finish the Newton-Raphson iteration.
    OCP_BOOL FinishNR(Reservoir& rs, OCPControl& ctrl);
    /// Finish the current time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl);
    /// GetNRsuite
    const OCPNRsuite& GetNRsuite() const;

private:
    /// current method
    OCPNLMethod  curMethod{ OCPNLMethod::none };
    OCPNLMethod  mainMethod{ OCPNLMethod::none };
    LinearSystem LSolver;
    IsoT_IMPEC   impec;
    IsoT_FIM     fim;
    IsoT_AIMc    aimc;
    IsoT_FIMddm  fim_ddm;
};

#endif /* end if __ISOTHERMALSOLVER_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/21/2021      Create file                          */
/*  Chensong Zhang      Jan/16/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/