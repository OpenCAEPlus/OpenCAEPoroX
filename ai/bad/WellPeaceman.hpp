/*! \file    WellPeaceman.hpp 
 *  \brief   WellPeacema class declaration 
 *  \author  Shizhe Li 
 *  \date    Aug/17/2023 
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __WELLPEACEMAN_HEADER__
#define __WELLPEACEMAN_HEADER__

// Standard header files
#include <cassert>

// OpenCAEPoroX header files
#include "Well.hpp"
#include "WellPerf.hpp"

using namespace std;

/// Peaceman Well Model class 
PeacemanWell : public Well{
public:
    /*! \brief Input the param of perforations.
     *  \param well The well parameters.
     *  \param domain The domain parameters.
     *  \param wId The well ID.
     */
    void InputPerfo(const WellParam& well, const Domain& domain, const USI& wId) override;

    /*! \brief Setup the well after Grid and Bulk finish setup.
     *  \param bk The bulk parameters.
     *  \param sols The vector of solvent parameters.
     */
    void Setup(const Bulk& bk, const vector<SolventINJ>& sols) override;
    
    /*! \brief Initialize Well Pressure.
     *  \param bk The bulk parameters.
     */
    void InitWellP(const Bulk& bk) override;

    /*! \brief Check if well operation mode would be changed.
     *  \param bk The bulk parameters.
     */
    void CheckOptMode(const Bulk& bk) override;

    /*! \brief Calculate Flux and initialize values which will not be changed during this time step.
     *  \param bk The bulk parameters.
     */
    void CalFluxInit(const Bulk& bk) override;

    /*! \brief Calculate Flux.
     *  \param bk The bulk parameters.
     */
    void CalFlux(const Bulk& bk) override;

    /*! \brief Check if abnormal Pressure occurs.
     *  \param bk The bulk parameters.
     *  \return The reservoir state.
     */
    ReservoirState CheckP(const Bulk& bk) override;

    /*! \brief Calculate flow rate of moles of phases for injection well and production well.
     *  \param bk The bulk parameters.
     *  \param dt The time step.
     */
    void CalIPRate(const Bulk& bk, const OCP_DBL& dt) override;

    /*! \brief Calculate max change of well pressure between two time step.
     *  \return The max change of well pressure.
     */
    OCP_DBL CalMaxChangeTime() const override;

    /*! \brief Calculate max change of well pressure between two NR step.
     *  \return The max change of well pressure.
     */
    OCP_DBL CalMaxChangeNR() override;

    /*! \brief Reset to last time step.
     *  \param bk The bulk parameters.
     */
    void ResetToLastTimeStep(const Bulk& bk) override;

    /*! \brief Update last time step.
     */
    void UpdateLastTimeStep() override;

protected:
    /*! \brief Calculate Well Index with Peaceman model.
     *  \param bk The bulk parameters.
     */
    void CalWI(const Bulk& bk);

    /*! \brief Calculate transmissibility for each phase in perforations.
     *  \param bk The bulk parameters.
     */
    void CalTrans(const Bulk& bk);

    /*! \brief Calculate the flux for each perforations.
     *  \param bk The bulk parameters.
     *  \param ReCalXi The boolean value indicating if Xi needs to be recalculated.
     */
    void CalFlux(const Bulk& bk, const OCP_BOOL ReCalXi);

    /*! \brief Calculate flow rate of moles of phases for injection well with maxBHP.
     *  \param bk The bulk parameters.
     *  \return The flow rate of moles of phases for injection well.
     */
    OCP_DBL CalInjRateMaxBHP(const Bulk& bk);

    /*! \brief Calculate flow rate of moles of phases for production well with minBHP.
     *  \param bk The bulk parameters.
     *  \return The flow rate of moles of phases for production well.
     */
    OCP_DBL CalProdRateMinBHP(const Bulk& bk);

    /*! \brief Calculate flow rate of moles of phases for injection well with calculated qi_lbmol.
     *  \param bk The bulk parameters.
     *  \param dt The time step.
     */
    void CalInjQj(const Bulk& bk, const OCP_DBL& dt);

    /*! \brief Calculate flow rate of moles of phases for production well with calculated qi_lbmol.
     *  \param bk The bulk parameters.
     *  \param dt The time step.
     */
    void CalProdQj(const Bulk& bk, const OCP_DBL& dt);

    /*! \brief Check if cross flow happens.
     *  \param bk The bulk parameters.
     *  \return The reservoir state.
     */
    ReservoirState CheckCrossFlow(const Bulk& bk);

    /*! \brief Calculate the production weight.
     *  \param bk The bulk parameters.
     */
    void CalFactor(const Bulk& bk) const;

    /*! \brief Calculate pressure difference between well and perforations.
     *  \param bk The bulk parameters.
     */
    void CaldG(const Bulk& bk);

    /*! \brief Calculate pressure difference between well and perforations for Injection.
     *  \param bk The bulk parameters.
     */
    void CalInjdG(const Bulk& bk);

    /*! \brief Calculate pressure difference between well and perforations for Production.
     *  \param bk The bulk parameters.
     */
    void CalProddG(const Bulk& bk);

    /*! \brief Calculate pressure difference between well and perforations for Production.
     *  \param bk The bulk parameters.
     */
    void CalProddG01(const Bulk& bk);

    /*! \brief Calculate pressure difference between well and perforations for Production.
     *  \param bk The bulk parameters.
     */
    void CalProddG02(const Bulk& bk);

    /*! \brief Calculate pressure of perforations.
     */
    void CalPerfP() { for (USI p = 0; p < numPerf; p++) perf[p].P = bhp + dG[p]; }

protected:
    /// difference of pressure between well and perforation: numPerf.
    vector<OCP_DBL> dG;

    /// components mole number -> target phase volume
    mutable vector<OCP_DBL> factor;
};

class PeacemanWellIsoT : public PeacemanWell{
public:
    /*! \brief Calculate residual for implicit time stepping.
     *  \param wId The well ID.
     *  \param res The residual parameters.
     *  \param bk The bulk parameters.
     *  \param dt The time step.
     */
    void CalResFIM(OCP_USI& wId, OCPNRresidual& res, const Bulk& bk, const OCP_DBL& dt) const override;

    /*! \brief Get the solution for implicit time stepping.
     *  \param u The vector of solution parameters.
     *  \param wId The well ID.
     */
    void GetSolutionFIM(const vector<OCP_DBL>& u, OCP_USI& wId) override;

    /*! \brief Assemble the matrix for implicit time stepping.
     *  \param ls The linear system parameters.
     *  \param bk The bulk parameters.
     *  \param dt The time step.
     */
    void AssembleMatFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const override;

protected:
    /*! \brief Assemble the matrix for injection well in implicit time stepping.
     *  \param ls The linear system parameters.
     *  \param bk The bulk parameters.
     *  \param dt The time step.
     */
    void AssembleMatInjFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;

    /*! \brief Assemble the matrix for production well in implicit time stepping.
     *  \param ls The linear system parameters.
     *  \param bk The bulk parameters.
     *  \param dt The time step.
     */
    void AssembleMatProdFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;

public:
    /*! \brief Get the solution for IMPEC time stepping.
     *  \param u The vector of solution parameters.
     *  \param wId The well ID.
     */
    void GetSolutionIMPEC(const vector<OCP_DBL>& u, OCP_USI& wId) override;

    /*! \brief Assemble the matrix for IMPEC time stepping.
     *  \param ls The linear system parameters.
     *  \param bk The bulk parameters.
     *  \param dt The time step.
     */
    void AssembleMatIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const override;

protected:
    /*! \brief Assemble the matrix for injection well in IMPEC time stepping.
     *  \param ls The linear system parameters.
     *  \param bk The bulk parameters.
     *  \param dt The time step.
     */
    void AssembleMatInjIMPEC(LinearSystem&

