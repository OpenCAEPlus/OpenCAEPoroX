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

/// Peaceman Well Model 
class PeacemanWell : public Well
{

public:
	/// Input the param of perforations.
	void InputPerfo(const WellParam& well, const Domain& domain, const USI& wId) override;
	/// Setup the well after Grid and Bulk finish setup.
	void Setup(const Bulk& bk, const vector<SolventINJ>& sols) override;

public:
	/// Initialize Well Pressure
	void InitWellP(const Bulk& bk) override;
	/// Check if well operation mode would be changed.
	void CheckOptMode(const Bulk& bk) override;
	/// Calculate Flux and initialize values which will not be changed during this time step
	void CalFluxInit(const Bulk& bk) override;
	/// Calculate Flux
	void CalFlux(const Bulk& bk) override;
	/// Check if abnormal Pressure occurs.
	OCP_INT CheckP(const Bulk& bk) override;
	/// Calculate flow rate of moles of phases for injection well and production well
	void CalIPRate(const Bulk& bk, const OCP_DBL& dt) override;
	/// Reset to last time step
	void ResetToLastTimeStep(const Bulk& bk) override;
	/// Update last time step
	void UpdateLastTimeStep() override;

protected:
	/// Calculate Well Index with Peaceman model.
	void CalWI(const Bulk& bk);
	/// Calculate transmissibility for each phase in perforations.
	void CalTrans(const Bulk& bk);
	/// Calculate the flux for each perforations.
	void CalFlux(const Bulk& bk, const OCP_BOOL ReCalXi);
	/// calculate flow rate of moles of phases for injection well with maxBHP.
	OCP_DBL CalInjRateMaxBHP(const Bulk& bk);
	/// calculate flow rate of moles of phases for production well with minBHP.
	OCP_DBL CalProdRateMinBHP(const Bulk& bk);
	/// Calculate flow rate of moles of phases for injection well with calculated qi_lbmol.
	void CalInjQj(const Bulk& bk, const OCP_DBL& dt);
	/// Calculate flow rate of moles of phases for production well with calculated qi_lbmol.
	void CalProdQj(const Bulk& bk, const OCP_DBL& dt);
	/// Check if cross flow happens.
	OCP_INT CheckCrossFlow(const Bulk& bk);
	/// Calculate the production weight
	void CalFactor(const Bulk& bk) const;
	/// Calculate pressure difference between well and perforations.
	void CaldG(const Bulk& bk);
	/// Calculate pressure difference between well and perforations for Injection.
	void CalInjdG(const Bulk& bk);
	/// Calculate pressure difference between well and perforations for Production.
	void CalProddG(const Bulk& bk);
	/// Calculate pressure difference between well and perforations for Production.
	void CalProddG01(const Bulk& bk);
	/// Calculate pressure difference between well and perforations for Production.
	void CalProddG02(const Bulk& bk);
	/// Calculate pressure of perforations
	void CalPerfP() { for (USI p = 0; p < numPerf; p++) perf[p].P = bhp + dG[p]; }

public:
	void GetSolutionFIM(const OCP_DBL* u) override;
	void GetSolutionIMPEC(const OCP_DBL* u) override; 

protected:
	/// difference of pressure between well and perforation: numPerf.
	vector<OCP_DBL>         dG;
	/// components mole number -> target phase volume
	mutable vector<OCP_DBL> factor;
};


class PeacemanWellIsoT : public PeacemanWell
{
public:
	void CalResFIM(OCP_USI& wId, OCPRes& res, const Bulk& bk, const OCP_DBL& dt) const override;
	void AssembleMatFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const override;
protected:
	void AssembleMatInjFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
	void AssembleMatProdFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;

public:
	void AssembleMatIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const override;
protected:
	void AssembleMatInjIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
	void AssembleMatProdIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;

};


class PeacemanWellT : public PeacemanWell
{
public:
	void CalResFIM(OCP_USI& wId, OCPRes& res, const Bulk& bk, const OCP_DBL& dt) const override;
	void AssembleMatFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const override;
protected:
	void AssembleMatInjFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
	void AssembleMatProdFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;

public:
	void AssembleMatIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const override { OCP_ABORT("NOT USED!"); }
};

#endif /* end if __WELLPEACEMAN_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/17/2023      Create file                          */
/*----------------------------------------------------------------------------*/