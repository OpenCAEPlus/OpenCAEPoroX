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
	/// Initialize Well Pressure
	void InitWellP(const Bulk& bk) { bhp = bk.vs.P[perf[0].location]; }
protected:
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

protected:

	/// difference of pressure between well and perforation: numPerf.
	vector<OCP_DBL> dG;
	/// Last dG
	vector<OCP_DBL> ldG;
};


class PeacemanWellIsoT : public PeacemanWell
{
public:
	void CalResFIM(OCP_USI& wId, OCPRes& res, const Bulk& bk, const OCP_DBL& dt) const;
	void AssembleMatFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
protected:
	void AssembleMatInjFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
	void AssembleMatProdFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;

public:
	void AssembleMatIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
protected:
	void AssembleMatInjIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
	void AssembleMatProdIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;

};


class PeacemanWellT : public PeacemanWell
{
public:
	void CalResFIM(OCP_USI& wId, OCPRes& res, const Bulk& bk, const OCP_DBL& dt) const;
	void AssembleMatFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;	
protected:
	void AssembleMatInjFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
	void AssembleMatProdFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;

public:
	void AssembleMatIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const { OCP_ABORT("NOT USED!"); }
protected:
	void AssembleMatInjIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const { OCP_ABORT("NOT USED!"); }
	void AssembleMatProdIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const { OCP_ABORT("NOT USED!"); }
};

#endif /* end if __WELLPEACEMAN_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/17/2023      Create file                          */
/*----------------------------------------------------------------------------*/