/*! \file    BulkAccumuModule.hpp
 *  \brief   BulkAccumuModule class declaration
 *  \author  Shizhe Li
 *  \date    Aug/30/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BULKACCUMUTERM_HEADER__
#define __BULKACCUMUTERM_HEADER__


// OpenCAEPoroX header files
#include "BulkVarSet.hpp"
#include "OptionalModules.hpp"

#include <vector>

using namespace std;


class BulkAccumuTerm
{
public:
	/// Calculate residual for FIM
	virtual const vector<OCP_DBL>& CalResFIM(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt) const = 0;
	/// Calculate dFdXp for FIM
	virtual const vector<OCP_DBL>& CaldFdXpFIM(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt) const = 0;
	/// Calculate diagnal value and rhs for IMPEC
	virtual void CalValRhsIMPEC(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt, OCP_DBL& valout, OCP_DBL& rhsout) const = 0;

protected:

	// for FIM
	/// dim of dFdXp
	USI                     dim;
	/// dF(accumulation term) / dXp (full derivatives)
	mutable vector<OCP_DBL> dFdXp;
	/// residual (accumulation term)
	mutable vector<OCP_DBL> res;
					      
	// dependent module
	const OptionalModules*  optM;
};


/// for Isothermal model
class BulkAccumuTerm01 : public BulkAccumuTerm
{
public:
	/// constructor
	BulkAccumuTerm01(const BulkVarSet& bvs, const OptionalModules* opt);
	/// Calculate residual for FIM
	const vector<OCP_DBL>& CalResFIM(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt) const override;
	/// Calculate dFdXp for FIM
	const vector<OCP_DBL>& CaldFdXpFIM(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt) const override;
	/// Calculate diagnal value and rhs for IMPEC
	void CalValRhsIMPEC(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt, OCP_DBL& valout, OCP_DBL& rhsout) const override;
};


/// for thermal model
class BulkAccumuTerm02 : public BulkAccumuTerm
{
public:
	/// constructor
	BulkAccumuTerm02(const BulkVarSet& bvs, const OptionalModules* opt);
	/// Calculate residual for FIM
	const vector<OCP_DBL>& CalResFIM(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt) const override;
	/// Calculate dFdXp for FIM
	const vector<OCP_DBL>& CaldFdXpFIM(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt) const override;
	/// Calculate diagnal value and rhs for IMPEC
	void CalValRhsIMPEC(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt, OCP_DBL& valout, OCP_DBL& rhsout) const override { OCP_ABORT("NOT USED!"); }
};


// Bulk accumulation term
class BulkAccumuModule
{
public:
	/// Setup accumulation module
	void Setup(const ParamReservoir& param, const BulkVarSet& bvs, const OptionalModules& opt);
	/// Get accumulation term
	const auto GetAccumuTerm() const { return bacT; }
protected:
	/// accumulation term
	BulkAccumuTerm* bacT;
};



#endif /* end if __BULKACCUMUMODULE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/30/2023      Create file                          */
/*----------------------------------------------------------------------------*/
