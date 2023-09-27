/*! \file    OCPFuncTable.hpp
 *  \brief   Table Functions in OCP
 *  \author  Shizhe Li
 *  \date    Jul/11/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPFUNCTABLE_HEADER__
#define __OCPFUNCTABLE_HEADER__

// OpenCAEPoroX header files
#include "OCPTable.hpp"


using namespace std;


/////////////////////////////////////////////////////
// OCPFuncTable
/////////////////////////////////////////////////////


class OCPFuncTable
{
public:
	OCPFuncTable() = default;
	void Setup(const vector<vector<OCP_DBL>>& src) {
		table.Setup(src);
		data.resize(table.GetColNum());
		cdata.resize(table.GetColNum());
	}
	OCP_BOOL IsEmpty() const { return table.IsEmpty(); }

protected:
	OCPTable                  table;
	mutable vector<OCP_DBL>   data;
	mutable vector<OCP_DBL>   cdata;
};


/////////////////////////////////////////////////////
// OCPFuncTable2
/////////////////////////////////////////////////////


class OCPFuncTable2
{
public:
	OCPFuncTable2() = default;
	void Setup(const Table2& tab) {
		table.Setup(tab);
		data.resize(table.GetColNum());
		cdata1.resize(table.GetColNum());
		cdata2.resize(table.GetColNum());
	}
	OCP_BOOL IsEmpty() const { return table.IsEmpty(); }

protected:
	OCPTable2                 table;
	mutable vector<OCP_DBL>   data;
	mutable vector<OCP_DBL>   cdata1;
	mutable vector<OCP_DBL>   cdata2;
};



#endif // __OCPFUNCTABLE_HEADER__

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/11/2023      Create file                          */
/*----------------------------------------------------------------------------*/

