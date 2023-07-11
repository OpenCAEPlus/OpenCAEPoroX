/*! \file    OCPFuncTable.hpp
 *  \brief   Functions for Saturations in OCP
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
	OCPTable          table;
	vector<OCP_DBL>   data;
	vector<OCP_DBL>   cdata;
};



#endif // __OCPFUNCTABLE_HEADER__

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/11/2023      Create file                          */
/*----------------------------------------------------------------------------*/

