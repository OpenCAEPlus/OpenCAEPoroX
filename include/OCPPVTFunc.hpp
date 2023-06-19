/*! \file    OCPPVTFunc.hpp
 *  \brief   Funcs for PVT in OCP
 *  \author  Shizhe Li
 *  \date    Jun/18/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPPVTFUNC_HEADER__
#define __OCPPVTFUNC_HEADER__

// OpenCAEPoroX header files
#include "OCPTable.hpp"

using namespace std;

class PVTFunc
{

};

class PVTW_Table
{
public:
	PVTW_Table(const vector<vector<OCP_DBL>>& src);

protected:
	OCPTable table;

};





#endif // __OCPPVTFUNC_HEADER__

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jun/18/2023      Create file                          */
/*----------------------------------------------------------------------------*/

