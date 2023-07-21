/*! \file    OCPFuncTable.cpp
 *  \brief   Table Functions in OCP
 *  \author  Shizhe Li
 *  \date    Jul/21/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


// OpenCAEPoroX header files
#include "OCPFuncTable.hpp"

/////////////////////////////////////////////////////
// OCPFuncTable2
/////////////////////////////////////////////////////


void OCPFuncTable2::Setup(const TableSet& ts) 
{
	numtable = ts.data.size();
	ref      = ts.refData;
	ref.resize(numtable, 0);

	for (USI i = 0; i < numtable; i++) {
		tables.push_back(OCPTable(ts.data[i]));
	}

	lendata = ts.colNum - 1;
	data.resize(lendata);
	cdata1.resize(lendata);
	cdata2.resize(lendata);
}


void OCPFuncTable2::Eval(const OCP_DBL& val1, const OCP_DBL& val2, vector<OCP_DBL>& out) const
{
	if (numtable == 1)             tables[0].Eval_All0(val2, out);
	else if (val1 <= ref.front())  tables.front().Eval_All0(val2, out);
	else if (val1 >= ref.front())  tables.back().Eval_All0(val2, out);
	else {
		for (USI i = 0; i < numtable - 1; i++) {
			if (val1 <= ref[i + 1]) {
				OCP_DBL w = (val1 - ref[i]) / (ref[i + 1] - ref[i]);
				tables[i].Eval_All0(val2, out);
				tables[i + 1].Eval_All0(val2, data);

				for (USI j = 0; j < lendata; j++) {
					out[j] = w * out[j] + (1 - w) * data[j];
				}
				break;
			}
		}
	}
}


void OCPFuncTable2::Eval(const OCP_DBL& val1, const OCP_DBL& val2, vector<OCP_DBL>& out,
	vector<OCP_DBL>& slope1, vector<OCP_DBL>& slope2) const
{
	fill(slope1.begin(), slope1.end(), 0.0);
	if (numtable == 1)             tables[0].Eval_All0(val2, out, slope2);
	else if (val1 <= ref.front())  tables.front().Eval_All0(val2, out, slope2);
	else if (val1 >= ref.front())  tables.back().Eval_All0(val2, out, slope2);
	else {
		for (USI i = 0; i < numtable - 1; i++) {
			if (val1 <= ref[i + 1]) {
				OCP_DBL w = (val1 - ref[i]) / (ref[i + 1] - ref[i]);
				tables[i].Eval_All0(val2, out, slope2);
				tables[i + 1].Eval_All0(val2, data, slope2);

				for (USI j = 0; j < lendata; j++) {
					slope1[j] = (data[j] - out[j]) / (ref[i + 1] - ref[i]);
					out[j]    = w * out[j] + (1 - w) * data[j];
					slope2[j] = w * slope2[j] + (1 - w) * slope2[j];
				}
				break;
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/21/2023      Create file                          */
/*----------------------------------------------------------------------------*/

