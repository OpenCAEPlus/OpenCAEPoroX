/*! \file    OCPFlow.cpp
 *  \brief   OCPFlow class declaration
 *  \author  Shizhe Li
 *  \date    Oct/04/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPFlow.hpp"


/////////////////////////////////////////////////////
// OCPFlow
/////////////////////////////////////////////////////


OCPFlow::OCPFlow(const ParamReservoir& rs_param, const USI& i)
{
	if (rs_param.SGOF_T.data.size() > 0 && rs_param.SWOF_T.data.size() > 0) {
		pfMethod = new OCPFlowMethod_OGW01(rs_param.SGOF_T.data[i], rs_param.SWOF_T.data[i], 1, vs);
	}
	else if (rs_param.SOF3_T.data.size() > 0 &&
		rs_param.SGFN_T.data.size() > 0 &&
		rs_param.SWFN_T.data.size() > 0) {
		pfMethod = new OCPFlowMethod_OGW02(rs_param.SOF3_T.data[i], rs_param.SGFN_T.data[i], 
			rs_param.SWFN_T.data[i], 1, vs);
	}
	else if (rs_param.SWOF_T.data.size() > 0) {
		pfMethod = new OCPFlowMethod_OW01(rs_param.SWOF_T.data[i], vs);
	}
	else if (rs_param.SGOF_T.data.size() > 0) {
		pfMethod = new OCPFlowMethod_OG01(rs_param.SGOF_T.data[i], vs);
	}
	else if (rs_param.BCparam.size() > 0) {
		pfMethod = new OCPFlowMethod_GW01(rs_param.BCparam[i], vs);
	}
	else if (rs_param.SWGF_T.data.size() > 0) {
		pfMethod = new OCPFlowMethod_GW02(rs_param.SWGF_T.data[i], vs);
	}
	else {
		OCP_ABORT("NO MATCHED METHOD!");
	}
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/04/2023     Create file                          */
/*----------------------------------------------------------------------------*/