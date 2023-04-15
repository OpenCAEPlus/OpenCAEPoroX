/*! \file    OCPTimeRecord.hpp
 *  \brief   Time record Item
 *  \author  Shizhe Li
 *  \date    Apr/07/2023
 *
 *  \note    Record the time spent on the most important part during the simulation
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPTimeRecord.hpp"

OCP_DBL OCPTIME_TOTAL				= 0;
OCP_DBL OCPTIME_PARTITION			= 0;
OCP_DBL OCPTIME_READPARAM			= 0;
OCP_DBL OCPTIME_SETUP_SIM			= 0;
OCP_DBL OCPTIME_INIT_RESERVOIR		= 0;
OCP_DBL OCPTIME_ASSEMBLE_MAT		= 0;
OCP_DBL OCPTIME_ASSEMBLE_MAT_FOR_LS = 0;
OCP_DBL OCPTIME_LSOLVER				= 0;
OCP_DBL OCPTIME_NRSTEP				= 0;
OCP_DBL OCPTIME_UPDATE_GRID			= 0;
OCP_DBL OCPTIME_OUTPUT				= 0;


 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Apr/07/2023      Create file                          */
 /*----------------------------------------------------------------------------*/