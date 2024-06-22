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

OCP_DBL OCPTIME_TOTAL				  = 0;
OCP_DBL OCPTIME_PARTITION			  = 0;
OCP_DBL OCPTIME_PARMETIS              = 0;
OCP_DBL OCPTIME_READPARAM			  = 0;
OCP_DBL OCPTIME_SETUP_SIM			  = 0;
OCP_DBL OCPTIME_INIT_RESERVOIR		  = 0;
OCP_DBL OCPTIME_ASSEMBLE_MAT		  = 0;
OCP_DBL OCPTIME_CONVERT_MAT_FOR_LS_IF = 0;
OCP_DBL OCPTIME_LSOLVER				  = 0;
OCP_DBL OCPTIME_NRSTEP				  = 0;
OCP_DBL OCPTIME_NRSTEPC               = 0;
OCP_DBL OCPTIME_UPDATE_GRID			  = 0;
OCP_DBL OCPTIME_OUTPUT				  = 0;
OCP_DBL OCPTIME_COMM_COLLECTIVE       = 0;
OCP_DBL OCPTIME_COMM_1ALLREDUCE       = 0;
OCP_DBL OCPTIME_COMM_P2P              = 0;
OCP_DBL OCPTIME_COMM_TOTAL            = 0;
OCP_DBL OCPTIME_GROUPPROCESS          = 0;
OCP_DBL OCPTIME_LSOLVER_DDM           = 0;
OCP_DBL OCPITER_NR_DDM                = 0;
OCP_DBL OCPITER_LS_DDM                = 0;

 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Apr/07/2023      Create file                          */
 /*----------------------------------------------------------------------------*/