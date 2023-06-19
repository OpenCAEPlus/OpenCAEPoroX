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

#ifndef __OCPTIMERECORD_HEADER__
#define __OCPTIMERECORD_HEADER__

#include "OCPConst.hpp"


extern OCP_DBL OCPTIME_TOTAL;				///< Time for Total Simulation
extern OCP_DBL OCPTIME_PARTITION;			///< Time for Grid Inputting and Partition
extern OCP_DBL OCPTIME_READPARAM;			///< Time for Inputting Reservoir Params
extern OCP_DBL OCPTIME_SETUP_SIM;			///< Time for Setuping Simulator
extern OCP_DBL OCPTIME_INIT_RESERVOIR;		///< Time for Reservoir Initializing
extern OCP_DBL OCPTIME_ASSEMBLE_MAT;		///< Time for Matrix Assembling
extern OCP_DBL OCPTIME_ASSEMBLE_MAT_FOR_LS; ///< Time for Matrix Assembling used in external linear solver
extern OCP_DBL OCPTIME_LSOLVER;				///< Time for Linear Solver
extern OCP_DBL OCPTIME_NRSTEP;				///< Time for NR step
extern OCP_DBL OCPTIME_NRSTEPC;             ///< Time for main calculation in NRstep
extern OCP_DBL OCPTIME_UPDATE_GRID;			///< Time for Grid Updating
extern OCP_DBL OCPTIME_OUTPUT;				///< Time for Outputing Files
extern OCP_DBL OCPTIME_COMM_COLLECTIVE;     ///< Time for Collective Communication (Disposable calls will be ignored)
extern OCP_DBL OCPTIME_COMM_1ALLREDUCE;     ///< Time for Allreduce in OCP check
extern OCP_DBL OCPTIME_COMM_P2P;            ///< Time for P2P Communication (Disposable calls will be ignored)
extern OCP_DBL OCPTIME_COMM_TOTAL;          ///< Time for Communication (Disposable calls will be ignored)
extern OCP_DBL OCPTIME_PARMETIS;            ///< Time for ParMetis

#endif


 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Apr/07/2023      Create file                          */
 /*----------------------------------------------------------------------------*/
