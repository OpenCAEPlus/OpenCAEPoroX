/*! \file    OCPTimeRecord.hpp 
 *  \brief   时间记录项
 *  \author  Shizhe Li 
 *  \date    Apr/07/2023 
 *
 *  \note    记录模拟过程中最重要部分的耗时
 * 
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPTIMERECORD_HEADER__
#define __OCPTIMERECORD_HEADER__

#include "OCPConst.hpp"

/// 1 [s] = 1000 [ms]
const OCP_DBL TIME_S2MS = 1000;

extern OCP_DBL OCPTIME_TOTAL;				///< 总模拟时间
extern OCP_DBL OCPTIME_PARTITION;			///< 网格输入和划分时间
extern OCP_DBL OCPTIME_READPARAM;			///< 输入油藏参数时间
extern OCP_DBL OCPTIME_SETUP_SIM;			///< 模拟器设置时间
extern OCP_DBL OCPTIME_INIT_RESERVOIR;		///< 油藏初始化时间
extern OCP_DBL OCPTIME_ASSEMBLE_MAT;		///< 矩阵组装时间
extern OCP_DBL OCPTIME_CONVERT_MAT_FOR_LS_IF; ///< 外部线性求解器中的矩阵组装时间
extern OCP_DBL OCPTIME_LSOLVER;				///< 线性求解器时间
extern OCP_DBL OCPTIME_NRSTEP;				///< NR步骤时间
extern OCP_DBL OCPTIME_NRSTEPC;             ///< NR步骤中的主要计算时间
extern OCP_DBL OCPTIME_UPDATE_GRID;			///< 网格更新时间
extern OCP_DBL OCPTIME_OUTPUT;				///< 输出文件时间
extern OCP_DBL OCPTIME_COMM_COLLECTIVE;     ///< 集体通信时间（一次性调用将被忽略）
extern OCP_DBL OCPTIME_COMM_1ALLREDUCE;     ///< OCP检查中的Allreduce时间
extern OCP_DBL OCPTIME_COMM_P2P;            ///< 点对点通信时间（一次性调用将被忽略）
extern OCP_DBL OCPTIME_COMM_TOTAL;          ///< 通信时间（一次性调用将被忽略）
extern OCP_DBL OCPTIME_PARMETIS;            ///< ParMetis时间

#endif

/*----------------------------------------------------------------------------*/ 
/*  This File's Change History Brief                                          */ 
/*----------------------------------------------------------------------------*/ 
/*  Author              Date             Actions                              */ 
/*----------------------------------------------------------------------------*/ 
/*  Shizhe Li           Apr/07/2023      创建文件                            */ 
/*----------------------------------------------------------------------------*/

