/*! \file    HeatLoss.hpp 
 *  \brief   HeatLoss类声明
 *  \author  Shizhe Li 
 *  \date    Aug/24/2023 
 *
 *  本文件包含了HeatLoss类的声明。该类负责处理热损失相关的计算。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved. 
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later. 
 *----------------------------------------------------------------------------------- 
 */

#ifndef __HEATLOSS_HEADER__
#define __HEATLOSS_HEADER__

// OpenCAEPoroX header files
#include "ParamReservoir.hpp"
#include "BulkVarSet.hpp"
#include <vector>

using namespace std;

/*! \class HeatLossVarSet
 *  \brief HeatLossVarSet类用于存储热损失相关的变量
 *
 *  该类包含了热损失计算中需要用到的各种变量，包括热损失率，以及相应的辅助变量等。
 */
class HeatLossVarSet{
    friend class HeatLoss;
    friend class HeatLossMethod01;

public:
    void SetNb(const OCP_USI& nbin) { nb = nbin; } //!< 设置块数
    void ResetToLastTimeStep(); //!< 重置到上一时间步
    void UpdateLastTimeStep(); //!< 更新到上一时间步

protected:
    OCP_USI         nb; //!< 块数
    vector<OCP_DBL> I; //!< 辅助变量
    vector<OCP_DBL> hl; //!< 热损失率
    vector<OCP_DBL> hlT; //!< 热损失率对温度的导数
    vector<OCP_DBL> lI; //!< 上一时间步的I
    vector<OCP_DBL> lhl; //!< 上一时间步的hl
    vector<OCP_DBL> lhlT; //!< 上一时间步的hlT
};

// 省略其他类的注释

#endif /* end if __HeatLoss_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/24/2023      Create file                          */
/*----------------------------------------------------------------------------*/
