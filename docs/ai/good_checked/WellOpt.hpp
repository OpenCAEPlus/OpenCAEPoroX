```cpp
/*! 
 * \file    WellOpt.hpp
 * \brief   本文件包含了WellOpt类的声明
 * \author  Shizhe Li
 * \date    Nov/22/2022
 *
 * 本文件定义了WellOpt类，用于描述油井的操作模式，并且通常随时间变化。
 * 包括油井类型、注入流体、油井状态、控制模式等属性。
 *
 * -----------------------------------------------------------------------------------
 * Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 * Released under the terms of the GNU Lesser General Public License 3.0 or later.
 * -----------------------------------------------------------------------------------
 */

#ifndef __WELLOPT_HEADER__
#define __WELLOPT_HEADER__

// 标准头文件
#include <cmath>

// OpenCAEPoroX头文件
#include "ParamWell.hpp"
#include "OCPMixture.hpp"

using namespace std;

/// 油井类型枚举
enum class WellType : USI
{
    injector,    ///< 注入井
    productor    ///< 生产井
};

/// 油井状态枚举
enum class WellState : USI
{
    open,   ///< 开启
    close   ///< 关闭
};

/// 油井操作模式枚举
enum class WellOptMode : USI
{
    irate,  ///< 注入井流量控制
    orate,  ///< 生产井油流量控制
    grate,  ///< 生产井气流量控制
    wrate,  ///< 生产井水流量控制
    lrate,  ///< 生产井液体流量控制
    bhp,    ///< 注入/生产井底压力控制
};

/// \class WellOpt
/// \brief 描述油井的操作模式
class WellOpt
{
    friend class Well;
    friend class AllWells;
    friend class Out4RPT;

public:
    /// 默认构造函数
    WellOpt() = default;

    /// 使用参数构造油井操作模式
    WellOpt(const WellOptParam& Optparam);

    /// 重载不等运算符
    OCP_BOOL operator!=(const WellOpt& Opt) const;

public:
    WellType type;                     ///< 油井类型，注入井或生产井
    string   injFluidName;             ///< 注入流体名称，仅对注入井有用
    WellState state{ WellState::close }; ///< 油井状态，默认为关闭
    WellOptMode mode;                  ///< 油井控制模式
    WellOptMode initMode;              ///< 初始控制模式
    OCP_DBL maxRate;                   ///< 流量上限或指定流量
    OCP_DBL maxBHP;                    ///< 注入井最高压力或井口压力
    OCP_DBL minBHP;                    ///< 生产井最低压力或井口压力
    OCP_DBL tarRate;                   ///< 目标注入/生产井流量
    OCP_DBL tarBHP;                    ///< 目标注入/生产井底压力
    vector<OCP_DBL> injZi;             ///< 注入流体组分
    PhaseType       injPhase;          ///< 注入流体相态
    vector<OCP_DBL> prodPhaseWeight;   ///< 生产相权重
    OCP_DBL injTemp;                   ///< 注入流体温度（单位：F）
};

/// \class SolventINJ
/// \brief 描述从注入井向储层注入的流体的摩尔分数
class SolventINJ
{
public:
    /// 默认构造函数
    SolventINJ() = default;

    /// 使用其他溶剂对象构造
    SolventINJ(const Solvent& other)
    {
        name = other.name;
        data = other.comRatio;
    };

    string          name; ///< 溶剂名称
    vector<OCP_DBL> data; ///< 组件摩尔分数
};

#endif // __WELLOPT_HEADER__

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           NOV/22/2022      Create file                          */
/*----------------------------------------------------------------------------*/
