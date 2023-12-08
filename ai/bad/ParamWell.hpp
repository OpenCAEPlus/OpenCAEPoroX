```cpp
/*! 
 * \file    ParamWell.hpp
 * \brief   ParamWell类声明文件
 * \author  Shizhe Li
 * \date    Oct/01/2021
 *
 * 本文件包含ParamWell类及其相关类的声明，这些类用于存储和处理油井的参数信息。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARAMWELL_HEADER__
#define __PARAMWELL_HEADER__

// 标准头文件
#include <cassert>
#include <fstream>
#include <vector>

// OpenCAEPoroX头文件
#include "OCPConst.hpp"
#include "UtilInput.hpp"
#include "UtilOutput.hpp"

using namespace std;

/*!
 * \class   WellOptParam
 * \brief   包含油井操作所需的所有参数
 *
 * WellOptParam类存储了一个油井在不同时间可能改变的所有操作参数。
 */
class WellOptParam
{
public:
    /*!
     * \brief   构造函数，初始化油井操作参数
     * \param   intype 操作类型，注入或生产
     * \param   vbuf   参数值的字符串向量
     */
    WellOptParam(string intype, vector<string>& vbuf);
    
    string type;      ///< 油井类型，注入还是生产？
    string fluidType; ///< 注入井的注入流体类型（仅限注入井）
    string state;     ///< 油井状态，开启或关闭？
    string mode;      ///< 油井模式，按率控还是按底压控？
    OCP_DBL maxRate;  ///< 允许的最大流入/流出油井的流速
    OCP_DBL maxBHP;   ///< 注入井允许的最大井底压力
    OCP_DBL minBHP;   ///< 生产井允许的最小井底压力
    OCP_DBL injTemp;  ///< 注入流体的温度
};

/*!
 * \class   WellOptPair
 * \brief   包含油井操作模式及其应用的开始时间
 *
 * WellOptPair类包含两部分信息：一部分是油井的操作模式，另一部分是操作开始应用的时间，
 * 该时间通过在关键时间点的索引来表示。
 */
class WellOptPair
{
public:
    /*!
     * \brief   构造函数，初始化油井操作模式和开始时间
     * \param   i     关键时间点的索引
     * \param   type  操作类型
     * \param   vbuf  参数值的字符串向量
     */
    WellOptPair(USI i, string type, vector<string>& vbuf);
    
    USI          d;    ///< 操作开始应用的时间点的索引
    WellOptParam opt;  ///< 油井的操作参数
};

/*!
 * \class   WellParam
 * \brief   包含油井静态和动态参数信息
 *
 * WellParam类用于存储和处理油井的静态信息（如井的位置和尺寸）和动态信息（如操作模式）。
 */
class WellParam
{
public:
    // 函数声明略去，根据实际代码补充

    // 静态信息
    GridType gridType;    ///< 网格类型
    string   name;        ///< 油井名称
    string   group{ "FEILD" }; ///< 油井所属的组
    USI      I;           ///< 结构化网格中油井头的I索引
    USI      J;           ///< 结构化网格中油井头的J索引
    OCP_DBL  depth{ -1.0 }; ///< 油井头的深度
    OCP_DBL  X;           ///< 非结构化网格中油井头的x坐标
    OCP_DBL  Y;           ///< 非结构化网格中油井头的y坐标
    OCP_DBL  Z;           ///< 非结构化网格中油井头的z坐标

    // 动态信息
    vector<WellOptPair> optParam; ///< 油井的操作参数序列
};

/*!
 * \class   Solvent
 * \brief   描述从INJ注入到储层的流体组分的摩尔分数
 */
class Solvent
{
public:
    // 函数声明略去，根据实际代码补充

    string          name;     ///< 溶剂名称
    vector<OCP_DBL> comRatio; ///< 组件比例
};

/*!
 * \class   ParamWell
 * \brief   用于存储从输入文件中获取的油井信息的内部结构
 *
 * ParamWell类是一个中间接口，独立于主模拟器。在所有文件输入完成后，其中的参数将传递给相应的模块。
 */
class ParamWell
{
public:
    // 函数声明略去，根据实际代码补充

    OCP_BOOL          thermal{ OCP_FALSE }; ///< 是否使用热模型
    vector<WellParam> well;                 ///< 所有油井的信息
    vector<OCP_DBL>   criticalTime;         ///< 用户给定的关键时间记录
    vector<Solvent>   solSet;               ///< 溶剂集合
    OCP_DBL           Psurf;                ///< 表面压力
    OCP_DBL           Tsurf;                ///< 表面温度
};

#endif /* end if __PARAMWELL_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/
```

