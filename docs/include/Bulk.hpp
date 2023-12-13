/*! \file    Bulk.hpp
 *  \brief   Bulk类和BulkTypeAIM类的声明
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *  本文件包含Bulk类和BulkTypeAIM类的声明，用于描述油藏模拟中的物理信息和动态类型指示。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BULK_HEADER__
#define __BULK_HEADER__

// 标准头文件
#include <iostream>
#include <vector>

// OpenCAEPoroX头文件
#include "DenseMat.hpp"
#include "FlowUnit.hpp"
#include "LinearSystem.hpp"
#include "ParamReservoir.hpp"
#include "Domain.hpp"
#include "PreParamGridWell.hpp"
#include "BulkVarSet.hpp"
#include "PVTModule.hpp"
#include "SATModule.hpp"
#include "ROCKModule.hpp"
#include "BulkInitializer.hpp"
#include "BulkAccumuModule.hpp"
#include "OptionalModules.hpp"

using namespace std;

/**
 * @brief 动态类型指示器类
 * 
 * 用于存储和操作油藏动态类型的指示器。
 */
class BulkTypeAIM
{
public:
    /**
     * @brief 设置动态类型指示器的大小
     * @param nb 指示器的数量
     */
    void Setup(const OCP_USI& nb);

    /**
     * @brief 初始化动态类型指示器的值
     * @param flag 初始化标志
     */
    void Init(const OCP_INT& flag);

    /**
     * @brief 设置特定动态类型指示器的值
     * @param n 指示器索引
     * @param flag 指示器的值
     */
    void SetBulkType(const OCP_USI& n, const OCP_INT& flag);

    /**
     * @brief 判断是否为FIM类型的bulk
     * @param n 指示器索引
     * @return 是否为FIM类型
     */
    OCP_BOOL IfFIMbulk(const OCP_USI& n) const;

    /**
     * @brief 判断是否为IMPEC类型的bulk
     * @param n 指示器索引
     * @return 是否为IMPEC类型
     */
    OCP_BOOL IfIMPECbulk(const OCP_USI& n) const;

    /**
     * @brief 获取特定动态类型指示器的值
     * @param n 指示器索引
     * @return 指示器的值
     */
    OCP_INT GetBulkType(const OCP_USI& n) const;

    /**
     * @brief 获取FIM类型的bulk数量
     * @return FIM类型bulk的数量
     */
    OCP_USI GetNumFIMBulk() const;

protected:
    mutable OCP_USI numFIMBulk; ///< FIM类型bulk的数量
    vector<OCP_INT> indicator;  ///< 存储FIM bulk的索引，FIM bulk: >=0; IMPEC bulk: <0;
};

/**
 * @brief 油藏bulk的物理信息类
 * 
 * Bulk类包含活动网格的主要物理信息。它描述了用于模拟的实际几何域。
 * 变量按bulk存储，然后按相位，然后按组分。Bulk按字母顺序排列，即先按X轴索引，然后是Y轴和Z轴索引。
 * 这里还定义了对每个bulk的操作。
 */
class Bulk
{
    // 为了简化，以下仅展示部分友元类和函数声明
    friend class BulkConn;
    // ... 其他友元类和函数

    /////////////////////////////////////////////////////////////////////
    // 输入参数和设置
    /////////////////////////////////////////////////////////////////////
public:
    /**
     * @brief 从内部数据结构ParamReservoir输入参数
     * @param rs_param 油藏参数
     */
    void InputParam(const ParamReservoir& rs_param);

    /**
     * @brief 为等温模型的流体网格分配内存
     */
    void Setup();

    // ... 其他成员函数

protected:
    BulkVarSet        vs;       ///< 基础变量集合
    PVTModule         PVTm;     ///< PVT模块
    SATModule         SATm;     ///< 饱和度模块
    ROCKModule        ROCKm;    ///< 岩石模块
    BulkInitializer   INITm;    ///< 初始化器
    BulkAccumuModule  ACCm;     ///< 积累项
    OptionalModules   optMs;    ///< 可选模块

    // ... 其他成员变量和函数

public:
    /**
     * @brief 返回bulk的数量
     * @return bulk的数量
     */
    OCP_USI GetBulkNum() const;

    // ... 其他成员函数

protected:
    OCP_DBL          rsTemp;    ///< 油藏温度

    // ... 其他成员变量和函数

public:
    /**
     * @brief 添加一个元素到wellBulkId
     * @param n bulk的索引
     */
    void AddWellBulkId(const OCP_USI& n);

protected:
    vector<OCP_USI> wellBulkId; ///< 被井穿透的bulks的索引及其K-邻居
    BulkTypeAIM bulkTypeAIM;    ///< Bulk类型的AIM

    // ... 其他成员变量和函数
};

#endif /* end if __BULK_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/09/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/