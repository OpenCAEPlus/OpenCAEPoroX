/*! 
 * \file    PVTModule.hpp
 * \brief   PVT模块类声明
 * \author  Shizhe Li
 * \date    Aug/19/2023
 *
 * 本文件包含PVTModule类的声明，该类用于处理油藏的PVT（压力-体积-温度）数据。
 * PVTModule类负责管理和设置PVT区域，以及提供对PVT数据的访问。
 *
 *-----------------------------------------------------------------------------------
 * Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 * Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PVTMODULE_HEADER__
#define __PVTMODULE_HEADER__

// 标准头文件
#include <cassert>

// OpenCAEPoroX头文件
#include "MixtureUnit.hpp"
#include "BulkVarSet.hpp"

using namespace std;

/*!
 * \class PVTModule
 * \brief PVT模块类
 *
 * PVTModule类用于管理油藏的PVT区域，并提供对PVT数据的访问。
 */
class PVTModule
{
public:
    /*!
     * \brief 设置PVT模块
     *
     * 根据油藏参数、体积变量集合和可选模块初始化PVT模块。
     * \param rs_param 油藏参数
     * \param bvs 体积变量集合
     * \param opts 可选模块
     */
    void Setup(const ParamReservoir& rs_param, BulkVarSet& bvs, OptionalModules& opts);

    /*!
     * \brief 获取指定PVT区域
     *
     * 返回指定索引的PVT区域的指针。
     * \param n 索引
     * \return 指向PVT区域的指针
     */
    auto GetPVT(const OCP_USI& n) const;

    /*!
     * \brief 获取混合物
     *
     * 返回第一个PVT区域的混合物对象。
     * \return 混合物对象
     */
    auto GetMixture() const;

    /*!
     * \brief 获取PVTNUM
     *
     * 返回PVTNUM向量的引用。
     * \return PVTNUM向量的引用
     */
    auto& GetPVTNUM();

    /*!
     * \brief 获取混合类型
     *
     * 返回混合类型。
     * \return 混合类型
     */
    auto GetMixtureType() const;

    /*!
     * \brief 输出迭代信息
     *
     * 输出第一个PVT区域的混合物迭代信息。
     * \param i 索引
     */
    void OutputIters(const USI& i) const;

protected:
    /// 混合类型
    OCPMixtureType       mixType;
    /// PVT区域的数量
    USI                  NTPVT;
    /// 每个体积的PVT区域索引
    vector<USI>          PVTNUM;
    /// PVT模块集合
    vector<MixtureUnit>  PVTs;
};

#endif /* end if __PVTMODULE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/19/2023      Create file                          */
/*----------------------------------------------------------------------------*/