/*! 
 *  \file    BulkAccumuModule.hpp
 *  \brief   BulkAccumuModule 类的声明
 *  \author  Shizhe Li
 *  \date    Aug/30/2023
 *
 *  本文件包含 BulkAccumuModule 类及其派生类的声明，用于处理油藏模拟中的累积项计算。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BULKACCUMUTERM_HEADER__
#define __BULKACCUMUTERM_HEADER__

// OpenCAEPoroX header files
#include "BulkVarSet.hpp"
#include "OptionalModules.hpp"
#include <vector>

using namespace std;

/*!
 *  \class   BulkAccumuTerm
 *  \brief   累积项计算的抽象基类
 *
 *  该类定义了累积项计算的接口，并提供了一些基本的成员变量和函数。
 */
class BulkAccumuTerm
{
public:
    /// 计算 FIM 的残差
    virtual const vector<OCP_DBL>& CalResFIM(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt) const = 0;
    /// 计算 FIM 的 dFdXp
    virtual const vector<OCP_DBL>& CaldFdXpFIM(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt) const = 0;
    /// 计算 IMPEC 的对角线值和右手边
    virtual void CalValRhsIMPEC(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt, OCP_DBL& valout, OCP_DBL& rhsout) const = 0;

protected:
    // for FIM
    /// dFdXp 的维度
    USI                     dim;    ///< 类型: USI, 用途: 存储 dFdXp 的维度, 默认值: 未指定
    /// dF(累积项) / dXp (完全导数)
    mutable vector<OCP_DBL> dFdXp;  ///< 类型: vector<OCP_DBL>, 用途: 存储 dF(累积项) 对 dXp(未知数) 的导数, 默认值: 空向量
    /// 残差 (累积项)
    mutable vector<OCP_DBL> res;    ///< 类型: vector<OCP_DBL>, 用途: 存储累积项的残差, 默认值: 空向量

    // dependent module
    /// 依赖的模块
    const OptionalModules*  optM;   ///< 类型: const OptionalModules*, 用途: 指向依赖的模块, 默认值: nullptr
};

/*!
 *  \class   BulkAccumuTerm01
 *  \brief   等温模型的累积项计算类
 *
 *  用于等温条件下的累积项计算。
 */
class BulkAccumuTerm01 : public BulkAccumuTerm
{
public:
    /// 构造函数
    BulkAccumuTerm01(const BulkVarSet& bvs, const OptionalModules* opt);
    /// 计算 FIM 的残差
    const vector<OCP_DBL>& CalResFIM(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt) const override;
    /// 计算 FIM 的 dFdXp
    const vector<OCP_DBL>& CaldFdXpFIM(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt) const override;
    /// 计算 IMPEC 的对角线值和右手边
    void CalValRhsIMPEC(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt, OCP_DBL& valout, OCP_DBL& rhsout) const override;
};

/*!
 *  \class   BulkAccumuTerm02
 *  \brief   热力学模型的累积项计算类
 *
 *  用于热力学条件下的累积项计算。
 */
class BulkAccumuTerm02 : public BulkAccumuTerm
{
public:
    /// 构造函数
    BulkAccumuTerm02(const BulkVarSet& bvs, const OptionalModules* opt);
    /// 计算 FIM 的残差
    const vector<OCP_DBL>& CalResFIM(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt) const override;
    /// 计算 FIM 的 dFdXp
    const vector<OCP_DBL>& CaldFdXpFIM(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt) const override;
    /// 计算 IMPEC 的对角线值和右手边
    void CalValRhsIMPEC(const OCP_USI& bId, const BulkVarSet& bvs, const OCP_DBL& dt, OCP_DBL& valout, OCP_DBL& rhsout) const override { OCP_ABORT("NOT USED!"); }
};

/*!
 *  \class   BulkAccumuModule
 *  \brief   累积项模块类
 *
 *  用于设置和获取累积项计算相关的类和数据。
 */
class BulkAccumuModule
{
public:
    /// 设置累积模块
    void Setup(const ParamReservoir& param, const BulkVarSet& bvs, const OptionalModules& opt);
    /// 获取累积项
    const auto GetAccumuTerm() const { return bacT; }

protected:
    /// 累积项
    BulkAccumuTerm* bacT; ///< 类型: BulkAccumuTerm*, 用途: 指向累积项计算的对象, 默认值: nullptr
};

#endif /* end if __BULKACCUMUMODULE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/30/2023      Create file                          */
/*----------------------------------------------------------------------------*/
