/*! 
  \file    MixtureUnit.hpp
  \brief   MixtureUnit 类的声明
  \author  Shizhe Li
  \date    Oct/01/2021

  -----------------------------------------------------------------------------------
  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
  Released under the terms of the GNU Lesser General Public License 3.0 or later.
  -----------------------------------------------------------------------------------
*/

#ifndef __MIXTURE_HEADER__
#define __MIXTURE_HEADER__

// 标准头文件
#include <iostream>
#include <vector>

// OpenCAEPoroX 头文件
#include "OCPConst.hpp"
#include "OptionalModules.hpp"
#include "ParamReservoir.hpp"
#include "OCPMixture.hpp"
#include "BulkVarSet.hpp"

using namespace std;

/*!
  \class MixtureUnit
  \brief MixtureUnit 是一个抽象类，包含用于闪蒸计算的所有信息，包括变量和函数。
         它可以计算质量密度等相的属性，具有与bulks中的数据结构相同。
*/
class MixtureUnit
{
public:
    /*!
      \brief 构造函数
      \param rs_param 油藏参数
      \param i 索引
      \param opts 可选模块
    */
    MixtureUnit(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts);

    /*!
      \brief 获取混合类型
      \return 混合类型
    */
    auto GetMixtureType() const { return vs->mixtureType; }

    /*!
      \brief 获取混合对象
      \return 混合对象指针
    */
    OCPMixture* GetMixture() const { return mix; }

    /*!
      \brief 闪蒸计算
      \param Pin 输入压力
      \param Tin 输入温度
      \param Niin 输入组分摩尔数
    */
    void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) const;

    /*!
      \brief 以相的饱和状态进行闪蒸计算
      \param bId bulk的ID
      \param bvs bulk变量集
    */
    void InitFlashIMPEC(const OCP_USI& bId, const BulkVarSet& bvs) const;

    /*!
      \brief 以相的饱和状态进行闪蒸计算
      \param bId bulk的ID
      \param bvs bulk变量集
    */
    void InitFlashFIM(const OCP_USI& bId, const BulkVarSet& bvs) const;

    /*!
      \brief 以组分摩尔数进行闪蒸计算
      \param bId bulk的ID
      \param bvs bulk变量集
    */
    void FlashIMPEC(const OCP_USI& bId, const BulkVarSet& bvs) const;

    /*!
      \brief 以组分摩尔数进行闪蒸计算，并计算导数
      \param bId bulk的ID
      \param bvs bulk变量集
    */
    void FlashFIM(const OCP_USI& bId, const BulkVarSet& bvs) const;

    /*!
      \brief 计算相的质量密度
      \param Pin 输入压力
      \param Tin 输入温度
      \param Ziin 输入组分摩尔分数
      \param pt 相类型
      \return 相的质量密度
    */
    OCP_DBL XiPhase(const OCP_DBL& Pin, const OCP_DBL& Tin, const vector<OCP_DBL>& Ziin, const PhaseType& pt) const;

    /*!
      \brief 计算相的质量密度
      \param Pin 输入压力
      \param Pbb 气泡点压力
      \param Tin 输入温度
      \param Ziin 输入组分摩尔分数
      \param pt 相类型
      \return 相的质量密度
    */
    OCP_DBL RhoPhase(const OCP_DBL& Pin, const OCP_DBL& Pbb, const OCP_DBL& Tin,
                     const vector<OCP_DBL>& Ziin, const PhaseType& pt) const;

    /*!
      \brief 计算注入井的焓
      \param Tin 输入温度
      \param Ziin 输入组分摩尔分数
      \return 注入井的焓
    */
    OCP_DBL CalInjWellEnthalpy(const OCP_DBL& Tin, const OCP_DBL* Ziin) const;

    /*!
      \brief 输出混合迭代次数
    */
    void OutMixtureIters() const;

public:
    // 以下方法用于获取混合变量集合中的各种属性

    const auto GetVs() const;
    const auto& GetNt() const;
    const auto& GetNi(const USI& i) const;
    // ... 更多成员变量获取方法省略 ...

protected:
    /*!
      \brief 组件的混合
    */
    OCPMixture*             mix;

    /*!
      \brief 混合的变量集
    */
    const OCPMixtureVarSet* vs;

protected:
    /*!
      \brief 表面张力
    */
    SurfaceTension* surTen;
    USI             stMethodIndex;

    /*!
      \brief 可混合因子
    */
    MiscibleFactor* misFac;
    USI             mfMethodIndex;
};

#endif /* end if __MIXTURE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/
