/*! 
 *  \file    OCPFlow.hpp
 *  \brief   OCPFlow 类的声明
 *  \author  Shizhe Li
 *  \date    Jul/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPFLOW_HEADER__
#define __OCPFLOW_HEADER__

#include "ParamReservoir.hpp"
#include "OCPFlowMethod.hpp"
using namespace std;

/**
 * \class OCPFlow
 * \brief 模拟油气水三相流动的类，以油为参考相
 *
 * OCPFlow 类用于处理油气水三相流的计算，包括饱和度、相渗透率和毛细管压力的计算。
 */
class OCPFlow
{
public:
    /**
     * \brief 构造函数
     * \param rs_param 油藏参数
     * \param i 流动类型索引
     */
    OCPFlow(const ParamReservoir& rs_param, const USI& i);

    /**
     * \brief 获取变量集合
     * \return 变量集合的引用
     */
    auto& GetVarSet() { return vs; }

    /**
     * \brief 计算相渗透率和毛细管压力
     * \param S 饱和度数组
     */
    void CalKrPc(const OCP_DBL* S);

    /**
     * \brief 获取流动类型
     * \return 流动类型
     */
    auto FlowType() const { return vs.flowType; }

    /**
     * \brief 计算相渗透率和毛细管压力的导数
     * \param S 饱和度数组
     */
    void CalKrPcDer(const OCP_DBL* S);

    /**
     * \brief 获取油水临界饱和度
     * \return 油水临界饱和度
     */
    OCP_DBL GetSwco() const { return pfMethod->GetSwco(); }

    /**
     * \brief 获取最大毛细管压力（油水）
     * \return 最大毛细管压力（油水）
     */
    OCP_DBL GetMaxPcow() const { return pfMethod->GetMaxPcow(); }

    /**
     * \brief 获取最小毛细管压力（油水）
     * \return 最小毛细管压力（油水）
     */
    OCP_DBL GetMinPcow() const { return pfMethod->GetMinPcow(); }

    /**
     * \brief 通过水饱和度计算毛细管压力（油水）
     * \param Sw 水饱和度
     * \return 毛细管压力（油水）
     */
    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const;

    /**
     * \brief 通过毛细管压力（油水）计算水饱和度
     * \param Pcow 毛细管压力（油水）
     * \return 水饱和度
     */
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const;

    /**
     * \brief 通过气饱和度计算毛细管压力（油气）
     * \param Sg 气饱和度
     * \return 毛细管压力（油气）
     */
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const;

    /**
     * \brief 通过毛细管压力（油气）计算气饱和度
     * \param Pcgo 毛细管压力（油气）
     * \return 气饱和度
     */
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const;

    /**
     * \brief 通过毛细管压力（气水）计算水饱和度
     * \param Pcgw 毛细管压力（气水）
     * \return 水饱和度
     */
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const;

    /**
     * \brief 计算气相渗透率
     * \param Sg 气饱和度
     * \param dKrgdSg 气相渗透率对气饱和度的导数
     * \return 气相渗透率
     */
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const;

protected:
    /**
     * \brief 设置饱和度
     * \param S 饱和度数组
     */
    void SetSaturation(const OCP_DBL* S);

protected:
    OCPFlowVarSet vs; ///< 变量集合
    OCPFlowMethod* pfMethod; ///< 流动方法指针
};

#endif /* end if __OCPFLOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/