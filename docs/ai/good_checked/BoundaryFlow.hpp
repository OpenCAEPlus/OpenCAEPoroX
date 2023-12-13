/*! \file    BoundaryFlow.hpp
 *  \brief   边界流动类的声明
 *  \details 该文件包含了边界流动类的声明，用于模拟和管理边界流动的行为。
 *  \author  Shizhe Li
 *  \date    Sep/27/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BOUNDARYFLOW_HEADER__
#define __BOUNDARYFLOW_HEADER__

// OpenCAEPoroX header files
#include "ParamReservoir.hpp"
#include "BulkVarSet.hpp"
#include <vector>
using namespace std;

/*!
 * \class BoundaryFlowVarSet
 * \brief 边界流动变量集合类
 * \details 用于存储和管理边界流动过程中的变量集合。
 */
class BoundaryFlowVarSet
{
    friend class BoundaryFlow;
    friend class BoundaryFlowMethod01;

public:
    /*!
     * \brief 重置到上一个时间步
     */
    void ResetToLastTimeStep()
    {
    }

    /*!
     * \brief 更新到上一个时间步
     */
    void UpdateLastTimeStep()
    {
    }

protected:
};

/*!
 * \class BoundaryFlowMethod
 * \brief 边界流动方法基类
 * \details 提供了边界流动计算方法的基本接口。
 */
class BoundaryFlowMethod
{
public:
    /*!
     * \brief 默认构造函数
     */
    BoundaryFlowMethod() = default;
};

/*!
 * \class BoundaryFlowMethod01
 * \brief 常压流动方法类
 * \details 根据给定的边界参数和变量集合，实现常压流动计算。
 */
class BoundaryFlowMethod01 : public BoundaryFlowMethod
{
public:
    /*!
     * \brief 构造函数
     * \param bP 边界参数
     * \param bfvs 边界流动变量集合引用
     */
    BoundaryFlowMethod01(const BoundaryParam& bP, BoundaryFlowVarSet& bfvs)
    {
        name = bP.name;
        P    = bP.P;
    }

protected:
    string  name; /*!< 名称 */
    OCP_DBL P;    /*!< 压力 */
};

/*!
 * \class BoundaryFlow
 * \brief 边界流动类
 * \details 管理边界流动计算的主要类，包含了流动变量集合、计算方法等。
 */
class BoundaryFlow
{
    friend class BulkAccumuTerm01;

public:
    /*!
     * \brief 默认构造函数
     */
    BoundaryFlow() = default;

    /*!
     * \brief 判断是否使用此边界流动
     * \param n 边界索引
     * \return 是否使用
     */
    auto IfUse(const OCP_USI& n) const
    {
        if (!ifUse)
            return OCP_FALSE;
        else if (mIndex[n] < 0)
            return OCP_FALSE;
        else
            return OCP_TRUE;
    }

    /*!
     * \brief 设置边界流动
     * \param rs_param 油藏参数
     * \param bvs 整体变量集合
     * \param boundName 边界名称列表
     * \param boundIndex 边界索引列表
     */
    void Setup(const ParamReservoir& rs_param, const BulkVarSet& bvs, const vector<string>& boundName, const vector<USI>& boundIndex);

    /*!
     * \brief 重置到上一个时间步
     */
    void ResetToLastTimeStep() { if (ifUse)  vs.ResetToLastTimeStep(); }

    /*!
     * \brief 更新到上一个时间步
     */
    void UpdateLastTimeStep() { if (ifUse)  vs.UpdateLastTimeStep(); }

protected:
    OCP_BOOL                    ifUse{ OCP_FALSE }; /*!< 是否使用热损失 */
    BoundaryFlowVarSet          vs;                /*!< 边界流动变量集合 */
    vector<INT>                 mIndex;            /*!< 方法索引 */
    vector<BoundaryFlowMethod*> bfM;               /*!< 热损失计算方法集合 */
};

#endif /* end if __BOUNDARYFLOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/27/2023      Create file                          */
/*----------------------------------------------------------------------------*/
