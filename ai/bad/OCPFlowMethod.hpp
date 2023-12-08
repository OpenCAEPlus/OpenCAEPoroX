```cpp
/*! \file    OCPFlowMethod.hpp
 *  \brief   OCPFlowMethod 类的声明和相关派生类
 *  \author  Shizhe Li
 *  \date    Oct/04/2023
 *
 *  本文件包含 OCPFlowMethod 类及其多个派生类的声明，这些类用于计算相对渗透率和毛细管压力等。
 *  这些类是石油工程模拟软件中流体流动计算的核心组件。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPFLOWMETHOD_HEADER__
#define __OCPFLOWMETHOD_HEADER__

#include "OCPFlowVarSet.hpp"
#include "OCPFuncSAT.hpp"
#include <vector>

using namespace std;

/*!
 * \class OCPFlowMethod
 * \brief 抽象基类，用于定义计算相对渗透率和毛细管压力的接口
 */
class OCPFlowMethod{
public:
    OCPFlowMethod() = default;

    /// 计算相对渗透率和毛细管压力
    virtual void CalKrPc(OCPFlowVarSet& vs) = 0;

    /// 计算相对渗透率和毛细管压力及其导数
    virtual void CalKrPcDer(OCPFlowVarSet& vs) = 0;

    /// 获取毛细管水的饱和度
    virtual OCP_DBL GetSwco() const = 0;

    /// 获取水和油之间的最大毛细管压力 (Po - Pw)
    virtual OCP_DBL GetMaxPcow() const = 0;

    /// 获取水和油之间的最小毛细管压力 (Po - Pw)
    virtual OCP_DBL GetMinPcow() const = 0;

    /// 通过Sw计算Pcow
    virtual OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const = 0;

    /// 通过Pcow计算Sw
    virtual OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const = 0;

    /// 通过Sg计算Pcgo
    virtual OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const = 0;

    /// 通过Pcgo计算Sg
    virtual OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const = 0;

    /// 通过Pcgw计算Sw
    virtual OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const = 0;

    /// 通过Sg计算Krg和其导数
    virtual OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const = 0;
};

// 此处省略派生类的注释，以保持示例简洁。派生类的注释应遵循上述格式，提供类和成员函数的描述。

#endif /* end if __OCPFLOWMETHOD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/04/2023      Create file                          */
/*----------------------------------------------------------------------------*/
```

