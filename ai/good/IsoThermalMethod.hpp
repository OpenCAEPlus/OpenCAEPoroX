/*! \file    IsothermalMethod.hpp
 *  \brief   本文件包含了OpenCAEPoroX中流体部分解决方案方法的声明。
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __ISOTHERMALMETHOD_HEADER__
#define __ISOTHERMALMETHOD_HEADER__

// OpenCAEPoroX头文件
#include "LinearSystem.hpp"
#include "OCPControl.hpp"
#include "Reservoir.hpp"
#include "UtilOutput.hpp"
#include "UtilTiming.hpp"
#include "OCPTimeRecord.hpp"
#include "OCPNRsuite.hpp"

/*!
 * \class IsothermalMethod
 * \brief 本类是等温方法的基类。
 */
class IsothermalMethod{
public:
    /// 计算岩石
    void CalRock(Bulk& bk) const;
    /// 获取NRsuite
    const OCPNRsuite& GetNRsuite() const { return NR; }
    virtual void ExchangeSolutionP(Reservoir& rs) const;
    virtual void ExchangeSolutionNi(Reservoir& rs) const;
    void SetPreMethod(const OCP_BOOL& flag) { preM = flag; }
    void SetWorkLS(const USI& w, const USI& i);
    USI  GetWorkLS()const { return wls; }
protected:
    /// 如果用作其他方法的预处理器
    OCP_BOOL        preM{ OCP_FALSE }; ///< 预处理器标志，默认为假。
    /// 牛顿-拉夫森迭代套件
    OCPNRsuite      NR; ///< NR迭代套件对象。
    /// 线性求解器方法索引
    USI             wls; ///< 工作线性求解器索引。
};

/*!
 * \class IsoT_IMPEC
 * \brief IsoT_IMPEC是隐式压力显式饱和度(IMPEC)方法的类。
 */
class IsoT_IMPEC : virtual public IsothermalMethod{
public:
    /// 设置IMPEC
    void Setup(Reservoir& rs, const OCPControl& ctrl);
    /// 初始化储层
    void InitReservoir(Reservoir& rs) const;
    /// 准备组装矩阵
    void Prepare(Reservoir& rs, OCPControl& ctrl);
    /// 组装矩阵
    void AssembleMat(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// 解线性系统
    void SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl);
    /// 更新流体属性
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);
    /// 确定NR迭代是否完成
    OCP_BOOL FinishNR(const Reservoir& rs);
    void     FinishStep(Reservoir& rs, OCPControl& ctrl);
protected:
    /// 用Sj进行闪蒸并计算FIM所需值
    void InitFlash(Bulk& bk) const;
    /// 计算FIM所需的相对渗透率和毛细管压力
    void CalKrPc(Bulk& bk) const;
    /// 从闪蒸到体积传递FIM所需值
    void PassFlashValue(Bulk& bk, const OCP_USI& n) const;
private:
    /// 为储层分配内存
    void AllocateReservoir(Reservoir& rs);
    /// 用Ni进行闪蒸并计算FIM所需值
    void CalFlash(Bulk& bk);
    /// 计算体积和井之间的通量
    void CalFlux(Reservoir& rs) const;
    /// 计算体积之间的通量
    void CalBulkFlux(Reservoir& rs) const;
    /// 根据质量守恒更新每个体积的摩尔组成
    void MassConserve(Reservoir& rs, const OCP_DBL& dt) const;
    /// 为体积组装线性系统
    void AssembleMatBulks(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// 为井组装线性系统
    void AssembleMatWells(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// 解线性系统后更新P, Pj, BHP
    void GetSolution(Reservoir& rs, vector<OCP_DBL>& u);
    /// 将变量重置为上一时间步
    void ResetToLastTimeStep01(Reservoir& rs, OCPControl& ctrl);
    void ResetToLastTimeStep02(Reservoir& rs, OCPControl& ctrl);
    /// 更新FIM的上一步值
    void UpdateLastTimeStep(Reservoir& rs) const;
};

// 更多类定义...

#endif /* end if __ISOTHERMALMETHOD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/
