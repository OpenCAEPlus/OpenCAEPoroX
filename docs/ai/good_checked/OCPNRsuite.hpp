/*! 
 * \file    OCPNRsuite.hpp
 * \brief   非线性迭代中使用的数据结构
 * \author  Shizhe Li
 * \date    Oct/30/2021
 *
 *-----------------------------------------------------------------------------------
 * Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 * Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPNRSUITE_HEADER__
#define __OCPNRSUITE_HEADER__

// OpenCAEPoroX header files
#include "Reservoir.hpp"
using namespace std;

/*!
 * \enum OCPNRStateP
 * \brief NR物理检查状态
 */
enum class OCPNRStateP : USI
{
    continueSol, ///< 继续求解
    reset,       ///< 重置到上一个时间步
    resetCut,    ///< 重置到上一个时间步并减少时间
    resetCutCFL, ///< 由于CFL过大重置到上一个时间步并减少时间
};

/*!
 * \enum OCPNRStateC
 * \brief NR收敛检查状态
 */
enum class OCPNRStateC : USI
{
    converge,       ///< 收敛
    continueIter,   ///< 继续迭代
    not_converge,   ///< 未收敛
};

/*!
 * \class OCPNRsuite
 * \brief 非线性求解的NR数据集
 */
class OCPNRsuite
{
public:
    void Setup(const OCP_BOOL& ifthermal, const BulkVarSet& bvs, const OCP_USI& nw, const Domain& domain); ///< 设置使用NR的方法
    void Setup(const BulkVarSet& bvs, const Domain& domain); ///< 设置不使用NR的方法

protected:
    MPI_Comm myComm; ///< 通讯器
    OCP_INT numproc, myrank; ///< 进程数和当前进程排名
    OCP_BOOL ifUseNR{ OCP_FALSE }; ///< 是否使用NR
    OCP_BOOL ifThermal{ OCP_FALSE }; ///< 是否是热模型
    OCP_USI nb; ///< bulk数量
    USI np, nc; ///< 相数量和组分数量

public:
    OCPNRresidual res; ///< 残差

public:
    void InitStep(const BulkVarSet& bvs); ///< 初始化NR步骤
    void CalMaxChangeNR(const Reservoir& rs); ///< 计算热模型的最大变化
    OCP_DBL DP(const OCP_USI& n) const { return dP[n]; } ///< 获取压力变化
    OCP_DBL DN(const OCP_USI& n, const USI& i) const { return dN[n * nc + i]; } ///< 获取组分变化
    OCP_DBL DPBmaxNRc() const { return dPBmaxNR.back(); }; ///< 获取当前最大压力变化(bulk)
    OCP_DBL DSmaxNRc() const { return dSmaxNR.back(); }; ///< 获取当前最大饱和度变化
    const auto& DPBmaxNR() const { return dPBmaxNR; }; ///< 获取所有bulk的最大压力变化
    const auto& DPWmaxNR() const { return dPWmaxNR; }; ///< 获取所有well的最大压力变化
    const auto& DTmaxNR() const { return dTmaxNR; }; ///< 获取所有最大温度变化
    const auto& DNmaxNR() const { return dNmaxNR; }; ///< 获取所有最大组分变化
    const auto& DSmaxNR() const { return dSmaxNR; }; ///< 获取所有最大饱和度变化
    const auto& EVmaxNR() const { return eVmaxNR; }; ///< 获取所有最大体积误差

protected:
    vector<OCP_DBL> lP; ///< 上一个NR步骤的压力
    vector<OCP_DBL> lT; ///< 上一个NR步骤的温度
    vector<OCP_DBL> lN; ///< 上一个NR步骤的组分
    vector<OCP_DBL> lS; ///< 上一个NR步骤的饱和度
    vector<OCP_DBL> dP; ///< NR步骤之间的压力变化
    vector<OCP_DBL> dT; ///< NR步骤之间的温度变化
    vector<OCP_DBL> dN; ///< NR步骤之间的组分变化
    vector<OCP_DBL> dS; ///< NR步骤之间的饱和度变化
    vector<OCP_DBL> dPBmaxNR; ///< 时间步内所有NR步骤的最大压力差(bulk)
    vector<OCP_DBL> dPWmaxNR; ///< 时间步内所有NR步骤的最大压力差(well)
    vector<OCP_DBL> dTmaxNR; ///< 时间步内所有NR步骤的最大温度差(bulk)
    vector<OCP_DBL> dNmaxNR; ///< 时间步内所有NR步骤的最大组分差(bulk)
    vector<OCP_DBL> dSmaxNR; ///< 时间步内所有NR步骤的最大饱和度差(bulk)
    vector<OCP_DBL> eVmaxNR; ///< 时间步内所有NR步骤的最大体积误差(bulk)

public:
    void CalMaxChangeTime(const Reservoir& rs); ///< 计算时间步之间的最大变化
    OCP_DBL DPmaxT() const { return dPmaxT; } ///< 返回时间步的最大压力变化
    OCP_DBL DPBmaxT() const { return dPBmaxT; } ///< 返回时间步的最大压力变化(bulk)
    OCP_DBL DPWmaxT() const { return dPWmaxT; } ///< 返回时间步的最大压力变化(well)
    OCP_DBL DTmaxT() const { return dTmaxT; } ///< 返回时间步的最大温度变化
    OCP_DBL DNmaxT() const { return dNmaxT; } ///< 返回时间步的最大组分变化
    OCP_DBL DSmaxT() const { return dSmaxT; } ///< 返回时间步的最大饱和度变化
    OCP_DBL EVmaxT() const { return eVmaxT; } ///< 返回时间步的最大体积误差

protected:
    OCP_DBL dPmaxT; ///< 当前时间步的最大压力变化
    OCP_DBL dPBmaxT; ///< 当前时间步的最大压力变化(bulk)
    OCP_DBL dPWmaxT; ///< 当前时间步的最大压力变化(well)
    OCP_DBL dTmaxT; ///< 当前时间步的最大温度变化
    OCP_DBL dNmaxT; ///< 当前时间步的最大组分变化
    OCP_DBL dSmaxT; ///< 当前时间步的最大饱和度变化
    OCP_DBL eVmaxT; ///< 当前时间步的体积误差

public:
    void CalCFL(const Reservoir& rs, const OCP_DBL& dt, const OCP_BOOL& ifComm); ///< 计算CFL数
    OCP_DBL GetMaxCFL() const { return maxCFL; } ///< 获取最大CFL数
    const OCP_DBL& GetCFL(const OCP_USI& n, const USI& j) const { return cfl[n * np + j]; } ///< 获取CFL数

protected:
    ReservoirState CheckCFL(const OCP_DBL& cflLim) const; ///< 检查CFL数
    vector<OCP_DBL> cfl; ///< 每个bulk的CFL数
    OCP_DBL maxCFL{ 0 }; ///< 全局最大CFL数
    OCP_DBL maxCFL_loc{ 0 }; ///< 本地最大CFL数

public:
    OCP_BOOL CheckPhysical(Reservoir& rs, const initializer_list<string>& il) const; ///< 检查物理状态
    auto GetWorkState() const { return workState; } ///< 获取工作状态

protected:
    mutable OCPNRStateP workState; ///< 工作状态

public:
    void InitIter(); ///< 初始化迭代器
    void UpdateIter(const USI& lsIter); ///< 更新迭代器
    void ResetIter(); ///< 重置迭代器
    auto GetIterNR() const { return iterNR; } ///< 获取NR迭代次数
    auto GetIterLS() const { return iterLS; } ///< 获取LS迭代次数
    auto GetIterNRw() const { return iterNRw; } ///< 获取浪费的NR迭代次数
    auto GetIterLSw() const { return iterLSw; } ///< 获取浪费的LS迭代次数
    const auto& GetIterNRLS() const { return iterNRLS; }; ///< 获取每个NR步骤的LS迭代次数

protected:
    USI iterNR{ 0 }; ///< 时间步内的Newton-Raphson迭代次数
    USI iterLS{ 0 }; ///< 时间步内的线性求解器迭代次数
    USI iterNRw{ 0 }; ///< 时间步内浪费的Newton-Raphson迭代次数
    USI iterLSw{ 0 }; ///< 时间步内浪费的线性求解器迭代次数
    vector<USI> iterNRLS; ///< 每个NR步骤的线性求解器迭代次数
};

#endif  /* end if __OCPNRSUITE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/30/2021      Create file                          */
/*----------------------------------------------------------------------------*/
