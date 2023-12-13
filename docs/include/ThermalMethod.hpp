/*! \file    ThermalMethod.hpp
 *  \brief   热法模拟类的声明
 *  \author  Shizhe Li
 *  \date    Nov/10/2022
 *
 *  本文件包含T_FIM类的声明，用于石油工程中的热法模拟。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __THERMALMETHOD_HEADER__
#define __THERMALMETHOD_HEADER__

#include "LinearSystem.hpp"
#include "OCPControl.hpp"
#include "Reservoir.hpp"
#include "UtilOutput.hpp"
#include "UtilTiming.hpp"
#include "OCPTimeRecord.hpp"
#include "OCPNRsuite.hpp"

/*!
 * \class T_FIM
 * \brief 热法模拟类
 *
 * T_FIM类用于模拟热法回收过程中的多相流动和热传导。
 */
class T_FIM
{
public:
    /*!
     * \brief 设置模拟环境
     * \param rs 油藏对象
     * \param ctrl 控制参数对象
     */
    void Setup(Reservoir& rs, const OCPControl& ctrl);

    /*!
     * \brief 初始化油藏
     * \param rs 油藏对象
     */
    void InitReservoir(Reservoir& rs);

    /*!
     * \brief 准备模拟步骤
     * \param rs 油藏对象
     * \param ctrl 控制参数对象
     */
    void Prepare(Reservoir& rs, const OCPControl& ctrl);

    /*!
     * \brief 组装矩阵
     * \param ls 线性系统对象
     * \param rs 油藏对象
     * \param dt 时间步长
     */
    void AssembleMat(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt);

    /*!
     * \brief 解线性系统
     * \param ls 线性系统对象
     * \param rs 油藏对象
     * \param ctrl 控制参数对象
     */
    void SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl);

    /*!
     * \brief 更新属性
     * \param rs 油藏对象
     * \param ctrl 控制参数对象
     * \return 更新是否成功
     */
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);

    /*!
     * \brief 完成牛顿-拉夫森迭代
     * \param rs 油藏对象
     * \param ctrl 控制参数对象
     * \return 迭代是否完成
     */
    OCP_BOOL FinishNR(Reservoir& rs, OCPControl& ctrl);

    /*!
     * \brief 完成模拟步骤
     * \param rs 油藏对象
     * \param ctrl 控制参数对象
     */
    void FinishStep(Reservoir& rs, OCPControl& ctrl);

    /*!
     * \brief 获取牛顿-拉夫森迭代套件
     * \return 牛顿-拉夫森迭代套件对象的常量引用
     */
    const OCPNRsuite& GetNRsuite() const;

    /*!
     * \brief 设置线性求解器工作参数
     * \param w 线性求解器索引
     * \param i 迭代次数
     */
    void SetWorkLS(const USI& w, const USI& i);

    /*!
     * \brief 获取当前线性求解器工作索引
     * \return 线性求解器工作索引
     */
    USI GetWorkLS()const;

protected:
    /*!
     * \brief 分配油藏资源
     * \param rs 油藏对象
     */
    void AllocateReservoir(Reservoir& rs);

    /*!
     * \brief 初始化岩石属性
     * \param bk 岩石体积单元对象
     */
    void InitRock(Bulk& bk) const;

    /*!
     * \brief 计算岩石属性
     * \param bk 岩石体积单元对象
     */
    void CalRock(Bulk& bk) const;

    /*!
     * \brief 初始化闪蒸过程
     * \param bk 岩石体积单元对象
     */
    void InitFlash(Bulk& bk);

    /*!
     * \brief 计算闪蒸过程
     * \param bk 岩石体积单元对象
     */
    void CalFlash(Bulk& bk);

    /*!
     * \brief 传递闪蒸值
     * \param bk 岩石体积单元对象
     * \param n 组件数量
     */
    void PassFlashValue(Bulk& bk, const OCP_USI& n);

    /*!
     * \brief 计算相对渗透率和毛细管压力
     * \param bk 岩石体积单元对象
     */
    void CalKrPc(Bulk& bk) const;

    /*!
     * \brief 更新上一个时间步
     * \param rs 油藏对象
     */
    void UpdateLastTimeStep(Reservoir& rs) const;

    /*!
     * \brief 计算残余
     * \param rs 油藏对象
     * \param dt 时间步长
     * \param initRes0 是否为初始残余
     */
    void CalRes(Reservoir& rs, const OCP_DBL& dt, const OCP_BOOL& initRes0);

    /*!
     * \brief 组装岩石体积单元矩阵
     * \param ls 线性系统对象
     * \param rs 油藏对象
     * \param dt 时间步长
     */
    void AssembleMatBulks(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;

    /*!
     * \brief 组装井矩阵
     * \param ls 线性系统对象
     * \param rs 油藏对象
     * \param dt 时间步长
     */
    void AssembleMatWells(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;

    /*!
     * \brief 获取解决方案
     * \param rs 油藏对象
     * \param u 解向量
     * \param ctrlNR 牛顿-拉夫森控制对象
     */
    void GetSolution(Reservoir& rs, vector<OCP_DBL>& u, const ControlNR& ctrlNR);

    /*!
     * \brief 重置至上一个时间步
     * \param rs 油藏对象
     * \param ctrl 控制参数对象
     */
    void ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl);

protected:
    /// 如果作为其他方法的预处理器
    OCP_BOOL preM{ OCP_FALSE };
    /// 线性求解器方法索引
    USI wls;
    /// 牛顿-拉夫森迭代套件
    OCPNRsuite NR;
};

#endif /* end if __THERMALMETHOD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/
