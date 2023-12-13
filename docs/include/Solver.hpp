/*! \file    Solver.hpp
 *  \brief   求解器类的声明
 *  \author  Shizhe Li
 *  \date    Oct/21/2021
 *
 *  本文件包含Solver类的声明，用于实现油藏模拟的整体求解方法。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// OpenCAEPoroX头文件
#include "IsothermalSolver.hpp"
#include "ThermalSolver.hpp"
#include "OCPOutput.hpp"

#ifndef __SOLVER_HEADER__
#define __SOLVER_HEADER__

/*!
 * \class Solver
 * \brief 求解器类，用于油藏模拟的整体求解方法。
 */
class Solver
{
public:
    /*!
     * \brief 配置求解器
     * \param rs 油藏对象
     * \param ctrl 控制参数
     */
    void Setup(Reservoir& rs, const OCPControl& ctrl);

    /*!
     * \brief 初始化油藏
     * \param rs 油藏对象
     */
    void InitReservoir(Reservoir& rs);

    /*!
     * \brief 开始模拟
     * \param rs 油藏对象
     * \param ctrl 控制参数
     * \param output 输出对象
     */
    void RunSimulation(Reservoir& rs, OCPControl& ctrl, OCPOutput& output);

private:
    /*!
     * \brief 通用API，执行一个模拟步骤
     * \param rs 油藏对象
     * \param ctrl 控制参数
     * \return OCPNRsuite 对象，包含模拟步骤的结果
     */
    const OCPNRsuite& GoOneStep(Reservoir& rs, OCPControl& ctrl);

    /*!
     * \brief 配置等温模型的求解器
     * \param rs 油藏对象
     * \param ctrl 控制参数
     */
    void SetupIsoT(Reservoir& rs, const OCPControl& ctrl);

    /*!
     * \brief 对等温模型执行一个时间步骤
     * \param rs 油藏对象
     * \param ctrl 控制参数
     * \return OCPNRsuite 对象，包含模拟步骤的结果
     */
    const OCPNRsuite& GoOneStepIsoT(Reservoir& rs, OCPControl& ctrl);

    /*!
     * \brief 配置热模型的求解器
     * \param rs 油藏对象
     * \param ctrl 控制参数
     */
    void SetupT(Reservoir& rs, const OCPControl& ctrl);

    /*!
     * \brief 对热模型执行一个时间步骤
     * \param rs 油藏对象
     * \param ctrl 控制参数
     * \return OCPNRsuite 对象，包含模拟步骤的结果
     */
    const OCPNRsuite& GoOneStepT(Reservoir& rs, OCPControl& ctrl);

protected:
    /*!
     * \var model
     * \brief 模型类型
     * \details 默认值为OCPModel::none，表示没有指定模型类型
     */
    OCPModel model{ OCPModel::none };

    /*!
     * \var IsoTSolver
     * \brief 等温求解器
     * \details 用于固定温度的等温模型求解
     */
    IsothermalSolver IsoTSolver;

    /*!
     * \var TSolver
     * \brief 热求解器
     * \details 用于变温度的热模型求解
     */
    ThermalSolver TSolver;
};

#endif /* end if __SOLVER_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Shizhe Li           Oct/21/2021      Change from OCPMethod to Solver      */
/*  Chensong Zhang      Oct/27/2021      Rearrange and add comments           */
/*----------------------------------------------------------------------------*/