/*! \file    OCPControl.hpp
 *  \brief   OCPControl类的声明，用于管理和控制模拟过程中的各种参数。
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPCONTROL_HEADER__
#define __OCPCONTROL_HEADER__

// 标准头文件
#include <vector>

// OpenCAEPoroX头文件
#include "OCPConst.hpp"
#include "ParamControl.hpp"
#include "OCPControlTime.hpp"
#include "OCPControlNR.hpp"
#include "OCPControlMethod.hpp"

using namespace std;

/// \class OCPControl
/// \brief 用于解决控制的所有参数。
class OCPControl
{
public:
    /// \brief 输入控制参数。
    /// \param CtrlParam 控制参数对象。
    void InputParam(const ParamControl& CtrlParam);

    /// \brief 设置通信。
    /// \param domain 域对象。
    void Setup(const Domain& domain);

    /// \brief 应用第i时间步的控制。
    /// \param i 时间步索引。
    /// \param wellOptChange_loc 井优化更改的局部布尔标志。
    void ApplyControl(const USI& i, const OCP_BOOL& wellOptChange_loc);

    /// \brief 检查是否收敛。
    /// \param NRs NR套件对象。
    /// \param il 初始化列表，包含字符串。
    /// \return 返回NR状态。
    OCPNRStateC CheckConverge(const OCPNRsuite& NRs, const initializer_list<string>& il);

    /// \brief 计算下一个时间步。
    /// \param NRs NR套件对象。
    /// \param il 初始化列表，包含字符串。
    void CalNextTimeStep(const OCPNRsuite& NRs, const initializer_list<string>& il);

public:
    MPI_Comm         myComm; ///< MPI通信器。
    OCP_INT          numproc; ///< 进程数量。
    OCP_INT          myrank; ///< 当前进程的排名。

public:
    /// \brief 打印级别。
    USI              printLevel{0};
    /// \brief 时间控制对象。
    ControlTime      time;
    /// \brief NR控制对象。
    ControlNR        NR;
    /// \brief 求解方法控制对象。
    ControlMethod    SM;
    /// \brief 停止模拟标志。
    OCP_BOOL         StopSim{ OCP_FALSE };
    /// \brief 停止时间。
    OCP_DBL          MaxSimTime{ 1E20 };

public:
    /// \brief 获取工作目录名称。
    /// \return 返回工作目录字符串。
    auto GetWorkDir() const { return workDir; }

    /// \brief 获取OCP文件名称。
    /// \return 返回OCP文件字符串。
    auto GetOCPFile() const { return ocpFile; }

    /// \brief 设置快速控制。
    /// \param argc 参数数量。
    /// \param optset 选项集合。
    void SetupFastControl(const USI& argc, const char* optset[]);

    /// \brief 输出模型方法信息。
    void OutputModelMethodInfo() const;

protected:
    /// \brief 当前工作目录。
    string           workDir; ///< 工作目录字符串。
    /// \brief 当前文件名。
    string           ocpFile; ///< OCP文件字符串。
};

#endif /* end if __OCPControl_HEADER__ */

/*----------------------------------------------------------------------------
 *  Brief Change History of This File
 *----------------------------------------------------------------------------
 *
 *  Author              Date             Actions
 *
 *----------------------------------------------------------------------------
 *
 *  Shizhe Li           Oct/01/2021      Create file
 *  Chensong Zhang      Oct/15/2021      Format file
 *  Chensong Zhang      Jan/08/2022      Update Doxygen
 *
 *----------------------------------------------------------------------------
 */