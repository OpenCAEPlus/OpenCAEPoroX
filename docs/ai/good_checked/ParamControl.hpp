/*!
 *  \file    ParamControl.hpp
 *  \brief   ParamControl 类的声明
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *  本文件包含 ParamControl 类的声明，该类用于控制模拟过程中的参数设置，
 *  如离散方法的选择，线性求解器文件的使用，时间步长的变化等。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARAMCONTROL_HEADER__
#define __PARAMCONTROL_HEADER__

// 标准头文件
#include <cassert>
#include <fstream>
#include <vector>

// OpenCAEPoroX 头文件
#include "OCPConst.hpp"
#include "UtilInput.hpp"
#include "UtilOutput.hpp"

/// \brief Tuning 是一组控制参数，包含三个主要部分：
/// 1. 时间步进控制。
/// 2. 时间截断与收敛控制。
/// 3. 牛顿迭代和线性迭代控制。
/// 但是大多数参数尚未应用到当前程序。
typedef vector<vector<OCP_DBL>> TUNING;

/// \brief TuningPair 类用于存储一对时间和Tuning参数。
class TuningPair
{
public:
    /// \brief 构造函数
    /// \param t 时间
    /// \param tmp Tuning 参数
    TuningPair(const USI& t, const TUNING& tmp)
        : d(t)
        , Tuning(tmp){};

    USI    d; ///< 时间
    TUNING Tuning; ///< Tuning 参数
};

/// \brief ParamControl 类包含模拟控制相关的参数，例如使用哪种离散方法，
/// 使用哪个线性求解器文件，时间步长如何变化等。
class ParamControl
{
public:
    /// \brief 模型：热的或等温的。
    OCPModel           model{OCPModel::isothermal }; ///< 默认为等温模型

    /// \brief 当前工作目录。
    string             workDir; ///< 工作目录

    /// \brief 当前文件名。
    string             fileName; ///< 文件名

    /// \brief 流体方程的离散方法，求解器和预处理器。
    vector<string>     method{ "FIM" }; ///< 默认方法为 FIM

    /// \brief 方法的线性求解器输入文件。
    vector<string>     lsFile{ "bsr.fasp" }; ///< 默认线性求解器文件为 bsr.fasp

    /// \brief Tuning 集合。
    vector<TuningPair> tuning_T; ///< Tuning 对集合

    /// \brief Tuning 参数。
    TUNING             tuning; ///< Tuning 参数

    /// \brief 最大模拟时间。
    OCP_DBL            MaxSimTime{ 1E20 }; ///< 默认最大模拟时间

    /// \brief 关键时间记录模拟过程中的重要时间点，在这些时间点，模拟过程应该被仔细处理，
    /// 例如边界条件可能会改变。
    vector<OCP_DBL> criticalTime; ///< 关键时间点集合

    /// \brief 为参数分配默认值。
    void Init(string& indir, string& inFileName);

    /// \brief 初始化关键时间。
    void InitTime() { criticalTime.push_back(0); };

    /// \brief 确定默认的 Tuning 参数。
    void InitTuning();

    /// \brief 输入关键字：METHOD。
    void InputMETHOD(ifstream& ifs);

    /// \brief 输入最大模拟时间。
    void InpuMaxSimTime(ifstream& ifs);

    /// \brief 输入关键字：TUNING。
    void InputTUNING(ifstream& ifs);

    /// \brief 显示 Tuning 参数。
    void DisplayTuning() const;
};

#endif /* end if __ParamControl_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/