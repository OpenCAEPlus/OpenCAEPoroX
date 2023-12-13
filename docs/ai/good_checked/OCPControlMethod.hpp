/*! \file    ControlMethod.hpp
 *  \brief   ControlMethod 类声明
 *  \author  Shizhe Li
 *  \date    Dec/05/2023
 *
 *  本文件包含ControlMethod类的声明，用于控制求解器方法的使用。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPCONTROLMETHOD_HEADER__
#define __OCPCONTROLMETHOD_HEADER__

// 标准头文件
#include <vector>

// OpenCAEPoroX 头文件
#include "OCPConst.hpp"
#include "ParamControl.hpp"
#include "OCPControlFast.hpp"

using namespace std;

/// \class ControlMethod
/// \brief 控制求解器方法的使用
class ControlMethod
{
    friend class OCPControl;

public:
    /// \brief 设置控制参数
    /// \param CtrlParam 控制参数对象
    void SetCtrlParam(const ParamControl& CtrlParam);

    /// \brief 设置快速控制
    /// \param fCtrl 快速控制对象
    void SetFastControl(const FastControl& fCtrl);

    /// \brief 初始化方法调用顺序
    /// \return OCPNLMethod 非线性方法对象
    OCPNLMethod InitMethod() const;

    /// \brief 切换到主方法
    /// \return OCPNLMethod 非线性方法对象
    OCPNLMethod SwitchMethod() const;

    /// \brief 输出模型和方法信息
    /// \param pl 打印级别
    void OutputInfo(const USI& pl) const;

    /// \brief 获取模型
    /// \return 模型枚举类型
    auto GetModel() const { return model; }

    /// \brief 获取解决方案的类型
    /// \return 方法向量
    auto GetMethod() const { return method; }

    /// \brief 获取第i个线性求解器文件
    /// \param i 索引
    /// \return 线性求解器文件名
    auto GetLsFile(const USI& i) const { return lsFile[i]; }

    /// \brief 获取工作目录名
    /// \return 工作目录名字符串
    auto GetWorkDir() const { return workDir; }

protected:
    /// 工作目录
    string              workDir;

    /// \brief 模型：等温的，热的
    /// \details 定义模型的类型，可以是等温的或热的
    OCPModel            model{ OCPModel::none };

    /// \brief 工作方法的索引
    /// \details 可变成员，用于指示当前工作的方法索引
    mutable USI         wIndex;

    /// \brief 主方法在前，然后依次是预处理方法
    /// \details 包含非线性方法对象的向量
    vector<OCPNLMethod> method;

    /// \brief 线性求解器的文件名
    /// \details 包含线性求解器文件名的向量
    vector<string>      lsFile;
};

#endif /* end if __OCPControlMethod_HEADER__ */

/*----------------------------------------------------------------------------
 *  Brief Change History of This File
 *----------------------------------------------------------------------------
 *
 *  Author              Date             Actions
 *----------------------------------------------------------------------------
 *
 *  Shizhe Li           Dec/05/2023      Create file
 *----------------------------------------------------------------------------
 */
