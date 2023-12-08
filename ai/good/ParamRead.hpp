/*! 
 *  \file    ParamRead.hpp
 *  \brief   本文件包含ParamRead类的声明
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *  \note    OpenCAEPoroX中使用的参数大多与SLB的Eclipse兼容，
 *           但它也有一些自己的规则以便于使用。它是可扩展的和用户友好的。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team.
 *  All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARAMREAD_HEADER__
#define __PARAMREAD_HEADER__

// 标准头文件
#include <fstream>
#include <iostream>
#include <string>

// OpenCAEPoroX头文件
#include "ParamControl.hpp"
#include "ParamOutput.hpp"
#include "ParamReservoir.hpp"
#include "ParamWell.hpp"

using namespace std;

/*!
 *  \class  ParamRead
 *  \brief  ParamRead类用于从输入文件中读取参数的预处理单元。
 */
class ParamRead
{
public:
    string   inputFile; ///< 输入文件及其路径（绝对或相对）。
    string   workDir;   ///< 当前的工作目录。
    string   fileName;  ///< 输入文件的文件名。

    // ParamRead的主要工作负载：读取储层参数、井参数、
    // 控制参数和输出参数。这些工作负载分别分配给以下类。
    ParamReservoir paramRs;      ///< 读取储层参数的类实例。
    ParamWell      paramWell;    ///< 读取井参数的类实例。
    ParamControl   paramControl; ///< 读取控制参数的类实例。
    ParamOutput    paramOutput;  ///< 读取输出参数的类实例。

public:
    /// 从完整文件路径获取当前工作目录和输入文件名。
    void GetDirAndName();

    /// 初始化参数读取过程。
    void Init();

    /*!
     *  \brief  通用接口用于读取输入数据。
     *  \param  file 输入文件的字符串。
     */
    void ReadInputFile(const string& file);

    /// 读取输入文件。
    void ReadFile(const string& file);

    /*!
     *  \brief  处理INCLUDE关键字，它包含其他输入文件。
     *  \param  ifs 输入文件流的引用。
     */
    void ReadINCLUDE(ifstream& ifs);

    /// 检查参数是否包含错误。
    void CheckParam();
};

#endif /* end if __PARAMREAD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/16/2021      Format file                          */
/*----------------------------------------------------------------------------*/
