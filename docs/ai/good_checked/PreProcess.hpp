/*! \file    PreProcess.hpp
 *  \brief   预处理模块的头文件，用于OpenCAEPoroX模拟器的输入网格参数处理、网格划分和域的设置
 *  \author  Shizhe Li
 *  \date    Feb/15/2023
 *
 *  本文件包含PreProcess类的声明，该类负责OpenCAEPoroX模拟器中的预处理步骤。
 *  其中包括读取输入文件、网格和井参数的处理、利用Parmetis进行网格划分以及域的定义。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PREPROCESS_HEADER__
#define __PREPROCESS_HEADER__

#include "OCPConst.hpp"
#include "ParamRead.hpp"
#include "PreParamGridWell.hpp"
#include "Partition.hpp"
#include "Domain.hpp"
#include "UtilTiming.hpp"
#include "OCPTimeRecord.hpp"

using namespace std;

/// \class PreProcess
/// \brief 预处理类，负责输入网格参数的读取、网格划分以及域的定义。
class PreProcess
{
    friend class OpenCAEPoroX; ///< 声明OpenCAEPoroX为友元类
    friend class Reservoir;    ///< 声明Reservoir为友元类

public:
    /// \brief 构造函数
    /// \param myFile 输入文件的路径（绝对或相对）
    /// \param myRank 当前进程的MPI等级
    /// \param comm MPI通信器
    PreProcess(const string& myFile, const OCP_INT& myRank, MPI_Comm comm);

protected:
    /// \brief 读取输入文件
    /// \param myFile 输入文件的路径（绝对或相对）
    void GetFile(const string& myFile);

protected:
    string inputFile;    ///< 输入文件及其路径
    string workdir;      ///< 当前工作目录
    string filename;     ///< 输入文件的文件名
    PreParamGridWell preParamGridWell; ///< 网格和井的参数
    Partition        partition;        ///< 使用Parmetis进行的网格划分
    Domain           domain;           ///< 定义的域
};

#endif /* end if __PREPROCESS_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Feb/15/2023      Create file                          */
/*----------------------------------------------------------------------------*/
