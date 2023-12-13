/*! 
 *  \file    OCPControlFast.hpp
 *  \brief   快速控制类FastControl的声明
 *  \details 该文件包含了FastControl类的声明，用于处理来自命令行的快捷指令。
 *  \author  Shizhe Li
 *  \date    Dec/05/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPCONTROLFAST_HEADER__
#define __OCPCONTROLFAST_HEADER__

// 标准头文件
#include <vector>
// OpenCAEPoroX头文件
#include "OCPConst.hpp"

using namespace std;

/// \brief 命令行的快捷指令类
/// \details 用于解析命令行参数并设置模拟的控制参数。
class FastControl
{
public:
    /// \brief 构造函数
    /// \param argc 参数数量
    /// \param optset 参数字符串数组
    /// \details 解析命令行参数，设置时间步长和打印级别等
    FastControl(const USI& argc, const char* optset[])
    {
        ifUse = OCP_FALSE;
        timeInit = timeMax = timeMin = -1.0;
        std::stringstream buffer;
        string            tmp;
        string            key;
        string            value;
        for (USI n = 2; n < argc; n++) {
            buffer << optset[n];
            buffer >> tmp;
            string::size_type pos = tmp.find_last_of('=');
            if (pos == string::npos) OCP_ABORT("未知的使用方式！请查看-h");
            key = tmp.substr(0, pos);
            value = tmp.substr(pos + 1, tmp.size() - pos);
            switch (Map_Str2Int(&key[0], key.size())) {
            case Map_Str2Int("method", 6):
                if (value == "FIM") {
                    method = OCPNLMethod::FIM;
                }
                else if (value == "IMPEC") {
                    method = OCPNLMethod::IMPEC;
                }
                else if (value == "AIMc") {
                    method = OCPNLMethod::AIMc;
                }
                else {
                    OCP_ABORT("命令行中的方法参数错误！");
                }
                ifUse = OCP_TRUE;
                if (method == OCPNLMethod::FIM || method == OCPNLMethod::AIMc) {
                    if (timeInit <= 0) timeInit = 1;
                    if (timeMax <= 0) timeMax = 10.0;
                    if (timeMin <= 0) timeMin = 0.1;
                }
                else {
                    if (timeInit <= 0) timeInit = 0.1;
                    if (timeMax <= 0) timeMax = 1.0;
                    if (timeMin <= 0) timeMin = 0.1;
                }
                break;
            case Map_Str2Int("dtInit", 6):
                timeInit = stod(value);
                break;
            case Map_Str2Int("dtMin", 5):
                timeMin = stod(value);
                break;
            case Map_Str2Int("dtMax", 5):
                timeMax = stod(value);
                break;
            case Map_Str2Int("verbose", 7):
                printLevel = OCP_MIN(OCP_MAX(stoi(value), PRINT_NONE), PRINT_ALL);
                break;
            default:
                OCP_ABORT("未知选项：" + key + "   请查看-h");
                break;
            }
            buffer.clear();
        }
    }

public:
    /// \brief 是否使用快速控制
    /// \details 标记是否根据命令行参数启用快速控制
    OCP_BOOL    ifUse{ OCP_FALSE };
    
    /// \brief 数值方法
    /// \details 指定使用的数值方法，可能是IMPEC、FIM或AIM
    OCPNLMethod method;
    
    /// \brief 第一个时间步长的长度
    /// \details 下一个TSTEP开始的第一个时间步长
    OCP_DBL     timeInit;
    
    /// \brief 运行期间的最大时间步长
    /// \details 指定运行过程中使用的最大时间步长
    OCP_DBL     timeMax;
    
    /// \brief 运行期间的最小时间步长
    /// \details 指定运行过程中使用的最小时间步长
    OCP_DBL     timeMin;
    
    /// \brief 打印级别
    /// \details 控制日志打印的详细程度
    USI         printLevel{ 0 };
};

#endif /* end if __OCPControlFast_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/05/2023      Create file                          */
/*----------------------------------------------------------------------------*/
