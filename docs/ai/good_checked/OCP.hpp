/*! \file    OCP.hpp
 *  \brief   OpenCAEPoroX模拟器的主头文件
 *  \author  Shizhe Li
 *  \date    2021年10月1日
 *
 *-----------------------------------------------------------------------------------
 *  版权所有 (C) 2021至今 OpenCAEPoroX团队。保留所有权利。
 *  在GNU Lesser General Public License 3.0或更高版本的条款下发布。
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCP_HEADER__
#define __OCP_HEADER__

// OpenCAEPoroX头文件
#include "OCPControl.hpp"
#include "OCPOutput.hpp"
#include "PreProcess.hpp"
#include "ParamRead.hpp"
#include "Reservoir.hpp"
#include "Solver.hpp"
#include "OCPConst.hpp"
#include "UtilTiming.hpp"

#define OCPVersion "0.5.0" ///< 软件版本标签，用于git

using namespace std;

/// OpenCAEPoroX模拟器的顶层数据结构
class OpenCAEPoroX {
public:
    /// 输出OpenCAEPoroX版本信息
    void PrintVersion() const
    {
        cout << "=========================================" << endl
             << "OpenCAEPoroX版本-" << OCPVersion << endl
             << "=========================================" << endl;
    };

    /// 提供至少一个输入数据的InputFileName
    void PrintUsage(string cmdname) const
    {
        cout << "用法：" << endl
             << "  " << cmdname << " <InputFileName> [<options>]" << endl
             << endl;
        cout << "一个简单的例子是在默认设置下解决SPE1 Case A" << endl
             << "  " << cmdname << " examples/spe1a/spe1a.data" << endl;
        cout << endl
             << "另一个例子是使用FIM解决相同的问题" << endl
             << "  " << cmdname
             << " examples/spe1a/spe1a.data method=FIM dtInit=1 dtMax=10 dtMin=0.1"
             << endl
             << endl;
        cout << "您可以在输入文件后传递可选的命令行参数：" << endl
             << "     method = 使用的解决方法" << endl
             << "     dtInit = 初始时间步长" << endl
             << "      dtMax = 最大时间步长" << endl
             << "      dtMin = 最小时间步长" << endl
             << "    verbose = 屏幕上的打印级别" << endl
             << endl;
        cout << "注意：" << endl
             << "  - 只有设置了`method`，其他选项才会生效；" << endl
             << "  - 这些命令行选项将覆盖输入文件中的选项；" << endl
             << "  - 如果未设置(dtInit, dtMax, dtMin)，将使用默认值。"
             << endl
             << endl;
    }

    /// 从输入文件中读取参数
    void InputDistParam(const string& filename, PreProcess& prepro, const OCP_INT& myRank);

    /// 根据内部结构设置油藏
    void SetupSimulator(const USI& argc, const char* options[]);

    /// 初始化油藏或获取初始状态
    void InitReservoir();

    /// 运行动态模拟
    void RunSimulation();

    /// 输出后处理所需的信息
    void OutputResults() const;

protected:
    /// 输出主进程所需的信息
    void OutputTimeMain(streambuf* mysb) const;

    /// 输出每个进程所需的信息
    void OutputTimeProcess() const;

protected:
    /// 油藏的核心属性
    Reservoir reservoir;

    /// 包含离散方法和线性系统求解器
    Solver solver;

    /// 控制类处理算法参数和时间步长
    OCPControl control;

    /// 输出类处理程序的输出级别
    OCPOutput output;
};

#endif /* end if __OCP_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  此文件的简要更改历史                                                     */
/*----------------------------------------------------------------------------*/
/*  作者              日期             操作                                   */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li         2021年10月1日      创建文件                              */
/*  Chensong Zhang    2021年10月15日     格式化文件                            */
/*  Chensong Zhang    2022年1月8日       新的标签信息                          */
/*  Chensong Zhang    2022年9月21日      添加PrintUsage                        */
/*----------------------------------------------------------------------------*/

