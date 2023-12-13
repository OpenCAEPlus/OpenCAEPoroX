/*! 
 * \file    OCPConst.hpp
 * \brief   本文件包含了OpenCAEPoroX软件中使用的内建数据类型和常量的定义。
 * \author  Shizhe Li
 * \date    Oct/01/2021
 *
 * 本文件中定义了多种用于模拟过程中的常量、枚举类型以及宏定义，用于控制模拟过程和处理数据。
 * 这些定义对于软件的其他部分至关重要，因为它们提供了一致的参考和错误处理机制。
 */

#ifndef __OPENCAEPORO_CONSTS_HEADER__
#define __OPENCAEPORO_CONSTS_HEADER__

// 标准头文件
#include <iostream>
#include <math.h>

// OpenCAEPoroX头文件
#include "UtilError.hpp"
#include "OCPDataType.hpp"
#include "OCPUnits.hpp"

// 通用错误类型
const int OCP_SUCCESS         = 0;    ///< 成功完成，无错误
const int OCP_ERROR_NUM_INPUT = -1;   ///< 输入参数数量错误
const int OCP_ERROR           = -100; ///< 未识别错误

// 通用常量
const OCP_DBL  TINY         = 1E-8;        ///< 微小常数，用于避免除零等操作
const OCP_DBL  OCP_HUGE     = 1E16;        ///< 巨大常数，用于初始化极大值
const OCP_DBL  PI           = 3.141592653; ///< 圆周率常数
const OCP_BOOL OCP_TRUE     = 1;           ///< 布尔真值
const OCP_BOOL OCP_FALSE    = 0;           ///< 布尔假值

// 控制常量
const OCP_DBL MAX_TIME_STEP     = 365.0; ///< 最大时间步长
const OCP_DBL MIN_TIME_STEP     = 0.01;  ///< 最小时间步长
const OCP_DBL MIN_TIME_CURSTEP  = 1E-8;  ///< 当前步骤的最小时间步长，用于？
const OCP_DBL TIME_STEP_CUT     = 0.5;   ///< 时间步长减小比例
const OCP_DBL TIME_STEP_AMPLIFY = 2.0;   ///< 时间步长放大比例
const OCP_DBL MAX_VOLUME_ERR    = 0.01;  ///< 最大体积误差
const OCP_DBL MAX_DP_LIMIT      = 200;   ///< 最大压力变化限制
const OCP_DBL MAX_DS_LIMIT      = 0.1;   ///< 最大饱和度变化限制
const OCP_DBL TARGET_DP         = 50;    ///< 目标压力变化
const OCP_DBL TARGET_DS         = 0.01;  ///< 目标饱和度变化

// 网格类型枚举
enum class GridType : USI{
    structured,    ///< 结构化网格
    orthogonal,    ///< 正交网格
    corner,        ///< 角点网格
    unstructured,  ///< 非结构化网格
    gmsh           ///< Gmsh网格
};

// 模型类型枚举
enum class OCPModel : USI{
    none,        ///< 无模型
    isothermal,  ///< 等温模型
    thermal      ///< 热模型
};

// 解非线性方程的方法枚举
enum class OCPNLMethod : USI{
    none,    ///< 无方法
    IMPEC,   ///< 隐式方法
    FIM,     ///< 全隐式方法
    AIMc,    ///< 自适应隐式方法
    FIMddm   ///< 分布式全隐式方法
};

// 连接方向枚举
enum class ConnDirect : USI{
    n,    ///< 无方向
    x,    ///< x方向
    y,    ///< y方向
    z,    ///< z方向
    mf,   ///< 矩阵到裂缝
    fm,   ///< 裂缝到矩阵
    xp,   ///< (i,j,k)到(i+1,j,k)，仅限结构化网格
    xm,   ///< (i,j,k)到(i-1,j,k)，仅限结构化网格
    yp,   ///< (i,j,k)到(i,j+1,k)，仅限结构化网格
    ym,   ///< (i,j,k)到(i,j-1,k)，仅限结构化网格
    zp,   ///< (i,j,k)到(i,j,k+1)，仅限结构化网格
    zm,   ///< (i,j,k)到(i,j,k-1)，仅限结构化网格
    usg   ///< 非结构化网格
};

// 油藏状态枚举
enum class ReservoirState : USI{
    bulk_success,          ///< 批量成功
    bulk_negative_P,       ///< 批量负压力
    bulk_negative_T,       ///< 批量负温度
    bulk_negative_N,       ///< 批量负组分摩尔数
    bulk_large_EV,         ///< 批量大体积误差
    bulk_large_CFL,        ///< 批量大CFL条件
    well_success,          ///< 井成功
    well_negative_P,       ///< 井负压力
    well_switch_BHPm,      ///< 井转换到BHP模式
    well_crossflow         ///< 井交叉流
};

// 并行模块常量
const USI MASTER_PROCESS  = 0; ///< 主进程
const OCP_BOOL PRINTINPUT = OCP_TRUE; ///< 是否打印输入

/**
 * \brief 所有子程序的打印级别 -- 不包括DEBUG输出
 */
#define PRINT_NONE              0  /**< 完全沉默：不打印任何内容 */
#define PRINT_MIN               1  /**< 安静：只打印错误和重要警告 */
#define PRINT_SOME              2  /**< 一些：打印不太重要的警告 */
#define PRINT_MORE              4  /**< 更多：打印一些有用的调试信息 */
#define PRINT_MOST              8  /**< 最多：最大打印量，不包括文件 */
#define PRINT_ALL              10  /**< 全部：包括文件在内的所有打印内容 */

/**
 * \brief 定义max, min, abs宏
 */
#define OCP_MAX(a,b) (((a)>(b))?(a):(b))     /**< a和b中较大的一个 */
#define OCP_MIN(a,b) (((a)<(b))?(a):(b))     /**< a和b中较小的一个 */
#define OCP_ABS(a)   (((a)>=0.0)?(a):-(a))   /**< a的绝对值 */

#endif // __OPENCAEPORO_CONSTS_HEADER__

/*----------------------------------------------------------------------------
 * Brief Change History of This File
 *----------------------------------------------------------------------------
 * Author              Date             Actions
 *----------------------------------------------------------------------------
 * Shizhe Li           Oct/01/2021      Create file
 * Chensong Zhang      Oct/15/2021      Format file
 * Chensong Zhang      Oct/27/2021      Unify error check
 * Chensong Zhang      Jan/16/2022      Update Doxygen
 * Chensong Zhang      Sep/21/2022      Add error messages
 *----------------------------------------------------------------------------
 */
