/** 
 * \file    UtilTiming.hpp
 * \brief   本文件包含了计时功能的类和方法声明，用于测量程序运行的墙钟时间和CPU周期。
 * \author  Chensong Zhang
 * \date    Sep/24/2019
 *
 * \note    本文件由FASP++团队维护，遵循GNU Lesser General Public License 3.0或更晚版本发布。
 */

#ifndef __TIMING_HEADER__ // 允许多次包含
#define __TIMING_HEADER__ /**< 表示timing.hxx已被包含 */

typedef unsigned long long uint64; ///< 无符号长长整型

#include <chrono> // 高分辨率CPU时间
#include <iomanip>
#include <iostream>
#include <string> // 输出字符串

// 时间单位定义
const double CLOCK_USE_SEC = 5000;   ///< 以秒为单位显示时钟时间
const double CLOCK_USE_MIN = 200000; ///< 以分钟为单位显示时钟时间

/*! 
 * \class   GetWallTime
 * \brief   用于获取以毫秒为单位的经过的墙钟时间
 * 
 * 本类通过读取当前的墙钟时间，并返回从start()到stop()的持续时间。
 */
class GetWallTime
{
private:
    std::chrono::steady_clock::time_point timeStamp; ///< 当前CPU时间点

public:
#if defined(_CONSOLE) || defined(_WIN32) || defined(_WIN64)
    // Windows文件系统
    /// 开始计时器
    inline void Start() { timeStamp = std::chrono::steady_clock::now(); }
    /// 停止计时器并返回从start()到现在的持续时间（毫秒）
    inline double Stop() const
    {
        auto elapsedTime = std::chrono::steady_clock::now() - timeStamp;
        return std::chrono::duration<double, std::milli>(elapsedTime).count();
    }
#else
    /// 开始计时器
    __inline__ void Start() { timeStamp = std::chrono::steady_clock::now(); }
    /// 停止计时器并返回从start()到现在的持续时间（毫秒）
    __inline__ double Stop() const
    {
        auto elapsedTime = std::chrono::steady_clock::now() - timeStamp;
        return std::chrono::duration<double, std::milli>(elapsedTime).count();
    }
#endif
    /// 停止计时器并打印出持续时间
    void StopInfo(const std::string& info, std::ostream& out = std::cout) const;
};

#endif /* 结束 __TIMING_HEADER__ 的条件编译块 */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Sep/26/2022      Do not use RTC for x86               */
/*----------------------------------------------------------------------------*/