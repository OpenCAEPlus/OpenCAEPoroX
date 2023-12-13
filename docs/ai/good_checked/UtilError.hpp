/** 
 *  \file    UtilError.hpp 
 *  \brief   本文件包含了错误日志和警告消息的宏定义，用于在开发中快速定位和记录问题。
 *  \author  Ronghong Fan
 *  \date    Nov/01/2019
 *
 *  本文件定义了一系列宏来帮助开发者打印错误信息、警告、断言以及重要信息。
 *  它还包括了在调试模式下的断言检查。
 */

#ifndef __ERRORLOG_HXX__ /*-- 防止多重包含 --*/
#define __ERRORLOG_HXX__ /**< 表明ErrorLog.hxx已被包含 */

// 标准头文件
#include <iomanip>
#include <iostream>
#include <sstream>

/**
 *  \brief 打印出函数名、文件名和行号。
 */
#define OCP_LOCATION                                                                   \
    "\n    --> function: " << __PRETTY_FUNCTION__                                      \
                           << "\n    --> file:     " << __FILE__ << "::" << __LINE__

/**
 *  \brief 记录错误消息。
 *  \param msg 用户定义的错误消息。
 *  使用do-while结构以允许宏以分号结束。
 */
#define OCP_MESSAGE(msg)                                                               \
    do {                                                                               \
        std::ostringstream info;                                                       \
        info << std::setprecision(16);                                                 \
        info << msg << OCP_LOCATION << '\n';                                           \
        std::cerr << info.str().c_str();                                               \
    } while (false)

/**
 *  \brief 记录警告消息。
 *  \param msg 用户定义的警告消息。
 *  使用do-while结构以允许宏以分号结束。
 */
#define OCP_WARNING(msg)                                                               \
    do {                                                                               \
        OCP_MESSAGE("### WARNING: " << (msg));                                         \
    } while (false)

/**
 *  \brief 在发生严重错误时中止程序。
 *  \param msg 用户定义的中止消息。
 *  使用do-while结构以允许宏以分号结束。
 */
#define OCP_ABORT(msg)                                                                 \
    do {                                                                               \
        OCP_MESSAGE("### ABORT: " << (msg));                                           \
        std::abort();                                                                  \
    } while (false)

/**
 *  \brief 在DEBUG模式下断言条件并记录用户消息。
 *  \param cond 检查的条件。
 *  \param msg 用户定义的错误消息。
 *  使用do-while结构以允许宏以分号结束。
 */
#ifndef DEBUG
#define OCP_ASSERT(cond, msg)                                                          \
    do {                                                                               \
    } while (false)
#else
#define OCP_ASSERT(cond, msg)                                                          \
    do {                                                                               \
        if (!(cond)) {                                                                 \
            OCP_MESSAGE("### ASSERT: " << (msg) << " (" << #cond << ")");              \
            std::abort();                                                              \
        }                                                                              \
    } while (false)
#endif

/**
 *  \brief 打印重要信息。
 *  \param msg 用户定义的信息。
 *  使用do-while结构以允许宏以分号结束。
 */
#define OCP_INFO(msg)                                                               \
    do {                                                                            \
        std::cout << "###INFO: " << msg << "!\n";                                   \
    } while (false)                                                                    

/**
 *  \brief 打印函数名称。
 */
#ifndef OCPFUNCNAME
#define OCP_FUNCNAME                                                                   \
    do {                                                                               \
    } while (false)
#else
#define OCP_FUNCNAME                                                                   \
    do {                                                                               \
        std::cout << __FUNCTION__ << std::endl;                                        \
    } while (false)
#endif

#endif /* end if for __ERRORLOG_HXX__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Nov/15/2021      Test DEBUG mode                      */
/*----------------------------------------------------------------------------*/