/*! \file    UtilOutput.hpp
 *  \brief   用于输出文件的基础工具集。
 *  \details 本文件包含了用于输出文件的基础工具函数和模板。
 *  \author  Shizhe Li
 *  \date    Oct/11/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __UTILOUTPUT_HEADER__
#define __UTILOUTPUT_HEADER__

// 标准头文件
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// OpenCAEPoroX 头文件
#include "OCPConst.hpp"

using namespace std;

extern OCP_INT  CURRENT_RANK; ///< 当前进程的等级（类型：OCP_INT，默认值：无）

/*!
 * \brief 生成一个由特定字符组成的分隔线字符串。
 * \param n 分隔线的长度。
 * \return 返回一个由'-'字符组成的字符串。
 */
template <typename T>
constexpr auto OCP_SEP01(T n){
    return string(n, '-');
}

/*!
 * \brief 生成一个由特定字符组成的分隔线字符串。
 * \param n 分隔线的长度。
 * \return 返回一个由'='字符组成的字符串。
 */
template <typename T>
constexpr auto OCP_SEP02(T n){
    return string(n, '=');
}

/*!
 * \brief 获取格式化的索引字符串。
 * \param i 索引i的字符串表示。
 * \param j 索引j的字符串表示。
 * \param k 索引k的字符串表示。
 * \param s 间隔大小。
 * \return 返回格式化的索引字符串。
 */
string GetIJKformat(const string& i, const string& j, const string& k, const USI& s);

/*!
 * \brief 获取格式化的索引字符串。
 * \param i 索引i的无符号短整型表示。
 * \param j 索引j的无符号短整型表示。
 * \param k 索引k的无符号短整型表示。
 * \param s 间隔大小。
 * \return 返回格式化的索引字符串。
 */
string GetIJKformat(const USI& i, const USI& j, const USI& k, const USI& s);

#endif /* end if __UTILOUTPUT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/11/2022      Create file                          */
/*----------------------------------------------------------------------------*/
