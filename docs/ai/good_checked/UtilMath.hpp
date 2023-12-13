/*! \file    UtilMath.hpp
 *  \brief   数学实用工具类声明
 *  \author  Shizhe Li
 *  \date    Jul/23/2023
 *
 *  本文件包含了数学计算相关工具的类声明，主要用于开发石油模拟软件中的数学计算功能。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __UTILMATH_HEADER__
#define __UTILMATH_HEADER__

#include "OCPConst.hpp"
#include <vector>
#include <algorithm>

using namespace std;

/**
 *  \brief  计算立方根的函数。
 *  
 *  \param a        立方方程的二次项系数。
 *  \param b        立方方程的一次项系数。
 *  \param c        立方方程的常数项。
 *  \param NTflag   是否使用牛顿迭代法的标志。
 *  \param z        存储计算出的立方根的向量。
 *  \return         无返回值。
 */
USI CubicRoot(const OCP_DBL& a, const OCP_DBL& b, const OCP_DBL& c, const OCP_BOOL& NTflag, vector<OCP_DBL>& z);

/**
 *  \brief  使用牛顿迭代法计算立方根。
 *  
 *  \param root     计算得出的立方根。
 *  \param a        立方方程的二次项系数。
 *  \param b        立方方程的一次项系数。
 *  \param c        立方方程的常数项。
 *  \return         无返回值。
 */
void NTcubicroot(OCP_DBL& root, const OCP_DBL& a, const OCP_DBL& b, const OCP_DBL& c);

/**
 *  \brief  计算数值的符号函数。
 *  
 *  \param d        需要计算符号的数值。
 *  \return         返回数值的符号，类型为OCP_DBL。
 */
OCP_DBL signD(const OCP_DBL& d);

/**
 *  \brief  Kronecker delta函数。
 *  
 *  \param i        第一个索引。
 *  \param j        第二个索引。
 *  \return         如果i等于j，返回1.0，否则返回0.0。
 */
OCP_DBL delta(const USI& i, const USI& j);

#endif /* end if __UTILMATH_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/23/2023      Create file                          */
/*----------------------------------------------------------------------------*/