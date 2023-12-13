/*! \file    UtilInput.hpp
 *  \brief   本文件提供了用于输入文件的基础工具。
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __UTILINPUT_HEADER__
#define __UTILINPUT_HEADER__

// 标准头文件
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// OpenCAEPoroX 头文件
#include "OCPConst.hpp"

using namespace std;

/// \brief 宏定义：用于参数检查，替换为错误日志。
#define ParamCheck1(exp)                                                               \
    std::cout << exp << " in " << __func__ << "() in " << __LINE__ << " in "           \
              << __FILE__;

/// \brief 用于将字符串映射到整数，以便在输入文件的switch结构中高效匹配关键字。
/// \param mystr 输入的字符串
/// \param len 字符串长度
/// \return 映射后的长整型数
constexpr inline long long Map_Str2Int(const char* mystr, const USI& len){
    long long res = 0;
    long long t   = 100;
    for (USI i = 0; i < len; i++) {
        res += (int)mystr[len - 1 - i] * t;
        t *= 100;
    }
    return res;
}

/// \brief 读取文件时的核心函数。它会捕获下一行有意义的内容，
/// 例如，非空白行或注释，然后去掉一些无用字符，如空格、逗号。
/// 最后，剩余字符串的段将存储在结果中。如果返回OCP_FALSE，表示我们已达到文件末尾。
/// \param ifs 输入文件流
/// \param result 存储处理后的字符串段
/// \param no_slash 是否忽略斜杠
/// \return 是否成功读取到有意义的行
OCP_BOOL ReadLine(ifstream& ifs, vector<string>& result, const OCP_BOOL& no_slash = OCP_TRUE);

/// \brief 用于处理带有星号的表达式，例如
/// m*n  -> <n,...,n> size m ,  m* -> <DEFAULT,..., DEFAULT> size m.
/// \param result 存储处理后的字符串段
void DealDefault(vector<string>& result);

/// \brief 将一系列整数的乘积转换为两个数组。
/// 例如, 8*1  16*2  8*3  16*4  -> obj <8, 16, 8, 16> & val <1, 2, 3, 4>.
/// \param vbuf 包含乘积字符串的向量
/// \param obj 存储对象数量
/// \param region 存储区域值
template <typename T>
void DealData(const vector<string>& vbuf, vector<OCP_USI>& obj, vector<T>& region){
    obj.resize(0);
    region.resize(0);
    for (auto& str : vbuf) {
        auto pos = str.find('*');
        if (pos != string::npos) {
            USI     len = str.size();
            OCP_USI num = stoi(str.substr(0, pos));
            USI     val = stoi(str.substr(pos + 1, len - (pos + 1)));
            obj.push_back(num);
            region.push_back(val);
        }
    }
}

#endif /* end if __UTILINPUT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/
