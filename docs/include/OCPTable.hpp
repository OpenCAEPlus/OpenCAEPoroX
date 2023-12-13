/*! 
 * \file    OCPTable.hpp
 * \brief   OCPTable 类声明
 * \author  Shizhe Li
 * \date    Oct/01/2021
 *
 * 本文件包含了 OCPTable 类和 OCPTable_Interpolate 类的声明，这些类用于处理 OpenCAEPoroX 中
 * 与表格相关的所有事务，例如 PVT 表格和饱和度表格。
 *
 *-----------------------------------------------------------------------------------
 * Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 * Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPTable_HEADER__
#define __OCPTable_HEADER__

// 标准头文件
#include <iostream>
#include <vector>

// OpenCAEPoroX 头文件
#include "OCPConst.hpp"
#include "ParamReservoir.hpp"

using namespace std;

/// OCPTable_Interpolate 是一个用于表格插值的基类
class OCPTable_Interpolate{
public:
    /// 默认构造函数
    OCPTable_Interpolate() = default;

    /// 纯虚函数，用于表格插值
    /// \param table 二维数据表
    virtual void Interpolate(const vector<vector<OCP_DBL>>& table) = 0;
};

/// OCPTable 是一个表格类，用于处理 OpenCAEPoroX 中的表格，如 PVT 表格和饱和度表格。
class OCPTable{
public:
    /// 默认构造函数
    OCPTable() = default;

    /// 通过现有数据构造表格
    /// \param src 二维数据表
    OCPTable(const vector<vector<OCP_DBL>>& src);

    /// 从现有数据设置表格
    /// \param src 二维数据表
    void Setup(const vector<vector<OCP_DBL>>& src);

    /// 判断表格是否为空
    /// \return 表格是否为空的布尔值
    OCP_BOOL IsEmpty() const { return data.empty(); }

    /// 获取表格的列数
    /// \return 表格的列数
    USI GetColNum() const { return nCol; }

    /// 获取表格中的指定列
    /// \param j 列索引
    /// \return 指定列的引用
    vector<OCP_DBL>& GetCol(const USI& j) { return data[j]; }
    const vector<OCP_DBL>& GetCol(const USI& j) const { return data[j]; }

    /// 获取表格中的指定值
    /// \param row 行索引
    /// \param col 列索引
    /// \return 指定位置的值
    const auto& GetVal(const USI& row, const USI& col) { return data[col][row]; }

    /// 在表格中插值，评估所有列并返回斜率
    /// \param j 指定列索引
    /// \param val 插值的值
    /// \param outdata 输出数据
    /// \param slope 斜率
    /// \return 评估状态
    USI Eval_All(const USI& j, const OCP_DBL& val, vector<OCP_DBL>& outdata, vector<OCP_DBL>& slope) const;

    /// 在表格中插值，评估所有列
    /// \param j 指定列索引
    /// \param val 插值的值
    /// \param outdata 输出数据
    /// \return 评估状态
    USI Eval_All(const USI& j, const OCP_DBL& val, vector<OCP_DBL>& outdata) const;

    /// 在表格中插值，评估所有列，j = 0，返回数据索引从 1 开始
    /// \param val 插值的值
    /// \param outdata 输出数据
    /// \return 评估状态
    USI Eval_All0(const OCP_DBL& val, vector<OCP_DBL>& outdata) const;

    /// 在表格中插值，评估所有列并返回斜率，j = 0，返回数据索引从 1 开始
    /// \param val 插值的值
    /// \param outdata 输出数据
    /// \param slope 斜率
    /// \return 评估状态
    USI Eval_All0(const OCP_DBL& val, vector<OCP_DBL>& outdata, vector<OCP_DBL>& slope) const;

    /// 在表格中插值，评估目标列
    /// \param j 指定列索引
    /// \param val 插值的值
    /// \param destj 目标列索引
    /// \return 插值结果
    OCP_DBL Eval(const USI& j, const OCP_DBL& val, const USI& destj) const;

    /// 在表格中插值，评估目标列，并返回对应斜率
    /// \param j 指定列索引
    /// \param val 插值的值
    /// \param destj 目标列索引
    /// \param myK 斜率
    /// \return 插值结果
    OCP_DBL Eval(const USI& j, const OCP_DBL& val, const USI& destj, OCP_DBL& myK) const;

    /// 在表格中插值，评估单调递减的指定列
    /// \param j 指定列索引
    /// \param val 插值的值
    /// \param destj 目标列索引
    /// \return 插值结果
    OCP_DBL Eval_Inv(const USI& j, const OCP_DBL& val, const USI& destj) const;

    /// 返回距离特定值最近的行，给定列 j
    /// \param j 列索引
    /// \param val 特定值
    /// \param outdata 输出数据
    void GetCloseRow(const USI& j, const OCP_DBL& val, vector<OCP_DBL>& outdata) const;

    /// 在屏幕上显示表格数据
    void Display() const;

protected:
    USI                     nRow; ///< 表格的行数
    USI                     nCol; ///< 表格的列数
    mutable USI             bId;  ///< 插值时的起始行点
    vector<vector<OCP_DBL>> data; ///< 表格数据，data[i] 是第 i 列
};

/// 2D 表格类
class OCPTable2{
public:
    /// 设置表格
    /// \param tab 2D 表格数据
    void Setup(const Table2& tab);

    /// 在表格中插值，评估所有列
    /// \param val1 第一个值
    /// \param val2 第二个值
    /// \param j2 第二个列索引
    /// \param out 输出数据
    void Eval_All(const OCP_DBL& val1, const OCP_DBL& val2, const USI& j2, vector<OCP_DBL>& out) const;

    /// 在表格中插值，评估所有列，并返回两个斜率
    /// \param val1 第一个值
    /// \param val2 第二个值
    /// \param j2 第二个列索引
    /// \param out 输出数据
    /// \param slope1 第一个斜率
    /// \param slope2 第二个斜率
    void Eval_All(const OCP_DBL& val1, const OCP_DBL& val2, const USI& j2, vector<OCP_DBL>& out,
        vector<OCP_DBL>& slope1, vector<OCP_DBL>& slope2) const;

    /// 在表格中插值，评估所有列，ref 和 tables[i][0] 都是升序
    /// \param val1 第一个值
    /// \param val2 第二个值
    /// \param out 输出数据
    void Eval_All0(const OCP_DBL& val1, const OCP_DBL& val2, vector<OCP_DBL>& out) const;

    /// 在表格中插值，评估所有列，并返回两个斜率，ref 和 tables[i][0] 都是升序
    /// \param val1 第一个值
    /// \param val2 第二个值
    /// \param out 输出数据
    /// \param slope1 第一个斜率
    /// \param slope2 第二个斜率
    void Eval_All0(const OCP_DBL& val1, const OCP_DBL& val2, vector<OCP_DBL>& out,
        vector<OCP_DBL>& slope1, vector<OCP_DBL>& slope2) const;

    /// 在表格中插值，评估目标列
    /// \param val1 第一个值
    /// \param val2 第二个值
    /// \param j2 第二个列索引
    /// \param destj 目标列索引
    /// \return 插值结果
    OCP_DBL Eval(const OCP_DBL& val1, const OCP_DBL& val2, const USI& j2, const USI& destj) const;

    /// 在表格中插值，评估目标列，并返回对应斜率
    /// \param val1 第一个值
    /// \param val2 第二个值
    /// \param j2 第二个列索引
    /// \param destj 目标列索引
    /// \param myK1 第一个斜率
    /// \param myK2 第二个斜率
    /// \return 插值结果
    OCP_DBL Eval(const OCP_DBL& val1, const OCP_DBL& val2, const USI& j2, const USI& destj, OCP_DBL& myK1, OCP_DBL& myK2) const;

    /// 获取 1D 表格的列数
    /// \return 列数
    auto GetColNum() const { return nCol; }

    /// 判断表格是否为空
    /// \return 是否为空的布尔值
    auto IsEmpty() const { return numtable == 0; }

protected:
    USI                       numtable; ///< 表格数量
    vector<OCP_DBL>           ref;      ///< 参考值
    vector<OCPTable>          tables;   ///< 表格数据
    USI                       nCol;     ///< 列数
    mutable vector<OCP_DBL>   data1;    ///< 临时数据1
    mutable vector<OCP_DBL>   data2;    ///< 临时数据2
    mutable vector<OCP_DBL>   cdata1;   ///< 计算数据1
    mutable vector<OCP_DBL>   cdata2;   ///< 计算数据2
};

#endif /* end if __OCPTable_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/