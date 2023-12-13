/*! \file    OCPMatrix.hpp
 *  \brief   OCPMatrix类声明
 *  \author  Shizhe Li
 *  \date    Dec/06/2023
 *
 *-----------------------------------------------------------------------------------
 *  版权所有 (C) 2021--至今 OpenCAEPoroX团队。保留所有权利。
 *  根据GNU Lesser General Public License 3.0或更高版本的条款发布。
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPMATRIX_HEADER__
#define __OCPMATRIX_HEADER__

// 标准头文件
#include <vector>

// OpenCAEPoroX头文件
#include "OCPConst.hpp"
#include "Domain.hpp"

using namespace std;

/// OCP中的内部矩阵格式
class OCPMatrix{
public:
    /// 分配内存
    void Allocate(const Domain& domain, const USI& blockdim);
    /// 清除数据但不释放内存
    void ClearData();
    /// 设置块维度和块大小
    void SetBlockDim(const USI& blockdim);

protected:
    /// 为行分配内存
    void AllocateRowMem(const Domain& domain);
    /// 为列分配内存
    void AllocateColMem(const Domain& domain);

public:
    /// 设置维度
    OCP_USI AddDim(const OCP_USI& n);
    /// 第一次在对角线上添加新值
    void NewDiag(const OCP_USI& n, const OCP_DBL& v);
    void NewDiag(const OCP_USI& n, const vector<OCP_DBL>& v);
    /// 在对角线上添加值
    void AddDiag(const OCP_USI& n, const OCP_DBL& v);
    void AddDiag(const OCP_USI& n, const vector<OCP_DBL>& v);
    /// 在非对角线上添加值
    void NewOffDiag(const OCP_USI& bId, const OCP_USI& eId, const OCP_DBL& v);
    void NewOffDiag(const OCP_USI& bId, const OCP_USI& eId, const vector<OCP_DBL>& v);
    /// 在右手边添加值
    void AddRhs(const OCP_USI& n, const OCP_DBL& v);
    void AddRhs(const OCP_USI& n, const vector<OCP_DBL>& v);
    /// 将右手边复制到b
    void CopyRhs(const vector<OCP_DBL>& rhs);
    /// 设置一个猜测值
    void SetGuess(const OCP_USI& n, const OCP_DBL& v) { u[n] = v; }
    /// 返回解
    auto& GetSolution() { return u; }

public:
    /// 将A和b输出到文件
    void OutputLinearSystem(const string& dir, const string& fileA, const string& fileb) const;
    /// 将解输出到文件
    void OutputSolution(const string& dir, const string& fileU) const;

public:
    /// 检查线性系统中是否存在NaN值
    void CheckLinearSystem() const;
    /// 检查解中是否存在NaN值
    void CheckSolution() const;

public:
    /// 小块矩阵的维度
    USI nb;
    /// 小块矩阵的大小
    USI nb2;
    /// 最大的nb，用于矩阵的重用
    USI maxNb{ 0 };
    /// 矩阵的最大可能维度
    OCP_USI maxDim;
    /// 矩阵的实际维度
    OCP_USI dim;
    /// 最大的非零元素个数
    OCP_USI max_nnz;
    /// 非零元素的列索引
    vector<vector<OCP_USI>> colId;
    /// 非零元素的值
    vector<vector<OCP_DBL>> val;
    /// 线性系统的右手边
    vector<OCP_DBL> b;
    /// 线性系统的解
    vector<OCP_DBL> u;
};

#endif /* end if __OCPMATRIX_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/06/2023      Create file                          */
/*----------------------------------------------------------------------------*/