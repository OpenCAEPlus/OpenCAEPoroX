/*! \file    PardisoSolver.hpp
 *  \brief   PardisoSolver类的声明文件
 *  \author  Shizhe Li
 *  \date    Mar/30/2023
 *
 *  \note    OpenCAEPoroX中使用的参数与SLB的Eclipse大体兼容，
 *           但它有自己的一些规则以便于使用。它是可扩展的和友好的。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifdef WITH_PARDISO
#ifndef __PARDISOSOLVER_HEADER__
#define __PARDISOSOLVER_HEADER__

#include <math.h>
#include <mpi.h>
#include <mkl.h>
#include <mkl_cluster_sparse_solver.h>
#include <algorithm>
#include "LinearSolver.hpp"

#if OCPFLOATTYPEWIDTH == 64
#ifdef MKL_ILP64
#define MPI_DT MPI_LONG
#else
#define MPI_DT MPI_INT
#endif

#define MPI_REDUCE_AND_BCAST \
        MPI_Reduce(&err_mem, &error, 1, MPI_DT, MPI_SUM, 0, MPI_COMM_WORLD); \
        MPI_Bcast(&error, 1, MPI_DT, 0, MPI_COMM_WORLD);

using namespace std;

/**
 * \class PardisoSolver
 * \brief 一个统一的API，用于CSR和BSR矩阵（BSR目前不工作）
 *
 * PardisoSolver是一个线性求解器，用于解决稀疏矩阵的线性系统。
 * 它继承自LinearSolver类。
 */
class PardisoSolver : public LinearSolver
{
public:
    PardisoSolver() {};

    /// 设置参数。
    void SetupParam(const string& dir, const string& file) override;

    /// 初始化线性求解器的参数。
    void InitParam() override;

    /// 为pardiso求解器分配内存
    void Allocate(const OCPMatrix& mat) override;

    /// 计算通信所需的项
    void CalCommTerm(const USI& actWellNum, const Domain* domain) override;

    /// 组装系数矩阵。
    void AssembleMat(OCPMatrix& mat) override;

    /// 解决线性系统。
    OCP_INT Solve() override;

    /// 获取迭代求解器使用的迭代次数。
    USI GetNumIters() const override { return 1; }

protected:
    void*     pt[64]    = { 0 };  ///< 内部求解器内存指针pt, 不要修改它。
    MKL_INT   iparm[64] = { 0 };  ///< 集群稀疏求解器控制参数。
    MKL_INT   nrhs      = 1;      ///< 右手边的数量。
    MKL_INT   maxfct    = 1;      ///< 数值分解的最大数量。
    MKL_INT   mnum      = 1;      ///< 使用哪个分解。
    MKL_INT   mtype     = 11;     ///< 定义矩阵类型，影响选轴方法。
    MKL_INT   phase;              ///< 控制求解器的执行。
    MKL_INT   N;                  ///< 矩阵的维度。
    MKL_INT   msglvl    = 0;      ///< 消息级别信息。
    MKL_INT   error     = 0;      ///< 错误信息。
    MKL_INT   err_mem   = 0;      ///< 内存错误信息。
    double    ddum;               ///< 双精度虚拟变量。
    MKL_INT   idum;               ///< 整型虚拟变量。

    // CSR/BSR矩阵
    MKL_INT          nb;          ///< BSR矩阵的块大小。
    vector<MKL_INT>  iA;          ///< 行索引数组。
    vector<MKL_INT>  jA;          ///< 列索引数组。
    vector<double>   A;           ///< 非零元素数组。
    double*          b = nullptr; ///< 右手边向量。
    double*          x = nullptr; ///< 解向量。

    // 通信
    int              myComm = MPI_Comm_c2f(MPI_COMM_WORLD); ///< MPI通信域。
    const vector<OCP_ULL>* global_index; ///< 全局索引。
};

/**
 * \class VectorPardisoSolver
 * \brief 一个特殊的PardisoSolver，用于处理向量形式的数据。
 *
 * VectorPardisoSolver继承自PardisoSolver，并提供了特定于向量的方法。
 */
class VectorPardisoSolver : public PardisoSolver
{
public:
    VectorPardisoSolver() {};

    /// 为pardiso求解器分配内存。
    void Allocate(const OCPMatrix& mat) override;

    /// 计算通信所需的项。
    void CalCommTerm(const USI& actWellNum, const Domain* domain) override;

    /// 组装系数矩阵。
    void AssembleMat(OCPMatrix& mat) override;
};

#endif // OCPFLOATTYPEWIDTH == 64
#endif // __PARDISOSOLVER_HEADER__
#endif // WITH_PARDISO

 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Mar/30/2023      Create file                          */
 /*----------------------------------------------------------------------------*/