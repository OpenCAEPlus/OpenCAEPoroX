/*! \file    FaspSolver.hpp
 *  \brief   FASP求解器接口类声明
 *  \details 本文件包含了用于FASP求解器的类声明，这些类是与FASP库接口的封装。
 *  \author  Shizhe Li
 *  \date    Nov/22/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#if WITH_FASP

#ifndef __FASPSOLVER_HEADER__
#define __FASPSOLVER_HEADER__

// 标准头文件
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// OpenCAEPoroX头文件
#include "LinearSolver.hpp"

// faspsolver头文件
extern "C" {
#if OCPFLOATTYPEWIDTH == 64
#define SHORT  short
#define INT    int
#define LONG   long
#define REAL   double
#endif // OCPFLOATTYPEWIDTH == 64

#if OCPFLOATTYPEWIDTH == 128
#define SHORT  short
#define INT    int
#define LONG   long
#define REAL   long double
#endif // OCPFLOATTYPEWIDTH == 128

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

#undef SHORT
#undef INT
#undef LONG
#undef REAL
}

// faspcpr头文件
#if WITH_FASPCPR
extern "C" {
#include "faspcpr.h"
#include "faspcpr_functs.h"
}
#endif

// fasp4blkoil头文件
#if WITH_FASP4BLKOIL
extern "C" {
#if OCPFLOATTYPEWIDTH == 64
#define SHORT  short
#define INT    int
#define LONG   long
#define REAL   double
#endif // OCPFLOATTYPEWIDTH == 64

#if OCPFLOATTYPEWIDTH == 128
#define SHORT  short
#define INT    int
#define LONG   long
#define REAL   long double
#endif // OCPFLOATTYPEWIDTH == 128

#include "fasp4blkoil.h"
#include "fasp4blkoil_functs.h"

#undef SHORT
#undef INT
#undef LONG
#undef REAL
}
#endif

// fasp4cuda头文件
#if WITH_FASP4CUDA
#include "fasp4cuda.h"
#include "fasp4cuda_functs.h"
#endif

using namespace std;

// 标准预处理器类型
#define PC_NULL  60 ///< 无预处理器
#define PC_FASP1 61 ///< FASP1: MSP, 2020年FIM默认
#define PC_FASP2 62 ///< FASP2: MSP, 仅实验性使用
#define PC_FASP3 63 ///< FASP3: MSP, 整体预处理器
#define PC_FASP4 64 ///< FASP4: MSP, 2015年FIM默认
#define PC_FASP5 65 ///< FASP5: MSP, 仅实验性使用
#define PC_DIAG  68 ///< DIAG: 对角线预处理器
#define PC_BILU  69 ///< BILU: 块ILU预处理器

// 共享设置预处理器类型
#define PC_FASP1_SHARE 71 ///< FASP1共享设置阶段, 小心使用
#define PC_FASP4_SHARE 74 ///< FASP4共享设置阶段, 小心使用
#define RESET_CONST    35 ///< FASP1_SHARE与FASP4_SHARE的共享阈值

/// FASP求解器基类
class FaspSolver : public LinearSolver
{
public:
    /// 设置FASP参数
    void SetupParam(const string& dir, const string& file) override;

    /// 计算通信中使用的项
    void CalCommTerm(const USI& actWellNum, const Domain* domain) override {}

    /// 获取迭代求解器使用的迭代次数
    USI GetNumIters() const override { return itsParam.maxit; }

public:
    string      solveDir;  ///< 当前工作目录
    string      solveFile; ///< fasp文件的相对路径
    input_param inParam;   ///< 输入文件中的参数
    ITS_param   itsParam;  ///< 迭代方法的参数
    AMG_param   amgParam;  ///< AMG方法的参数
    ILU_param   iluParam;  ///< ILU方法的参数
    SWZ_param   swzParam;  ///< Schwarz方法的参数
};

/// FASP的CSR格式标量求解器
class ScalarFaspSolver : public FaspSolver
{
public:
    ScalarFaspSolver() {};

protected:
    /// 为线性系统分配内存
    void Allocate(const OCPMatrix& mat) override;

    /// 初始化线性求解器的参数
    void InitParam() override;

    /// 组装系数矩阵
    void AssembleMat(OCPMatrix& mat) override;

    /// 解线性系统
    OCP_INT Solve() override;

protected:
    dCSRmat A; ///< 标量问题的矩阵
    dvector b; ///< 标量问题的右侧项
    dvector x; ///< 标量问题的解
};

/// FASP的BSR格式向量求解器
class VectorFaspSolver : public FaspSolver
{
public:
    VectorFaspSolver() {};

protected:
    /// 为线性系统分配内存
    void Allocate(const OCPMatrix& mat) override;

    /// 初始化线性求解器的参数
    void InitParam() override;

    /// 组装系数矩阵
    void AssembleMat(OCPMatrix& mat) override;

    /// 解线性系统
    OCP_INT Solve() override;

    /// 对线性系统进行解耦
    void Decoupling(dBSRmat* Absr,
                    dvector* b,
                    dBSRmat* Asc,
                    dvector* fsc,
                    ivector* order,
                    OCP_DBL*  Dmatvec,
                    int      decouple_type);

protected:
    dBSRmat A; ///< 向量问题的矩阵
    dvector b; ///< 向量问题的右侧项
    dvector x; ///< 向量问题的解
    dBSRmat Asc;   ///< 向量问题的缩放矩阵
    dvector fsc;   ///< 向量问题的缩放右侧项
    ivector order; ///< 平滑过程中用户定义的排序
    vector<OCP_DBL> Dmat; ///< 解耦矩阵
};

#endif // __FASPSOLVER_HEADER__

#endif // WITH_FASP


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/22/2021      Create file                          */
/*  Chensong Zhang      Jan/08/2022      Update Doxygen                       */
/*  Chensong Zhang      Jan/19/2022      Set FASP4BLKOIL as optional          */
/*  Li Zhao             Apr/04/2022      Set FASP4CUDA   as optional          */
/*----------------------------------------------------------------------------*/
