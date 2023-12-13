/*! \file    SamgSolver.hpp 
 *  \brief   SANG Solver的API 
 *  \author  Shizhe Li 
 *  \date    Apr/07/2023 
 * 
 *  \note    使用SAMG作为线性求解器 
 * 
 *----------------------------------------------------------------------------------- 
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved. 
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later. 
 *----------------------------------------------------------------------------------- 
 */
#ifdef WITH_SAMG
#ifndef __SAMGSOLVER_HEADER__
#define __SAMGSOLVER_HEADER__

#include "LinearSolver.hpp"
#include "samg.h"
using namespace std;

/// SAMG Solver的API
// 注意：索引从1开始
class SamgSolver : public LinearSolver{
public:
    /// 设置参数
    void SetupParam(const string& dir, const string& file) override;
    /// 初始化线性求解器参数
    void InitParam() override;
    /// 为Pardiso求解器分配内存
    void Allocate(const OCPMatrix& mat) override;
    /// 计算通信中使用的项
    void CalCommTerm(const USI& actWellNum, const Domain* domain) override;
    /// 解线性系统
    OCP_INT Solve() override;
    /// 获取迭代求解器使用的迭代次数
    USI GetNumIters() const override { return 1; }
protected:
    // CSR矩阵
    SAMG_INT            nb;            ///< 矩阵的块数
    vector<SAMG_INT>    iA;            ///< CSR格式的行指针
    vector<SAMG_INT>    jA;            ///< CSR格式的列指针
    vector<SAMG_REAL>   A;             ///< CSR格式的非零元素
    SAMG_REAL*          b = nullptr;   ///< 右端项
    SAMG_REAL*          x = nullptr;   ///< 解向量

    // 物理信息
    SAMG_INT           nsys;            ///< 未知数的数量（物理变量类型）
    SAMG_INT           nnu;             ///< （本地）变量的数量
    SAMG_INT           nna;             ///< 存储在（本地）向量a（本地nnz）中的矩阵条目的数量
    SAMG_INT           ndiu;            ///< iu的大小
    SAMG_INT           ndip;            ///< ip的大小
    vector<SAMG_INT>   iu;              ///< 变量到未知数指针：nnu
    vector<SAMG_INT>   ip;              ///< 变量到点指针：nnu
    vector<SAMG_INT>   iu_tmp;          ///< iu的模板：nb
    vector<SAMG_INT>   nunknown_description;        // 通信
    SAMG_INT            myComm = (SAMG_INT)MPI_Comm_c2f(MPI_COMM_WORLD);
    SAMG_INT            npsnd;          ///< 发送的邻居处理器的总数
    vector<SAMG_INT>    iranksnd;       ///< 发送的邻居处理器的排名
    SAMG_INT            nshalo;         ///< 发送的变量的总数
    vector<SAMG_INT>    ipts;           ///< 在isndlist中发送的变量的范围
    vector<SAMG_INT>    isndlist;       ///< 发送的变量的索引
    SAMG_INT            nprec;          ///< 接收的邻居处理器的总数
    vector<SAMG_INT>    irankrec;       ///< 接收的邻居处理器的排名
    SAMG_INT            nrhalo;         ///< 接收的变量的总数
    vector<SAMG_INT>    iptr;           ///< 在isndlist中接收的变量的范围
    vector<SAMG_INT>    ireclist;       ///< 接收的变量的索引
    const vector<OCP_USI>* global_index;

protected:
    // SAMG参数
    SAMG_INT          noil_approach            = 12;    ///< SAMG的正交插值方法
    SAMG_INT          noil_cyc                 = 19;    ///< SAMG的正交插值方法
    SAMG_REAL         noil_preparation         = 19.4;  ///< SAMG的正交插值方法
    SAMG_INT          ierr                     = 0;     ///< 表示错误或警告的代码编号
    SAMG_INT          samg_matrix              = 220;   ///< 矩阵A的类型
    SAMG_INT          ncyc_done                = 0;     ///< 执行的总循环（迭代）次数
    SAMG_REAL         res_in                   = 0;     ///< 初始猜测的残差
    SAMG_REAL         res_out                  = 0;     ///< 最终近似值的残差
    SAMG_INT          nsolve                   = 2;     ///< 指定SAMG的求解策略
    SAMG_INT          ncyc                     = 13050; ///< 循环和加速策略
    SAMG_INT          ifirst                   = 1;     ///< 选择初始猜测的方式
    SAMG_INT          iswtch                   = 51;    ///< 内存扩展开关
    SAMG_INT          iout                     = -1;    ///< 在解决阶段打印输出
    SAMG_INT          idump                    = -1;    ///< 在设置阶段打印输出
    SAMG_REAL         eps                      = 1.0e-3;///< AMG迭代的标准停止准则
    SAMG_REAL         chktol                   = -1.0;  ///< 输入矩阵的检查
    SAMG_REAL         a_cmplx                  = 2.5;   ///< 用于分配SAMG的初始内存
    SAMG_REAL         g_cmplx                  = 1.8;   ///< 用于分配SAMG的初始内存
    SAMG_REAL         p_cmplx                  = 2.0;   ///< 用于分配SAMG的初始内存
    SAMG_REAL         w_avrge                  = 2.5;   ///< 用于分配SAMG的初始内存
};

// 将内部的矩阵（类似CSR）转换为CSR矩阵
class ScalarSamgSolver : public SamgSolver{
public:
    ScalarSamgSolver(const OCPModel& model) {
        nb       = 1;
        nsys     = 1;
        ndiu     = 1;
        ndip     = 1;
        ifirst   = 0;   // 使用上一个解作为初始猜测
    }
    /// 组装系数矩阵
    void AssembleMat(OCPMatrix& mat) override;
};

// 将内部的矩阵（类似BSR）转换为CSR矩阵
class VectorSamgSolver : public SamgSolver{
public:
    VectorSamgSolver(const USI& blockDim, const OCPModel& model) {
        nb = blockDim;
        nsys     = nb;
        iu_tmp.resize(nb);
        ifirst   = 1;   // 使用零解作为初始猜测
        for (USI i = 0; i < nb; i++)  iu_tmp[i] = i + 1;
        if (model == OCPModel::isothermal) {
            nunknown_description.resize(nsys, 2);
            nunknown_description[0] = 0;          // 压力
        }
        else if (model == OCPModel::thermal) {
            nunknown_description.resize(nsys, 2); // 浓度（初始）
            nunknown_description.front() = 0;     // 压力
            nunknown_description.back()  = 100;   // 温度
        }
        else
            OCP_ABORT("Wrong Model for SAMG Solver!");
    }
    /// 组装系数矩阵
    void AssembleMat(OCPMatrix& mat) override;
};

#endif
#endif // WITH_SAMG

 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Apr/07/2023      Create file                          */
 /*----------------------------------------------------------------------------*/