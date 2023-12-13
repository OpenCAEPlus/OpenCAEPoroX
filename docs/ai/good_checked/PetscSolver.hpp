/*! 
 * \file    PetscSolver.hpp 
 * \brief   Petsc Solver的API
 * \author  Shizhe Li 
 * \date    2023年5月4日
 * \note    使用Petsc作为线性求解器
 * \details
 * 使用Petsc作为线性求解器。这个文件主要包括了PetscSolver类，ScalarPetscSolver类和VectorPetscSolver类。
 * 这些类主要用于解决线性问题。
 * 
 * -----------------------------------------------------------------------------------
 * 版权所有 (C) 2021--至今 OpenCAEPoroX团队。保留所有权利。
 * 根据GNU Lesser General Public License 3.0或更高版本的条款发布。
 * -----------------------------------------------------------------------------------
 */

#ifdef WITH_PETSCSOLVER
#ifndef __PETSCSOLVER_HEADER__
#define __PETSCSOLVER_HEADER__

#include "LinearSolver.hpp"
#include "PETScBSolverPS.h"
#include <vector>

using namespace std;

/*!
 * \class PetscSolver
 * \brief Petsc求解器的API
 * \note  索引是从0开始的
 */
class PetscSolver : public LinearSolver
{
public:
    PetscSolver() = default;
    
    /*!
     * \brief 设置参数
     * \param dir 目录
     * \param file 文件
     */
    void SetupParam(const string& dir, const string& file) override {};
    
    /// 初始化线性求解器的参数
    void InitParam() override {};
    
    /// 为pardiso求解器分配内存
    void Allocate(const OCPMatrix& mat) override;
    
    /// 计算通信中使用的术语
    void CalCommTerm(const USI& actWellNum, const Domain* domain) override;
    
    /// 组装系数矩阵
    void AssembleMat(OCPMatrix& mat) override;
    
    /// 获取迭代求解器使用的迭代次数
    USI GetNumIters() const override { return 1; }

protected:
    // CSR/BSR
    OCP_INT                 nb;         ///< 块维度
    vector<OCP_DBL>         A;          ///< 值
    vector<OCP_SLL>         iA;         ///< 行指针
    vector<OCP_SLL>         jA;         ///< 列索引
    OCP_DBL*                b;          ///< 右侧项
    OCP_DBL*                x;          ///< 解
    // 通信
    MPI_Comm                myComm;       ///< 通信器
    OCP_INT                 numproc;      ///< 进程数
    OCP_INT                 myrank;       ///< 当前等级
    const vector<OCP_ULL>*  global_index; ///< 全局索引
    vector<OCP_SLL>         allBegin;     ///< 所有进程（包括自身）的开始
    vector<OCP_SLL>         allEnd;       ///< 所有进程（包括自身）的结束
    vector<OCP_INT>         allEle;       ///< 每个进程的元素数量
};

/*!
 * \class ScalarPetscSolver 
 * \brief 继承自PetscSolver的类，用于解决标量问题
 */
class ScalarPetscSolver : public PetscSolver
{
public:
    ScalarPetscSolver() {};

    /// 解决线性系统
    OCP_INT Solve() override;
};

/*!
 * \class VectorPetscSolver 
 * \brief 继承自PetscSolver的类，用于解决向量问题
 */
class VectorPetscSolver : public PetscSolver
{
public:
    VectorPetscSolver(const Domain* domain) {
        myComm  = domain->myComm;
        numproc = domain->numproc;
        myrank  = domain->myrank;
        allBegin.resize(numproc);
        allEnd.resize(numproc);
        allEle.resize(numproc);
    }

    /// 解决线性系统
    OCP_INT Solve() override;
};

#endif 
#endif // WITH_PETSCSOLVER

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           May/04/2023      Create file                          */
/*----------------------------------------------------------------------------*/