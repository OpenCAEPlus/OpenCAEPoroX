/*! \file    LinearSolver.hpp
*  \brief   LinearSolver类声明
*  \author  Shizhe Li
*  \date    Nov/22/2021
*
*-----------------------------------------------------------------------------------
*  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
*  Released under the terms of the GNU Lesser General Public License 3.0 or later.
*-----------------------------------------------------------------------------------
*/

#ifndef __LINEARSOLVER_HEADER__
#define __LINEARSOLVER_HEADER__

// 标准头文件
#include <string>
#include <vector>

// OpenCAEPoroX头文件
#include "OCPConst.hpp"
#include "Domain.hpp"
#include "OCPMatrix.hpp"

using namespace std;

// OCP线性求解器类型
enum class OCPLStype : USI{
    none,    ///< 无
    fasp,    ///< 使用FASP方法
    pardiso, ///< 使用Pardiso方法
    petsc,   ///< 使用Petsc方法
    samg     ///< 使用SAMG方法
};

/// 线性求解器的虚基类
class LinearSolver{
public:
    /// 从输入文件中读取线性求解器的参数
    virtual void SetupParam(const string& dir, const string& file) = 0;

    /// 为线性求解器分配最大内存
    virtual void Allocate(const OCPMatrix& mat) = 0;

    /// 计算通信中使用的项
    virtual void CalCommTerm(const USI& actWellNum, const Domain* domain) = 0;

    /// 从内部矩阵数据组装线性求解器的矩阵
    virtual void AssembleMat(OCPMatrix& mat) = 0;

    /// 解决线性系统并返回迭代次数
    virtual OCP_INT Solve() = 0;

    /// 获取迭代次数
    virtual USI GetNumIters() const = 0;

protected:
    /// 初始化线性求解器的参数
    virtual void InitParam() = 0;
};


#endif // __LINEARSOLVER_HEADER__

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/22/2021      Create file                          */
/*  Chensong Zhang      Jan/18/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/