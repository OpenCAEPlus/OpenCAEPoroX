/*! \file    LinearSystem.hpp
 *  \brief   线性系统求解类的声明
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *  本文件包含了线性系统求解类的声明，用于离散系统的线性求解。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __LINEARSYSTEM_HEADER__
#define __LINEARSYSTEM_HEADER__

// 标准头文件
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// OpenCAEPoroX 头文件
#include "OCPConst.hpp"
#include "OCPControlMethod.hpp"
#include "OCPMatrix.hpp"
#include "PardisoSolver.hpp"
#include "SamgSolver.hpp"
#include "FaspSolver.hpp"
#include "PetscSolver.hpp"
#include "LinearSolver.hpp"

using namespace std;

/**
 * \class LinearSystem
 * \brief 离散系统的线性求解器。
 *
 * 注意：矩阵内部以行分割的CSR格式存储。
 */
class LinearSystem
{
public:
    /**
     * \brief 设置线性系统
     * \param model 模型对象
     * \param dir 工作目录
     * \param file 文件名
     * \param d 领域
     * \param nb 块维度
     * \return 成功设置的维度
     */
    USI Setup(const OCPModel& model, const string& dir, const string& file, const Domain& d, const USI& nb);

    /**
     * \brief 设置工作线性求解器
     * \param i 线性求解器的索引
     */
    void SetWorkLS(const USI& i);

    /**
     * \brief 清除标量问题的内部矩阵数据。
     */
    void ClearData() { mat.ClearData(); }

    /**
     * \brief 计算全局开始
     * \param actWellNum 活跃井数量
     */
    void CalCommTerm(const USI& actWellNum) { LS[wIndex]->CalCommTerm(actWellNum, domain); }

    /**
     * \brief 为线性求解器组装矩阵。
     */
    void AssembleMatLinearSolver();

    /**
     * \brief 解线性系统。
     * \return 解算结果的状态码
     */
    OCP_INT Solve();

protected:
    /**
     * \brief 设置线性求解器。
     * \param model 模型对象
     * \param lsFile 线性求解器配置文件
     */
    void SetupLinearSolver(const OCPModel& model, const string& lsFile);

public:
    /**
     * \brief 设置维度。
     * \param n 新增的维度
     * \return 新增维度的数量
     */
    OCP_USI AddDim(const OCP_USI& n) { return mat.AddDim(n); }

    /**
     * \brief 在首位位置添加一个对角线值。
     * \param n 维度
     * \param v 对角线值
     */
    void NewDiag(const OCP_USI& n, const OCP_DBL& v) { mat.NewDiag(n, v); }
    void NewDiag(const OCP_USI& n, const vector<OCP_DBL>& v) { mat.NewDiag(n, v); }

    /**
     * \brief 在对角线值上添加一个值。
     * \param n 维度
     * \param v 要添加的值
     */
    void AddDiag(const OCP_USI& n, const OCP_DBL& v) { mat.AddDiag(n, v); }
    void AddDiag(const OCP_USI& n, const vector<OCP_DBL>& v) { mat.AddDiag(n, v); }

    /**
     * \brief 添加一个非对角线值。
     * \param bId 开始索引
     * \param eId 结束索引
     * \param v 值
     */
    void NewOffDiag(const OCP_USI& bId, const OCP_USI& eId, const OCP_DBL& v) { mat.NewOffDiag(bId, eId, v); }
    void NewOffDiag(const OCP_USI& bId, const OCP_USI& eId, const vector<OCP_DBL>& v) { mat.NewOffDiag(bId, eId, v); }

    /**
     * \brief 在b[n]上添加一个值。
     * \param n 维度
     * \param v 值
     */
    void AddRhs(const OCP_USI& n, const OCP_DBL& v) { mat.AddRhs(n, v); }
    void AddRhs(const OCP_USI& n, const vector<OCP_DBL>& v) { mat.AddRhs(n, v); }

    /**
     * \brief 将rhs复制到b。
     * \param rhs 右手边向量
     */
    void CopyRhs(const vector<OCP_DBL>& rhs) { mat.CopyRhs(rhs); }

    /**
     * \brief 在u[n]上设置一个初始值。
     * \param n 维度
     * \param v 初始值
     */
    void SetGuess(const OCP_USI& n, const OCP_DBL& v) { mat.SetGuess(n, v); }

    /**
     * \brief 返回解。
     * \return 解向量
     */
    vector<OCP_DBL>& GetSolution() { return mat.GetSolution(); }

public:
    /**
     * \brief 将矩阵和rhs输出到文件fileA和fileb。
     * \param fileA 矩阵文件
     * \param fileb 右手边向量文件
     */
    void OutputLinearSystem(const string& fileA, const string& fileb) const { mat.OutputLinearSystem(solveDir, fileA, fileb); }

    /**
     * \brief 将解输出到文件。
     * \param fileX 解向量文件
     */
    void OutputSolution(const string& fileX) const { mat.OutputSolution(solveDir, fileX); }

protected:
    /// 工作索引
    USI                   wIndex{ 0 };
    /// 当前工作目录
    string                solveDir;
    /// 领域分解信息
    const Domain*         domain;
    /// 内部矩阵
    OCPMatrix             mat;
    /// 块维度集
    vector<USI>           nbs;
    /// 线性求解器类型集
    vector<OCPLStype>     LStype;
    /// 线性求解器集
    vector<LinearSolver*> LS;
};

#endif /* end if __LINEARSOLVER_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      创建文件                              */
/*  Chensong Zhang      Oct/15/2021      格式化文件                            */
/*  Chensong Zhang      Nov/09/2021      移除解耦方法                          */
/*  Chensong Zhang      Nov/22/2021      文件重命名为LinearSystem              */
/*----------------------------------------------------------------------------*/
