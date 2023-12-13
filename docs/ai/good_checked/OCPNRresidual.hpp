/*! 
 *  \file    OCPNRresidual.hpp
 *  \brief   本文件包含了在非线性迭代中使用的数据结构
 *  \author  Shizhe Li
 *  \date    Oct/10/2023
 *
 *  本文件定义了一个用于存储和管理非线性迭代中残差信息的类。
 *  它包含了各种残差的度量和相关的最大值信息。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPNRRESIDUAL_HEADER__
#define __OCPNRRESIDUAL_HEADER__

// 标准头文件
#include <vector>

// OpenCAEPoroX头文件
#include "OCPConst.hpp"

using namespace std;

/*!
 *  \class   OCPNRresidual
 *  \brief   非线性迭代残差数据结构类
 *
 *  该类包含了在非线性迭代过程中计算的各种残差值，以及相关的最大值信息。
 */
class OCPNRresidual
{
public:
    /*!
     *  \brief  初始化等温残差存储结构
     *  \param  nb 模型中的块(block)数量
     *  \param  nw 模型中的水相(water phase)数量
     *  \param  nc 模型中的组分(component)数量
     */
    void SetupIsoT(const OCP_USI& nb, const OCP_USI& nw, const OCP_USI& nc);
    
    /*!
     *  \brief  初始化温度相关残差存储结构
     *  \param  nb 模型中的块(block)数量
     *  \param  nw 模型中的水相(water phase)数量
     *  \param  nc 模型中的组分(component)数量
     */
    void SetupT(const OCP_USI& nb, const OCP_USI& nw, const OCP_USI& nc);
    
    /*!
     *  \brief  将所有残差值设置为零
     */
    void SetZero();
    
    /// \brief 每个块(block)的所有方程的绝对残差
    vector<OCP_DBL> resAbs;
    
    /// \brief 相对于孔隙体积的所有方程的2-范数相对残差
    vector<OCP_DBL> resRelV;
    
    /// \brief 相对于总摩尔数的质量守恒方程的2-范数相对残差
    vector<OCP_DBL> resRelN;
    
    /// \brief 相对于总能量的能量守恒方程的2-范数相对残差
    vector<OCP_DBL> resRelE;
    
    /// \brief 每个块(block)的初始最大相对残差（相对于孔隙体积）
    OCP_DBL maxRelRes0_V;
    
    /// \brief 每个块(block)的迭代中的最大相对残差（相对于孔隙体积）
    OCP_DBL maxRelRes_V;
    
    /// \brief 每个块(block)的最大相对残差（相对于总摩尔数）
    OCP_DBL maxRelRes_N;
    
    /// \brief 每个块(block)的最大相对残差（相对于总能量）
    OCP_DBL maxRelRes_E;
    
    /// \brief 每口井的最大相对残差（相对于总摩尔数）
    OCP_DBL maxWellRelRes_mol;
    
    // 使用负数表示井号（待办事项）
    
    /// \brief 拥有最大相对残差maxRelRes_V的块(block)的索引
    OCP_INT maxId_V;
    
    /// \brief 拥有最大相对残差maxRelRes_N的块(block)的索引
    OCP_INT maxId_N;
    
    /// \brief 拥有最大相对残差maxRelRes_E的块(block)的索引
    OCP_INT maxId_E;
};

#endif  /* end if __OCPNRRESIDUAL_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/
