/*! \file    BulkConnVarSet.hpp
 *  \brief   BulkConnVarSet类的声明
 *  \author  Shizhe Li
 *  \date    Aug/23/2023
 *
 *  本文件包含BulkConnVarSet类及相关辅助类的声明，用于表示和处理储层网格单元之间的连接。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BULKCONNVARSET_HEADER__
#define __BULKCONNVARSET_HEADER__

// OpenCAEPoroX头文件
#include "OCPConst.hpp"
#include <vector>
using namespace std;

/*!
 * \class BulkConnPair
 * \brief 代表两个储层网格单元之间的连接。
 *
 * BulkConnPair类用于存储两个储层网格单元之间的连接信息，通常索引bId大于eId。
 * 注意：网格单元是活跃的网格单元。
 */
class BulkConnPair
{
    friend class BulkConnTransMethod01;

public:
    /// 默认构造函数。
    BulkConnPair() = default;

    /*!
     * \brief 使用bId和eId初始化BulkPair。
     *
     * \param BId        起始网格单元索引。
     * \param EId        结束网格单元索引。
     * \param direct     连接方向。
     * \param AreaB      起始网格单元的交叉面积。
     * \param AreaE      结束网格单元的交叉面积。
     * \param transM     渗透率乘数。
     */
    BulkConnPair(const OCP_USI& BId,
                 const OCP_USI& EId,
                 const ConnDirect& direct,
                 const OCP_DBL& AreaB,
                 const OCP_DBL& AreaE,
                 const OCP_DBL& transM) 
        : bId(BId)
        , eId(EId)
        , direction(direct)
        , areaB(AreaB)
        , areaE(AreaE)
        , transMult(transM) {};

    void SetFluxNum(const USI& i) { fluxnum = i; }
    const auto& BId() const { return bId; }
    const auto& EId() const { return eId; }
    const auto& Trans() const { return trans; }
    const auto& FluxNum() const { return fluxnum; }
    const auto& AreaB() const { return areaB; }
    const auto& AreaE() const { return areaE; }
    const auto& Direction() const { return direction; }

protected:
    OCP_USI    bId;     //!< 第一个网格单元的索引。
    OCP_USI    eId;     //!< 第二个网格单元的索引。
    OCP_DBL    trans;   //!< 有效面积。
    USI        fluxnum{ 0 }; //!< 连接类型。
    ConnDirect direction; //!< 连接方向。
    OCP_DBL    areaB;   //!< 第一个网格单元的交叉面积。
    OCP_DBL    areaE;   //!< 第二个网格单元的交叉面积。
    OCP_DBL    transMult{ 1.0 }; //!< 渗透率乘数。
};

/*!
 * \class BulkConnVarSet
 * \brief 用于存储储层网格单元连接的变量集合。
 *
 * BulkConnVarSet类包含了一系列向量，这些向量存储了与储层网格单元连接相关的物理和流体动力学参数。
 */
class BulkConnVarSet
{
    /////////////////////////////////////////////////////////////////////
    // 一般信息
    /////////////////////////////////////////////////////////////////////
public:
    vector<OCP_USI> upblock; //!< 每个相的上行网格单元索引。
    vector<OCP_DBL> rho;     //!< 连接处每个相的质量密度。
    vector<OCP_DBL> flux_vj; //!< 从上行网格单元的相体积流量。
    vector<OCP_DBL> flux_ni; //!< 组件的摩尔流量。
    vector<OCP_USI> lupblock; //!< 最后的上行网格单元。
    vector<OCP_DBL> lrho;     //!< 最后的质量密度。
};

#endif /* end if __BulkConnVarSet_HEADER__ */

/*----------------------------------------------------------------------------
 *  Brief Change History of This File
 *----------------------------------------------------------------------------
 *  Author              Date             Actions
 *----------------------------------------------------------------------------
 *  Shizhe Li           Aug/23/2023      Create file
 *----------------------------------------------------------------------------
 */