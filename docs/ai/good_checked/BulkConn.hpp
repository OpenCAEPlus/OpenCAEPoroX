/*! \file    BulkConn.hpp 
 *  \brief   BulkConn类的声明
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  版权所有 (C) 2021--至今 OpenCAEPoroX团队。保留所有权利。
 *  根据GNU Lesser General Public License 3.0或更高版本的条款发布。
 *-----------------------------------------------------------------------------------
 */

#ifndef __BULKCONN_HEADER__
#define __BULKCONN_HEADER__

// OpenCAEPoroX头文件
#include "Bulk.hpp"
#include "OCPFlux.hpp"

using namespace std;

/// 关于bulk之间连接的属性和操作（活动网格之间的连接）。
//  注意：BulkConn是储层的核心组件，它包含了关于bulk之间连接的所有属性和操作。
//  通过一个有效的迭代器，您可以遍历所有的连接。
//  包括活动bulk之间的流量计算和仅来自bulk的矩阵装配的贡献。
class BulkConn{
    friend class Reservoir;
    friend class NRsuite;
    // temp
    friend class IsoT_FIM;
    friend class IsoT_IMPEC;
    friend class IsoT_AIMc;
    friend class IsoT_FIMddm;
    friend class T_FIM;

public:
    /////////////////////////////////////////////////////////////////////
    // 通用变量
    /////////////////////////////////////////////////////////////////////

public:
    /// 输入参数
    void InputParam(const ParamReservoir& rs_param, const Bulk& bk);

    /// 获取变量集合
    auto& GetVarSet() const { return vs; }

    /// 获取连接数量
    auto GetNumConn() const { return numConn; }

protected:
    OCP_USI numConn; ///< bulk之间连接的数量

    /// 所有bulk之间的连接（索引对）：numConn。
    //  注意：在每对中，第一个bulk的索引小于第二个bulk的索引。
    //  iteratorConn中的数据是从neighbor生成的。
    vector<BulkConnPair> iteratorConn;

    /// 每个bulk的邻居信息：activeGridNum。
    //  注意：第i个条目存储第i个bulk的邻居，按升序排序。
    vector<vector<OCP_USI>> neighbor;

    /////////////////////////////////////////////////////////////////////
    // 物理变量
    /////////////////////////////////////////////////////////////////////

    BulkConnVarSet          vs;   ///< 连接之间的值

    /////////////////////////////////////////////////////////////////////
    // 流量
    /////////////////////////////////////////////////////////////////////

    vector<OCPFlux*>        flux;    ///< 流量项
};

#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/17/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/