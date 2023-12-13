/*! \file    WellPerf.hpp
 *  \brief   本文件包含了Perforation类的声明。
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *  本文件定义了Perforation类，该类描述了油井与储层之间的连接关系。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PERFORATION_HEADER__
#define __PERFORATION_HEADER__

// 标准头文件
#include <vector>

// OpenCAEPoroX头文件
#include "OCPConst.hpp"
#include "WellOpt.hpp"

using namespace std;

/// \enum PerfDirection
/// \brief 描述井筒穿透网格块的方向。
enum class PerfDirection : USI {
    x,              ///< x方向
    y,              ///< y方向
    z,              ///< z方向
    usg             ///< 非结构化网格
};

/// \class Perforation
/// \brief 描述油井与储层之间的连接关系。
class Perforation {
public:
    /// 默认构造函数。
    Perforation() = default;

    /// \brief 设置射孔的状态。
    /// \param flag 射孔的状态。
    void SetState(const WellState& flag) { state = flag; };

    /// \brief 获取射孔的位置：储层的索引。
    /// \return 返回储层的索引。
    OCP_USI Location() const { return location; }

public:
    WellState state{ WellState::close }; ///< 射孔状态，默认为关闭。
    OCP_USI location;                    ///< 连接到当前射孔的储层的索引。
    OCP_DBL depth;                       ///< 连接到当前射孔的储层的深度。
    OCP_DBL P;                           ///< 当前射孔的压力。
    OCP_DBL WI;                          ///< 连接的可传递性因子，可以由用户直接提供。
    OCP_DBL radius;                      ///< 井的半径。
    OCP_DBL kh;                          ///< 连接的有效渗透率乘以净厚度。
    OCP_DBL skinFactor;                  ///< 皮肤因子。
    PerfDirection direction;             ///< 井穿透网格块的方向。
    OCP_DBL multiplier;                  ///< 当前射孔的可传递性乘数。目前等于0（关闭）或1（打开）。
    mutable OCP_DBL xi;                  ///< 当前射孔中流体的摩尔密度。用于注入井，其中流体只含有单相。
    vector<OCP_DBL> qi_lbmol;            ///< 从当前射孔进出的组分的摩尔流量。
    vector<OCP_DBL> transj;              ///< 当前射孔中相的可传递性。
    OCP_DBL transINJ;                    ///< 注入可传递性。
    vector<OCP_DBL> qj_ft3;              ///< 从当前射孔进出的相的体积流量。
    OCP_DBL qt_ft3;                      ///< 从当前射孔进出的流体的体积流量。
};

#endif /* end if __PERFORATION_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/