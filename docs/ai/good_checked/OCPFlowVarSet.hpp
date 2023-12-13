/*! \file    OCPFlowVarSet.hpp
 *  \brief   OCPFlowVarSet类的声明
 *  \author  Shizhe Li
 *  \date    Oct/04/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPFLOWVARSET_HEADER__
#define __OCPFLOWVARSET_HEADER__

#include "OCPConst.hpp"
#include <vector>

using namespace std;

/// 流体类型枚举
enum class OCPFlowType : USI {
    SP,     ///< 单相
    OG,     ///< 油气两相
    OW,     ///< 油水两相
    GW,     ///< 气水两相
    OGW     ///< 油气水三相
};

/// 流体变量集合类
class OCPFlowVarSet {
    /// 以油相为参考相
public:
    OCPFlowVarSet() { Init0(); }

    /// 初始化为零的函数
    void Init0() {
        Swco = 0;
        krocw = 0;
        krog = 0; krow = 0;
        dKrogdSo = 0; dKrogdSg = 0;
        dKrowdSo = 0; dKrowdSw = 0;
    }

    /// 根据流体类型、相数和组分数初始化变量集合
    void Init(const OCPFlowType& fType, const USI& numPhase, const USI& numCom) {
        flowType = fType;
        np = numPhase;
        nc = numCom;
        switch (flowType) {
        case OCPFlowType::SP:
            o = 0;
            g = -1;
            w = -1;
            break;
        case OCPFlowType::OG:
            o = 0;
            g = 1;
            w = -1;
            break;
        case OCPFlowType::OW:
            o = 0;
            w = 1;
            g = -1;
            break;
        case OCPFlowType::GW:
            w = 0;
            g = 1;
            o = -1;
            break;
        case OCPFlowType::OGW:
            o = 0;
            g = 1;
            w = 2;
            break;
        default:
            OCP_ABORT("不可用的流体类型！");
            break;
        }
        oo = o * np + o;
        og = o * np + g;
        ow = o * np + w;
        go = g * np + o;
        gg = g * np + g;
        gw = g * np + w;
        wo = w * np + o;
        wg = w * np + g;
        ww = w * np + w;
        // 初始化变量
        Init0();
        S.resize(np);
        kr.resize(np);
        Pc.resize(np);
        dKrdS.resize(np * np);
        dPcdS.resize(np * np);
    }

public:
    OCPFlowType flowType;          ///< 流体类型
    INT o, g, w;                   ///< 油、气、水的索引
    INT oo, og, ow, go, gg, gw, wo, wg, ww; ///< 辅助索引
    USI np, nc;                    ///< 相数和组分数
    vector<OCP_DBL> S;             ///< 饱和度
    vector<OCP_DBL> kr;            ///< 相对渗透率
    vector<OCP_DBL> Pc;            ///< 比参考相的毛细管压力, Pj - Pr
    vector<OCP_DBL> dKrdS;         ///< 相对渗透率对饱和度的导数
    vector<OCP_DBL> dPcdS;         ///< 毛细管压力对饱和度的导数
    OCP_DBL Swco;                  ///< 连通水的饱和度
    OCP_DBL krocw;                 ///< 仅存在连通水时的油相相对渗透率
    OCP_DBL krog;                  ///< 存在油气和连通水时相应的油相相对渗透率
    OCP_DBL krow;                  ///< 仅存在油水时相应的油相相对渗透率
    OCP_DBL dKrogdSo, dKrogdSg;    ///< 渗透率对饱和度的导数
    OCP_DBL dKrowdSo, dKrowdSw;    ///< 渗透率对饱和度的导数
};

#endif /* end if __OCPFLOWVARSET_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/04/2023      Create file                          */
/*----------------------------------------------------------------------------*/