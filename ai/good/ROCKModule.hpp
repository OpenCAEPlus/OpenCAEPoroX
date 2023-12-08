/*! \file    ROCKModule.hpp
 *  \brief   ROCK模块类的声明
 *  \author  Shizhe Li
 *  \date    Aug/21/2023
 *
 *  本文件包含ROCKModule类的声明，用于设置和获取岩石模块的相关参数。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __ROCKMODULE_HEADER__
#define __ROCKMODULE_HEADER__

// 标准头文件
#include <cassert>

// OpenCAEPoroX头文件
#include "OCPRock.hpp"
#include "BulkVarSet.hpp"
#include "OptionalModules.hpp"

using namespace std;

/*!
 * \class ROCKModule
 * \brief 岩石模块类
 *
 * 该类负责根据参数设置岩石模型，并提供接口以获取岩石模型和相关数据。
 */
class ROCKModule
{
public:
    /*!
     * \brief 设置岩石模型
     *
     * 根据油藏参数、块数量和可选模块设置岩石模型。
     * \param rs_param 油藏参数
     * \param nb 块的数量
     * \param opts 可选模块
     */
    void Setup(const ParamReservoir& rs_param, const OCP_USI& nb, OptionalModules& opts)
    {
        NTROCK = rs_param.NTROOC;
        if (rs_param.thermal) {
            for (USI i = 0; i < NTROCK; i++) {
                if (rs_param.rockSet[i].type == "LINEAR") {
                    ROCKs.push_back(new OCPRockT_Linear(rs_param.rockSet[i]));
                }
                else {
                    ROCKs.push_back(new OCPRockT_Exp(rs_param.rockSet[i]));
                }
            }
        }
        else {
            for (USI i = 0; i < NTROCK; i++) {
                ROCKs.push_back(new OCPRockIsoT_Linear(rs_param.rockSet[i]));
            }
        }
        if (ROCKNUM.empty() || NTROCK == 1) {
            ROCKNUM.clear();
            ROCKNUM.resize(nb, 0);
        }
    }

    /*!
     * \brief 获取指定索引的岩石模型
     *
     * \param n 岩石模型的索引
     * \return 返回指定索引的岩石模型对象指针
     */
    auto GetROCK(const OCP_USI& n) const { return ROCKs[ROCKNUM[n]]; }

    /*!
     * \brief 获取岩石编号的引用
     *
     * \return 返回岩石编号的向量引用
     */
    auto& GetROCKNUM() { return ROCKNUM; }

    /*!
     * \brief 获取岩石区域的数量
     *
     * \return 返回岩石区域的数量
     */
    auto GetNTROCK() { return NTROCK; }

protected:
    /// 岩石区域的数量
    USI              NTROCK;
    /// 每个块的岩石区域索引
    vector<USI>      ROCKNUM;
    /// 岩石模块的向量
    vector<OCPRock*> ROCKs;
};

#endif /* end if __ROCKMODULE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/21/2023      Create file                          */
/*----------------------------------------------------------------------------*/
