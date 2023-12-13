/*! \file    SATModule.hpp
 *  \brief   SATModule类声明
 *  \author  Shizhe Li
 *  \date    Aug/21/2023
 *
 *-----------------------------------------------------------------------------------
 *  版权所有 (C) 2021--至今 OpenCAEPoroX团队。保留所有权利。
 *  根据GNU Lesser General Public License 3.0或更高版本的条款发布。
 *-----------------------------------------------------------------------------------
 */

#ifndef __SATMODULE_HEADER__
#define __SATMODULE_HEADER__

// 标准头文件
#include <cassert>

// OpenCAEPoroX头文件
#include "FlowUnit.hpp"
#include "BulkVarSet.hpp"

using namespace std;

class SATModule {
public:
    /**
     * \brief 设置函数，用于设置饱和度模块的参数
     * \param rs_param 油藏参数
     * \param bvs 模块参数
     * \param opts 可选模块
     */
    void Setup(const ParamReservoir& rs_param, const BulkVarSet& bvs, OptionalModules& opts) {
        NTSFUN = rs_param.NTSFUN;
        
        // 设置饱和度函数
        for (USI i = 0; i < NTSFUN; i++)
            SATs.push_back(FlowUnit(rs_param, i, opts));
        
        if (SATNUM.empty() || NTSFUN == 1) {
            SATNUM.clear();
            SATNUM.resize(bvs.nb, 0);
        }
    }
    
    /**
     * \brief 获取饱和度模块
     * \param n bulk的索引
     * \return 返回指向饱和度模块的指针
     */
    auto GetSAT(const OCP_USI& n) const { return &SATs[SATNUM[n]]; }
    
    /**
     * \brief 获取饱和度模块索引
     * \return 返回饱和度模块索引的引用
     */
    auto& GetSATNUM() { return SATNUM; }
    
    /**
     * \brief 获取饱和度模块索引
     * \param n bulk的索引
     * \return 返回饱和度模块索引
     */
    auto GetSATNUM(const OCP_USI& n) const { return SATNUM[n]; }
    
    /**
     * \brief 获取饱和度函数的数量
     * \return 返回饱和度函数的数量
     */
    auto GetNTSFUN() { return NTSFUN; }
    
protected:
    /// 饱和度区域的数量
    USI               NTSFUN;
    
    /// 每个bulk的饱和度区域索引
    vector<USI>       SATNUM;
    
    /// 饱和度模块
    vector<FlowUnit>  SATs;
};

#endif /* end if __SATMODULE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/21/2023      Create file                          */
/*----------------------------------------------------------------------------*/