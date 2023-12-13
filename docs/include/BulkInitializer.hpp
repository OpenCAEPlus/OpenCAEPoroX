/*! \file    BulkInitializer.hpp
 *  \brief   BulkInitializer类声明
 *  \details 该文件包含了BulkInitializer类的声明。
 *  \author  Shizhe Li
 *  \date    Aug/26/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BULKINITIALIZER_HEADER__
#define __BULKINITIALIZER_HEADER__

#include "BulkVarSet.hpp"
#include "PVTModule.hpp"
#include "SATModule.hpp"
#include "Domain.hpp"
#include <vector>

using namespace std;

/*!
 *  \class Equil
 *  \brief 初始油藏信息，用于计算初始平衡。
 *  \note 该类包括参考深度和其压力，相间接触面的深度，以及相接触面的毛管压力。
 */
class Equil{
    friend class BulkInitializer;
    friend class Bulk;

protected:
    OCP_DBL  Dref;     ///< 参考深度
    OCP_DBL  Pref;     ///< 参考深度的压力
    OCP_DBL  DOWC;     ///< 油水接触面的深度
    OCP_DBL  PcOW;     ///< 油水接触面的毛管压力 Pcow = Po - Pw
    OCP_DBL  DGOC;     ///< 气油接触面的深度
    OCP_DBL  PcGO;     ///< 气油接触面的毛管压力 Pcgo = Pg - Po
    OCPTable PBVD;     ///< PBVD表：气泡点压力与深度
};

/*!
 *  \class BulkInitializer
 *  \brief 初始化大体
 */
class BulkInitializer{
public:
	void Setup(const ParamReservoir& rs_param, const OCPMixtureType& mixType);
    void Initialize(BulkVarSet& bvs, const PVTModule& PVTm, const SATModule& SATm, const OptionalModules& optMs, const Domain& domain);
    auto& GetSwat() { return swat; } ///< 获取初始油藏水饱和度

protected:
    void InitHydroEquil(BulkVarSet& bvs, const PVTModule& PVTm, const SATModule& SATm, const Domain& domain); ///< 仅用水静力平衡初始化油藏
    void InitHydroEquilW(BulkVarSet& bvs, const PVTModule& PVTm, const SATModule& SATm, const OptionalModules& optMs, const Domain& domain); ///< 用水和水静力平衡初始化油藏

protected:
    vector<Equil>    EQUIL;    ///< 平衡数据规范
    vector<OCPTable> initZi_Tab;     ///< 初始组分摩尔比与深度，表集
    vector<OCPTable> initT_Tab;    ///< 初始温度与深度，表集
    USI              numNodes{ 50 };    ///< P vs. Z表的节点数
    vector<OCP_DBL>  swat;    ///< 初始油藏水饱和度
};

#endif /* end if __BULKINITIALIZER_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/26/2023      Create file                          */
/*----------------------------------------------------------------------------*/
