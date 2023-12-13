/*! \file    Well.hpp
 *  \brief   本文件包含Well类的声明
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *  本文件定义了Well类，包含了与油井相关的所有操作。
 *  Well类通过射孔连接到储层单元，射孔作为资源和接收端。
 *  由于生产实践中的困难，对油井进行良好处理非常重要，
 *  优秀的处理将使井中的流量更加稳定。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __WELL_HEADER__
#define __WELL_HEADER__

// 标准头文件
#include <cassert>

// OpenCAEPoroX头文件
#include "Bulk.hpp"
#include "DenseMat.hpp"
#include "LinearSystem.hpp"
#include "OCPConst.hpp"
#include "ParamWell.hpp"
#include "WellOpt.hpp"
#include "WellPerf.hpp"
#include "OCPMixture.hpp"
#include "OCPNRresidual.hpp"

using namespace std;

/**
 * \class Well
 * \brief 定义了油井及其相关操作。
 *
 * Well类通过射孔与储层单元连接，射孔作为资源和接收端。
 * 在生产实践中，由于存在诸多困难，对油井进行良好的处理非常重要，
 * 以确保井内流量的稳定性。
 */
class Well
{
    friend class AllWells;
    friend class Out4RPT;

public:
    Well() = default;

    /////////////////////////////////////////////////////////////////////
    // 输入参数和设置
    /////////////////////////////////////////////////////////////////////

public:
    /**
     * \brief 输入射孔参数。
     * \param well 射孔参数
     * \param domain 域
     * \param wId 井ID
     */
    virtual void InputPerfo(const WellParam& well, const Domain& domain, const USI& wId) = 0;

    /**
     * \brief 在Grid和Bulk设置完成后设置井。
     * \param bk Bulk对象
     * \param sols 溶剂注入信息数组
     */
    virtual void Setup(const Bulk& bk, const vector<SolventINJ>& sols) = 0;

    /**
     * \brief 应用第i次操作。
     * \param i 操作次数
     * \return 操作是否成功
     */
    OCP_BOOL ApplyOpt(const USI& i);

protected:
    /// 设置井操作
    void SetupOpts(const vector<SolventINJ>& sols);
    /// 设置注入井操作
    void SetupOptsInj(WellOpt& opt, const vector<SolventINJ>& sols);
    /// 设置生产井操作
    void SetupOptsProd(WellOpt& opt);

public:
    /// 初始化井压力
    virtual void InitWellP(const Bulk& bk) = 0;
    /// 检查井操作模式是否改变
    virtual void CheckOptMode(const Bulk& bk) = 0;
    /// 计算通量并初始化在此时间步骤中不会改变的值
    virtual void CalFluxInit(const Bulk& bk) = 0;
    /// 计算通量
    virtual void CalFlux(const Bulk& bk) = 0;
    /// 检查是否发生异常压力
    virtual ReservoirState CheckP(const Bulk& bk) = 0;
    /// 计算注入井和生产井的相的摩尔流量
    virtual void CalIPRate(const Bulk& bk, const OCP_DBL& dt) = 0;
    /// 计算两个时间步之间井压力的最大变化
    virtual OCP_DBL CalMaxChangeTime() const = 0;
    /// 计算两个NR步之间井压力的最大变化
    virtual OCP_DBL CalMaxChangeNR() = 0;
    /// 重置到上一个时间步
    virtual void ResetToLastTimeStep(const Bulk& bk) = 0;
    /// 更新上一个时间步
    virtual void UpdateLastTimeStep() = 0;

public:
    /// 显示井操作模式和射孔状态
    void     ShowPerfStatus(const Bulk& bk) const;
    /// 获取射孔数量
    USI      PerfNum() const { return numPerf; }
    /// 返回井的状态，开启或关闭
    OCP_BOOL IsOpen() const { return (opt.state == WellState::open); }
    /// 返回井的类型，注入或生产
    auto     WellType() const { return opt.type; }
    /// 返回射孔状态
    auto     PerfState(const USI& p) const { return perf[p].state; }
    /// 返回射孔位置（储层索引）
    USI      PerfLocation(const USI& p) const { return perf[p].location; }
    /// 返回从射孔j的组分i的摩尔流量
    OCP_DBL  PerfQi_lbmol(const USI& p, const USI& i) const { return perf[p].qi_lbmol[i]; }
    /// 返回从射孔j的组分i的体积流量
    OCP_DBL  PerfProdQj_ft3(const USI& p, const USI& j) const { return perf[p].qj_ft3[j]; }

protected:
    /// 井名称
    string              name;
    /// 井所属组，如果需要，可以移到opt中
    string              group;
    /// 井的参考深度
    OCP_DBL             depth;
    /// 参考深度处的井压力
    OCP_DBL             bhp;
    /// 井控制参数，包含当前控制参数
    WellOpt             opt;
    /// 井控制参数集合
    vector<WellOpt>     optSet;
    /// 射孔数量
    USI                 numPerf;
    /// 射孔信息
    vector<Perforation> perf;
    /// 井表面压力，psia
    OCP_DBL             Psurf;
    /// 井表面温度，F
    OCP_DBL             Tsurf;
    /// 混合模型
    OCPMixture*         mixture;
    /// 相数量
    USI                 np;
    /// 组分数量
    USI                 nc;
    /// 储层温度
    OCP_DBL             rsTemp;
    /// 组件流入/流出的摩尔流量
    vector<OCP_DBL>     qi_lbmol;
    /// 是否使用未加权渗透率
    OCP_BOOL            ifUseUnweight{ OCP_FALSE };
    /// 上一个时间步的bhp
    OCP_DBL             lbhp;
    /// 上一个NR步的bhp
    OCP_DBL             NRbhp;
    /// 井石油生产率
    OCP_DBL             WOPR{0};
    /// 井石油总生产量
    OCP_DBL             WOPT{0};
    /// 井天然气生产率
    OCP_DBL             WGPR{0};
    /// 井天然气总生产量
    OCP_DBL             WGPT{0};
    /// 井水生产率
    OCP_DBL             WWPR{0};
    /// 井水总生产量
    OCP_DBL             WWPT{0};
    /// 井天然气注入率
    OCP_DBL             WGIR{0};
    /// 井天然气总注入量
    OCP_DBL             WGIT{0};
    /// 井水注入率
    OCP_DBL             WWIR{0};
    /// 井水总注入量
    OCP_DBL             WWIT{0};

    /////////////////////////////////////////////////////////////////////
    // 单位
    /////////////////////////////////////////////////////////////////////

protected:
    /// 设置单位
    void SetupUnit();
    /// 储层单位 -> 表面单位，针对特定相
    OCP_DBL UnitConvertR2S(const PhaseType& pt, const OCP_DBL& val) const;
    /// 储层单位 -> 表面单位，针对特定相索引
    OCP_DBL UnitConvertR2S(const USI& j, const OCP_DBL& val) const;

protected:
    /// 单位转换 (m3,ft3 -> m3,stb,Mscf)
    vector<OCP_DBL> unitConvert;

    /////////////////////////////////////////////////////////////////////
    // 方法
    /////////////////////////////////////////////////////////////////////

public:
    /// 为FIM方法计算残差
    virtual void CalResFIM(OCP_USI& wId, OCPNRresidual& res, const Bulk& bk, const OCP_DBL& dt) const = 0;
    /// 为FIM方法组装矩阵
    virtual void AssembleMatFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const = 0;
    /// 获取FIM方法的解
    virtual void GetSolutionFIM(const vector<OCP_DBL>& u, OCP_USI& wId) = 0;
    /// 为IMPEC方法组装矩阵
    virtual void AssembleMatIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const = 0;
    /// 获取IMPEC方法的解
    virtual void GetSolutionIMPEC(const vector<OCP_DBL>& u, OCP_USI& wId) = 0;
};

#endif /* end if __WELL_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/