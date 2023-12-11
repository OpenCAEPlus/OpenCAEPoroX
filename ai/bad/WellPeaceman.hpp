/*! \file    WellPeaceman.hpp 
 *  \brief   WellPeacema class declaration 
 *  \author  Shizhe Li 
 *  \date    Aug/17/2023 
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __WELLPEACEMAN_HEADER__
#define __WELLPEACEMAN_HEADER__

// Standard header files
#include <cassert>

// OpenCAEPoroX header files
#include "Well.hpp"
#include "WellPerf.hpp"

using namespace std;

/// Peaceman井模型模块
PeacemanWell : public Well{
public:
    /// 输入射孔信息
    void InputPerfo(const WellParam& well, const Domain& domain, const USI& wId) override;

    /// 组装井
    void Setup(const Bulk& bk, const vector<SolventINJ>& sols) override;
    
    /// 初始化井压力
    void InitWellP(const Bulk& bk) override;

    /// 检查井控制模式
    void CheckOptMode(const Bulk& bk) override;

    /// 时间步开始计算井流量
    void CalFluxInit(const Bulk& bk) override;

    /// 计算井流量
    void CalFlux(const Bulk& bk) override;

    /// 检查井压力
    ReservoirState CheckP(const Bulk& bk) override;

    /// 计算产量注入量
    void CalIPRate(const Bulk& bk, const OCP_DBL& dt) override;

    /// 计算时间步内最大变化
    OCP_DBL CalMaxChangeTime() const override;

    /// 计算牛顿步内最大变化
    OCP_DBL CalMaxChangeNR() override;

    /// 重置为上一个时间步井状态
    void ResetToLastTimeStep(const Bulk& bk) override;

    /// 保存当前时间步井状态
    void UpdateLastTimeStep() override;

protected:
    /// 计算井指数
    void CalWI(const Bulk& bk);

    /// 计算井传导率
    void CalTrans(const Bulk& bk);

    /// 计算井流量
    void CalFlux(const Bulk& bk, const OCP_BOOL ReCalXi);

    /// 计算最大注入速率
    OCP_DBL CalInjRateMaxBHP(const Bulk& bk);

    /// 计算最大注入速率
    OCP_DBL CalProdRateMinBHP(const Bulk& bk);

    /// 计算注入流量
    void CalInjQj(const Bulk& bk, const OCP_DBL& dt);

    /// 计算生产流量
    void CalProdQj(const Bulk& bk, const OCP_DBL& dt);

    /// 检查交叉流
    ReservoirState CheckCrossFlow(const Bulk& bk);

    /// 计算生产注入因子
    void CalFactor(const Bulk& bk) const;

    /// 计算射孔压差
    void CaldG(const Bulk& bk);

    /// 计算注入井射孔压差
    void CalInjdG(const Bulk& bk);

    /// 计算生产井射孔压差
    void CalProddG(const Bulk& bk);

    /// 计算生产井射孔压差方法1
    void CalProddG01(const Bulk& bk);

    /// 计算生产井射孔压差方法2
    void CalProddG02(const Bulk& bk);

    /// 计算射孔压力
    void CalPerfP() { for (USI p = 0; p < numPerf; p++) perf[p].P = bhp + dG[p]; }

protected:

    vector<OCP_DBL> dG;               ///< 射孔压差
    mutable vector<OCP_DBL> factor;   ///< 注入生产因子
};

class PeacemanWellIsoT : public PeacemanWell
{
public:
    /// FIM方法计算井残差
    void CalResFIM(OCP_USI& wId, OCPNRresidual& res, const Bulk& bk, const OCP_DBL& dt) const override;
    /// FIM方法获取井解
    void GetSolutionFIM(const vector<OCP_DBL>& u, OCP_USI& wId) override;
    /// FIM方法装配井矩阵
    void AssembleMatFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const override;
protected:
    /// FIM方法装配注入井矩阵
    void AssembleMatInjFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
    /// FIM方法装配生产井矩阵
    void AssembleMatProdFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;

public:
    /// IMPEC方法获取井解
    void GetSolutionIMPEC(const vector<OCP_DBL>& u, OCP_USI& wId) override;
    /// IMPEC方法装配井矩阵
    void AssembleMatIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const override;
protected:
    /// IMPEC方法装配注入井矩阵
    void AssembleMatInjIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
    /// IMPEC方法装配生产井矩阵
    void AssembleMatProdIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;

};


class PeacemanWellT : public PeacemanWell
{
public:
    /// FIM方法计算热流井残差
    void CalResFIM(OCP_USI& wId, OCPNRresidual& res, const Bulk& bk, const OCP_DBL& dt) const override;
    /// FIM方法获取热流井解
    void GetSolutionFIM(const vector<OCP_DBL>& u, OCP_USI& wId) override;
    /// FIM方法装配热流井矩阵
    void AssembleMatFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const override;
protected:
    /// FIM方法装配热流注入井矩阵
    void AssembleMatInjFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;
    /// FIM方法装配热流生产井矩阵
    void AssembleMatProdFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const;

public:
    /// 不可用
    void GetSolutionIMPEC(const vector<OCP_DBL>& u, OCP_USI& wId) override { OCP_ABORT("NOT USED!"); }
    /// 不可用
    void AssembleMatIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const override { OCP_ABORT("NOT USED!"); }
};

#endif /* end if __WELLPEACEMAN_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/17/2023      Create file                          */
/*----------------------------------------------------------------------------*/