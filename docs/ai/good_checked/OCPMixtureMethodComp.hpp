/*! \file    OCPMixtureMethodComp.hpp
 *  \brief   OCPMixtureMethodComp 类声明
 *  \author  Shizhe Li
 *  \date    Jul/31/2023
 *
 *  本文件包含 OCPMixtureMethodComp 类的声明，该类是组成模型中使用的基础类，
 *  其中的变量涉及参与相平衡计算的组分和相，这些组分和相应当首先被排序。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPMIXTURECOMPMETHOD_HEADER__
#define __OCPMIXTURECOMPMETHOD_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPFuncPVT.hpp"
#include "OCPMixtureVarSet.hpp"
#include "OCPPhaseEquilibrium.hpp"
#include "BulkVarSet.hpp"
#include <vector>
using namespace std;

/**
 *  \class OCPMixtureMethodComp
 *  \brief 组分方法基类
 *
 *  OCPMixtureMethodComp 是在组成模型中使用的基础类，所有变量都与参与相平衡计算的组分和相有关。
 */
class OCPMixtureMethodComp
{
public:
    /**
     *  \brief 默认构造函数
     */
    OCPMixtureMethodComp() = default;

    /**
     *  \brief 设置变量集
     *  \param bId 块ID
     *  \param bvs 块变量集
     *  \param mvs 混合物变量集
     */
    virtual void SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const = 0;

    /**
     *  \brief 设置变量集
     *  \param P 压力
     *  \param T 温度
     *  \param Ni 组件摩尔数数组
     *  \param mvs 混合物变量集
     */
    virtual void SetVarSet(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, OCPMixtureVarSet& mvs) const = 0;

    /**
     *  \brief 仅执行闪蒸计算
     *  \param vs 变量集
     */
    virtual void Flash(OCPMixtureVarSet& vs) = 0;

    /**
     *  \brief 执行闪蒸计算，并仅计算 VfP, Vfi
     *  \param Vp 体积百分率
     *  \param vs 变量集
     */
    virtual void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;

    /**
     *  \brief 执行闪蒸计算，并仅计算 VfP, Vfi
     *  \param vs 变量集
     *  \param ftype 闪蒸类型
     */
    virtual void Flash(OCPMixtureVarSet& vs, const USI& ftype) = 0;

    /**
     *  \brief 执行闪蒸计算，并计算 VfP, Vfi, dXsdXp
     *  \param Vp 体积百分率
     *  \param vs 变量集
     */
    virtual void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;

    /**
     *  \brief 执行闪蒸计算，并计算 VfP, Vfi, dXsdXp
     *  \param vs 变量集
     *  \param ftype 闪蒸类型
     */
    virtual void FlashDer(OCPMixtureVarSet& vs, const USI& ftype) = 0;

    /**
     *  \brief 计算目标相的摩尔密度
     *  \param P 压力
     *  \param T 温度
     *  \param z 组件摩尔分数数组
     *  \param pt 相类型
     *  \return 目标相的摩尔密度
     */
    virtual OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;

    /**
     *  \brief 标准条件下的闪蒸计算
     *  \param vs 变量集
     */
    virtual void CalVStd(OCPMixtureVarSet& vs) = 0;

    /**
     *  \brief 计算标准条件下目标相的摩尔密度
     *  \param P 压力
     *  \param T 温度
     *  \param z 组件摩尔分数数组
     *  \param pt 相类型
     *  \return 标准条件下目标相的摩尔密度
     */
    virtual OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;

    /**
     *  \brief 计算目标相的质量密度
     *  \param P 压力
     *  \param T 温度
     *  \param z 组件摩尔分数数组
     *  \param pt 相类型
     *  \return 目标相的质量密度
     */
    virtual OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;

    /**
     *  \brief 输出模拟期间总闪蒸迭代次数
     */
    virtual void OutIters() const { PE.OutMixtureIters(); }

    /**
     *  \brief 判断当前混合物是否适合井
     *  \return 是否适合井的布尔值
     */
    virtual OCP_BOOL IfWellFriend() const = 0;

    /**
     *  \brief 输入现有相的总数，返回 PEC 中存在相的数目
     *  \param np 现有相的总数
     *  \return PEC 中存在相的数目
     */
    virtual USI GetNumPhasePE(const USI& np) const = 0;

    /**
     *  \brief 获取 PE 中组分的数量
     *  \return 组分数量
     */
    const auto& GetNC() const { return NC; }

    /**
     *  \brief 获取 PE 中相的数量
     *  \return 相数量
     */
    const auto& GetNP() const { return NP; }

    /**
     *  \brief 获取 PE 中允许的最大相数
     *  \return 允许的最大相数
     */
    const auto& GetNPmax() const { return NPmax; }

    /**
     *  \brief 获取 EoS
     *  \return EoS 指针
     */
    const auto GetEoS() const { return &eos; }

    /**
     *  \brief 获取 Ftype
     *  \return Ftype
     */
    const auto& GetFtype() const { return PE.GetFtype(); }

    /**
     *  \brief 获取 zi
     *  \return 组件摩尔分数
     */
    const auto& GetZi() const { return zi; }

    /**
     *  \brief 获取 Nt
     *  \return 组件总摩尔数
     */
    const auto& GetNt() const { return Nt; }

    // 基础变量
protected:
    USI             NPmax; ///< 允许的最大相数
    USI             NC;    ///< 组件数
    USI             NP;    ///< 现有相数
    OCP_DBL         Nt;    ///< 组件总摩尔数
    vector<OCP_DBL> zi;    ///< 组件摩尔分数

    // 基础组分属性
protected:
    vector<OCP_DBL> Tc;    ///< 组件临界温度
    vector<OCP_DBL> Pc;    ///< 组件临界压力
    vector<OCP_DBL> Vc;    ///< 组件临界体积
    vector<OCP_DBL> MWC;   ///< 组件分子量

    // EoS 和相平衡计算
protected:
    void CopyPhaseFromPE(OCPMixtureVarSet& vs); ///< 从相平衡复制相属性

protected:
    OCPPhaseEquilibrium PE; ///< 相平衡计算
    EoSCalculation      eos; ///< EoS 计算

    // 基础相属性
protected:
    void CalMW(OCPMixtureVarSet& vs); ///< 计算分子量
    void CalVmVj(OCPMixtureVarSet& vs); ///< 计算摩尔体积和相体积
    void CalProperty(OCPMixtureVarSet& vs); ///< 计算相属性
    void CalXiRhoMu(OCPMixtureVarSet& vs); ///< 计算相的摩尔密度、质量密度和粘度
    void CalPropertyDer(OCPMixtureVarSet& vs); ///< 计算相属性及其导数
    void CalXiRhoMuDer(OCPMixtureVarSet& vs); ///< 计算相的摩尔密度、质量密度和粘度及其导数

protected:
    vector<OCP_DBL>         MW;  ///< 相的分子量
    vector<OCP_DBL>         vm;  ///< 相的摩尔体积
    vector<OCP_DBL>         vmP; ///< 摩尔体积对压力的导数
    vector<vector<OCP_DBL>> vmx; ///< 摩尔体积对组分摩尔分数的导数
    ViscosityCalculation    visCal; ///< 粘度计算

    // 相识别
protected:
    void IdentifyPhase(OCPMixtureVarSet& vs); ///< 识别相
    void ReOrderPhase(OCPMixtureVarSet& vs); ///< 重排序相
    void ReOrderPhaseDer(OCPMixtureVarSet& vs); ///< 重排序相及其导数

protected:
    vector<PhaseType> phaseLabel; ///< 相标签，用于识别相
    vector<USI>       epIndex;    ///< 所有存在相的索引
    vector<OCP_DBL>   rowork;     ///< 重排序工作空间

    // 特殊导数计算
protected:
    virtual void CalVfiVfp_full01(OCPMixtureVarSet& vs); ///< 计算完整导数 dVf / dP, dVf / dNi
    void AssembleMatVfiVfp_full01(); ///< 组装矩阵用于计算完整导数
    void AssembleRhsVfiVfp_full01(); ///< 组装右手边用于计算完整导数
    virtual void CaldXsdXp01(OCPMixtureVarSet& vs); ///< 计算导数 d (Sj xij) / d (P , Ni)

    // 当相数 <=2 时的特殊导数计算
    virtual void CalVfiVfp_full02(OCPMixtureVarSet& vs); ///< 计算完整导数 dVf / dP, dVf / dNi
    void AssembleMatVfiVfp_full02(); ///< 组装矩阵用于计算完整导数
    void AssembleRhsVfiVfp_full02(); ///< 组装右手边用于计算完整导数
    virtual void CaldXsdXp02(OCPMixtureVarSet& vs); ///< 计算导数 d (Sj xij) / d (P , Ni)

protected:
    vector<vector<OCP_DBL>> lnfugP; ///< 导数 d ln fij / d P
    vector<vector<OCP_DBL>> lnfugN; ///< 导数 d ln fij / d nkj
    vector<OCP_DBL>         JmatDer; ///< 计算导数的雅可比矩阵
    vector<OCP_DBL>         rhsDer;  ///< 计算导数 d nij / d Nk, d nij / dP 以及 dXs / dXp 的右手边
    vector<OCP_INT>         pivot;   ///< 用于 lapack 线性求解的 pivot 数组
};

// OCPMixtureMethodComp01 类的其他声明和定义应在此处继续...

#endif /* end if __OCPMIXTURECOMPMETHOD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/31/2023      Create file                          */
/*----------------------------------------------------------------------------*/