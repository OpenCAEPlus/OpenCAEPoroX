/*! \file    OCPMiscible.hpp
 *  \brief   OCPMiscible文件包含了与石油工程中的可混溶性相关的类的声明。
 *  \author  Shizhe Li
 *  \date    Dec/26/2022
 *
 *  本文件定义了与可混溶性因子计算、曲线校正相关的类。这些类用于模拟油气藏工程中流体的可混溶行为，
 *  包括渗透率和毛细管压力的校正。
 */

#ifndef __OCPMISCIBLE_HEADER__
#define __OCPMISCIBLE_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPSurfaceTension.hpp"
#include "OCPFlow.hpp"
#include <vector>
using namespace std;

/*! \class MisFacVarSet
 *  \brief 用于存储可混溶性因子计算结果的类。
 */
class MisFacVarSet
{
public:
    void SetNb(const OCP_USI& nbin) { nb = nbin; }
    
public:
    /// \brief 储存区块数量。
    OCP_USI         nb;
    /// \brief 渗透率的可混溶性因子。
    vector<OCP_DBL> Fk;
    /// \brief 毛细管压力的可混溶性因子。
    vector<OCP_DBL> Fp;
};

/*! \class MisFacMethod
 *  \brief 可混溶性因子计算方法的抽象基类。
 */
class MisFacMethod
{
public:
    MisFacMethod() = default;
    virtual void CalculateMiscibleFactor(const OCP_USI& bId, const SurTenVarSet& stvs, MisFacVarSet& mfvs) const = 0;
};

/*! \class MisFacMethod01
 *  \brief Coats表达式的实现，用于计算可混溶性因子。
 */
class MisFacMethod01 : public MisFacMethod
{
public:
    MisFacMethod01(const Miscstr& param, MisFacVarSet& mfvs);
    void CalculateMiscibleFactor(const OCP_USI& bId, const SurTenVarSet& stvs, MisFacVarSet& mfvs) const override;

public:
    /// \brief 参考表面张力。
    OCP_DBL stref;
    /// \brief Fk计算的指数。
    OCP_DBL fkExp;
    /// \brief 毛细管压力计算的最大表面张力与参考表面张力的比值。
    OCP_DBL stPcf;
};

/*! \class MiscibleFactor
 *  \brief 用于管理可混溶性因子计算和存储的类。
 */
class MiscibleFactor
{
public:
    /// \brief 默认构造函数。
    MiscibleFactor() = default;
    USI Setup(const ParamReservoir& param, const USI& i, const OCP_USI& nb, const SurfaceTension* st);
    void CalMiscibleFactor(const OCP_USI& bId, const USI& mIndex) {
        if (ifUse) {
            mfMethod[mIndex]->CalculateMiscibleFactor(bId, surTen->GetVS(), vs);
        }       
    }
    const auto& GetVarSet() const { return vs; }
    const auto IfUse()const { return ifUse; }
    const auto GetFk(const OCP_USI& bId) const { return vs.Fk[bId]; }
    const auto GetFp(const OCP_USI& bId) const { return vs.Fp[bId]; }
    void ResetTolastTimeStep() { }
    void UpdateLastTimeStep() { }

protected:
    /// \brief 是否计算可混溶性因子。
    OCP_BOOL              ifUse{ OCP_FALSE };
    /// \brief 可混溶性因子计算变量集。
    MisFacVarSet          vs;
    /// \brief 可混溶性因子计算方法。
    vector<MisFacMethod*> mfMethod;
    // 依赖模块
    /// \brief 表面张力计算。
    const SurfaceTension* surTen;
};

/*! \class MisCurveMethod
 *  \brief 渗透率和毛细管压力曲线校正方法的抽象基类。
 */
class MisCurveMethod
{
public:
    MisCurveMethod() = default;
    virtual void CurveCorrect(const OCP_USI& bId, const MisFacVarSet& mfvs) = 0;
    virtual void CurveCorrectDer(const OCP_USI& bId, const MisFacVarSet& mfvs) = 0;
};

/*! \class MisCurveMethod01
 *  \brief CMG中的SIGMA方法的实现，用于曲线校正。
 */
class MisCurveMethod01 : public MisCurveMethod
{
public:
    MisCurveMethod01() = default;
    MisCurveMethod01(OCPFlow* flowin) {
        switch (flowin->FlowType())
        {
        case OCPFlowType::OG:
        case OCPFlowType::OGW:
            break;
        default:
            OCP_ABORT("Wrong FlowType for ScalePcow!");
            break;
        }

        flow = flowin;
    }
    void CurveCorrect(const OCP_USI& bId, const MisFacVarSet& mfvs) override;
    void CurveCorrectDer(const OCP_USI& bId, const MisFacVarSet& mfvs) override;

protected:
    // 依赖模块
    /// \brief 流动计算。
    OCPFlow* flow;
};

/*! \class MiscibleCurve
 *  \brief 用于管理可混溶曲线校正的类。
 */
class MiscibleCurve
{
public:
    MiscibleCurve() = default;
    USI Setup(OCPFlow* flowin, const MiscibleFactor* misfactor);
    void CorrectCurve(const OCP_USI& bId, const USI& mIndex) {
        if (ifUse) {
            mcMethod[mIndex]->CurveCorrect(bId, misFac->GetVarSet());
        }
    }
    void CorrectCurveDer(const OCP_USI& bId, const USI& mIndex) {
        if (ifUse) {
            mcMethod[mIndex]->CurveCorrectDer(bId, misFac->GetVarSet());
        }
    }
    void ResetTolastTimeStep() { }
    void UpdateLastTimeStep() { }

protected:
    /// \brief 是否使用可混溶曲线校正。
    OCP_BOOL                ifUse{ OCP_FALSE };
    /// \brief 可混溶曲线校正方法。
    vector<MisCurveMethod*> mcMethod;
    // 依赖模块
    /// \brief 可混溶性因子。
    const MiscibleFactor*   misFac;
};

#endif /* end if __OCPMISCIBLE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/26/2022      Create file                          */
/*  Shizhe Li           Jul/02/2022      Rename to OCPMiscible                */
/*----------------------------------------------------------------------------*/