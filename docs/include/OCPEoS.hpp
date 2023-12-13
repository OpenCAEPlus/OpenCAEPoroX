/*! \file    OCPEoS.hpp
 *  \brief   OCPEoS 类声明和相关的状态方程计算功能
 *  \author  Shizhe Li
 *  \date    Jul/23/2023
 *
 *  本文件包含了OCPEoS类的声明，以及与状态方程（Equation of State, EoS）相关的计算功能。
 *  这些功能包括计算逸度、摩尔体积等热力学性质，以及它们对应的导数。这些计算对于模拟
 *  石油工程中的流体行为至关重要。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPEOS_HEADER__
#define __OCPEOS_HEADER__

#include "OCPConst.hpp"
#include "UtilMath.hpp"
#include "DenseMat.hpp"
#include "ParamReservoir.hpp"
#include <vector>

using namespace std;

/**
 * \class EoS
 * \brief 抽象基类，定义了状态方程必须实现的接口
 *
 * EoS类是一个抽象基类，为不同的状态方程提供了统一的接口。通过继承这个类，
 * 可以实现具体的状态方程，如Peng-Robinson方程。
 */
class EoS{
public:
    EoS() = default; ///< 默认构造函数

    /// 计算逸度
    virtual void CalFug(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* fug) const = 0;

    /// 计算逸度和逸度系数
    virtual void CalFugPhi(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* fug, OCP_DBL* phi) const = 0;

    /// 计算d(lnfug) / dx
    virtual void CalLnFugX(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* lnfugx) const = 0;

    /// 计算d(lnfug) / dn
    virtual void CalLnFugN(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, const OCP_DBL& nt, OCP_DBL* lnfugn) const = 0;

    /// 计算d(lnfug) / dP
    virtual void CalLnFugP(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* lnfugP) const = 0;

    /// 计算d(lnphi) / dn
    virtual void CalLnPhiN(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, const OCP_DBL& nt, OCP_DBL* lnphin) const = 0;

public:
    /// 计算摩尔体积
    virtual OCP_DBL CalVm(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x) const = 0;

    /// 计算摩尔体积和其导数
    virtual OCP_DBL CalVmDer(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL& vmP, OCP_DBL* vmx) const = 0;
};

/////////////////////////////////////////////////
// EoS_PR
/////////////////////////////////////////////////

/**
 * \class EoS_PR
 * \brief 采用Peng-Robinson状态方程的具体实现类
 *
 * EoS_PR类继承自EoS类，实现了Peng-Robinson状态方程的具体计算方法。
 * 它包括了计算逸度、摩尔体积及其导数等热力学性质的功能。
 */
class EoS_PR : public EoS{
public:
    EoS_PR(const ComponentParam& param, const USI& tarId); ///< 构造函数，初始化状态方程需要的参数

    // 以下方法实现了EoS类中定义的虚函数
    void CalFug(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* fug) const override;
    void CalFugPhi(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* fug, OCP_DBL* phi) const override;
    void CalLnFugX(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* lnfugx) const override;
    void CalLnFugN(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, const OCP_DBL& nt, OCP_DBL* lnfugn) const override;
    void CalLnFugP(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* lnfugP) const override;
    void CalLnPhiN(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, const OCP_DBL& nt, OCP_DBL* lnphin) const override;

public:
    // 以下方法实现了EoS类中定义的虚函数
    OCP_DBL CalVm(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x) const override;
    OCP_DBL CalVmDer(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL& vmP, OCP_DBL* vmx) const override;

protected:
    // 以下方法用于计算状态方程中的各种中间变量
    void CalAiBi(const OCP_DBL& P, const OCP_DBL& T) const;
    void CalAjBj( const OCP_DBL* x) const;
    void CalZj(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x) const;
    void CalAjBjZj(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x) const;
    void CalAxBxZx(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x) const;
    void CalAnBnZn(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, const OCP_DBL& nt) const;
    void CalApBpZp(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x) const;

protected:
    USI             nc; ///< 组件数目
    vector<OCP_DBL> Pc; ///< 组件的临界压力
    vector<OCP_DBL> Tc; ///< 组件的临界温度
    vector<OCP_DBL> OmegaA; ///< 组件的OMEGA_A
    vector<OCP_DBL> OmegaB; ///< 组件的OMEGA_B
    vector<OCP_DBL> acf; ///< 组件的偏心因子
    vector<OCP_DBL> BIC; ///< 组件间的二元相互作用系数
    vector<OCP_DBL> Vshift; ///< 组件的体积偏移

    // PR-EoS 参数
    const OCP_DBL delta1 = 2.41421356237; ///< PR-EoS 参数 delta1
    const OCP_DBL delta2 = -0.41421356237; ///< PR-EoS 参数 delta2

    // 辅助变量
    mutable vector<OCP_DBL> Ai; ///< 组件的辅助变量 A
    mutable vector<OCP_DBL> Bi; ///< 组件的辅助变量 B
    mutable OCP_DBL         Aj; ///< 辅助变量 Aj
    mutable OCP_DBL         Bj; ///< 辅助变量 Bj
    mutable OCP_DBL         Zj; ///< 可压缩因子 Zj

    // 根容器
    mutable vector<OCP_DBL> Z; ///< 可压缩因子的根

    // 导数
    mutable vector<OCP_DBL> Ax; ///< dAj / dx 的导数
    mutable vector<OCP_DBL> Bx; ///< dBj / dx 的导数
    mutable vector<OCP_DBL> Zx; ///< dZj / dx 的导数
    mutable vector<OCP_DBL> An; ///< dAj / dn 的导数
    mutable vector<OCP_DBL> Bn; ///< dBj / dn 的导数
    mutable vector<OCP_DBL> Zn; ///< dZj / dn 的导数
    mutable OCP_DBL         Ap; ///< dAj / dP 的导数
    mutable OCP_DBL         Bp; ///< dBj / dP 的导数
    mutable OCP_DBL         Zp; ///< dZj / dP 的导数
};

/**
 * \class EoSCalculation
 * \brief 提供状态方程计算的功能类
 *
 * EoSCalculation类包装了状态方程的对象，提供了计算逸度、摩尔体积及其导数等功能的接口。
 * 通过这个类可以方便地进行状态方程相关的计算。
 */
class EoSCalculation{
public:
    EoSCalculation() = default; ///< 默认构造函数

    void Setup(const ComponentParam& param, const USI& tarId); ///< 设置状态方程计算需要的参数

    // 以下方法封装了状态方程对象的相关计算功能
    void CalFug(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* fug) const {
        eos->CalFug(P, T, x, fug);
    }

    void CalFugPhi(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* fug, OCP_DBL* phi) const {
        eos->CalFugPhi(P, T, x, fug, phi);
    }

    void CalLnFugX(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* lnfugx) const {
        eos->CalLnFugX(P, T, x, lnfugx);
    }

    void CalLnFugN(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x,  const OCP_DBL& nt, OCP_DBL* lnfugn) const {
        eos->CalLnFugN(P, T, x, nt, lnfugn);
    }

    void CalLnFugP(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL* lnfugP) const {
        eos->CalLnFugP(P, T, x, lnfugP);
    }

    void CalLnPhiN(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, const OCP_DBL& nt, OCP_DBL* lnphin) const {
        eos->CalLnPhiN(P, T, x, nt, lnphin);
    }

public:
    OCP_DBL CalVm(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x) const {
        return eos->CalVm(P, T, x);
    }

    OCP_DBL CalVmDer(const OCP_DBL& P, const OCP_DBL& T,  const OCP_DBL* x, OCP_DBL& vmP, OCP_DBL* vmx) const {
        return eos->CalVmDer(P, T, x, vmP, vmx);
    }

protected:
    EoS*  eos; ///< 状态方程对象指针
};

#endif /* end if __OCPEOS_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/23/2023      Create file                          */
/*----------------------------------------------------------------------------*/
