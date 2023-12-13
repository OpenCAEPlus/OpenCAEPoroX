/*! \file    OCPFlux.hpp
 *  \brief   OCPFlux 类的声明
 *  \author  Shizhe Li
 *  \date    May/10/2023
 *
 *  OCPFlux 类用于计算多相流体流动模拟中的通量，包括组分和相的流量计算，
 *  以及用于不同数值方法的矩阵组装。
 *
 *---------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------
 */

#ifndef __OCPFLUX_HEADER__
#define __OCPFLUX_HEADER__

// 标准头文件
#include <vector>

// OpenCAEPoroX 头文件
#include "Bulk.hpp"
#include "BulkConnVarSet.hpp"
#include "BulkConnTrans.hpp"

using namespace std;

/*!
 * \class OCPFlux
 * \brief 通量计算的抽象基类
 *
 * OCPFlux 类是一个抽象基类，用于定义通量计算的接口。它包含了一些基本的属性和
 * 用于计算和组装通量的虚函数。
 */
class OCPFlux
{
public:
    /*!
     * \brief 默认构造函数
     */
    OCPFlux() = default;

    /*!
     * \brief 计算组分和相的通量
     * \param bp 连接对
     * \param bk 批量对象
     */
    virtual void CalFlux(const BulkConnPair& bp, const Bulk& bk) = 0;

    /*!
     * \brief 为 FIM 方法组装矩阵
     * \param bp 连接对
     * \param c 组件索引
     * \param bcvs 批量连接变量集
     * \param bk 批量对象
     */
    virtual void AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) = 0;

    /*!
     * \brief 为 AIM 方法组装矩阵
     * \param bp 连接对
     * \param c 组件索引
     * \param bcvs 批量连接变量集
     * \param bk 批量对象
     */
    virtual void AssembleMatAIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) = 0;

    /*!
     * \brief 为 IMPEC 方法组装矩阵
     * \param bp 连接对
     * \param c 组件索引
     * \param bcvs 批量连接变量集
     * \param bk 批量对象
     */
    virtual void AssembleMatIMPEC(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) = 0;

    // Getters
    const vector<OCP_USI>& GetUpblock() const { return upblock; }
    const vector<OCP_DBL>& GetRho() const { return rho; }
    const vector<OCP_DBL>& GetFluxVj() const { return flux_vj; }
    const vector<OCP_DBL>& GetFluxNi() const { return flux_ni; }
    const vector<OCP_DBL>& GetFluxHj() const { return flux_Hj; }
    OCP_DBL GetConductH() const { return conduct_H; }
    const vector<OCP_DBL>& GetdFdXpB() const { return dFdXpB; }
    const vector<OCP_DBL>& GetdFdXpE() const { return dFdXpE; }
    const vector<OCP_DBL>& GetdFdXsB() const { return dFdXsB; }
    const vector<OCP_DBL>& GetdFdXsE() const { return dFdXsE; }
    OCP_DBL GetValbb() const { return valbb; }
    OCP_DBL GetValee() const { return valee; }
    OCP_DBL GetRhsb() const { return rhsb; }
    OCP_DBL GetRhse() const { return rhse; }

protected:
    /*!
     * \brief 设置通量类的基本属性
     * \param npin 相的数量
     * \param ncin 组件的数量
     */
    void Setup(const USI& npin, const USI& ncin) {
        Allocate(npin, ncin);
        bcT.Setup();
    }

    /*!
     * \brief 分配内存空间
     * \param npin 相的数量
     * \param ncin 组件的数量
     */
    void Allocate(const USI& npin, const USI& ncin) {
        np = npin;
        nc = ncin;
        upblock.resize(np, 0.0);
        rho.resize(np, 0.0);
        flux_vj.resize(np, 0.0);
        flux_ni.resize(nc, 0.0);
    }

protected:
    /// 相的数量
    USI              np;
    /// 组件的数量
    USI              nc;
    /// 每个相的上游批量索引
    vector<OCP_USI>  upblock;
    /// 连接处相的质量密度
    vector<OCP_DBL>  rho;
    /// 从上游批量流出的相的体积流量
    vector<OCP_DBL>  flux_vj;
    /// 组件的摩尔流量
    vector<OCP_DBL>  flux_ni;
    /// 从上游批量流出的相的焓流量
    vector<OCP_DBL>  flux_Hj;
    // 热传导项
    OCP_DBL          conduct_H;
    // FIM 方法使用
    /// bId 批量的 dF / dXp
    vector<OCP_DBL>  dFdXpB;
    /// eId 批量的 dF / dXp
    vector<OCP_DBL>  dFdXpE;
    /// bId 批量的 dF / dXs
    vector<OCP_DBL>  dFdXsB;
    /// eId 批量的 dF / dXs
    vector<OCP_DBL>  dFdXsE;
    // IMPEC 方法使用
    /// b-b 的 val，b-e 的 -val
    OCP_DBL          valbb;
    /// e-e 的 val，e-b 的 -val
    OCP_DBL          valee;
    /// b 中的 rhs
    OCP_DBL          rhsb;
    /// e 中的 rhs
    OCP_DBL          rhse;

public:
    /*!
     * \brief 计算批量连接区域
     * \param bp 连接对
     * \param bk 批量对象
     */
    void CalBulkConnArea(BulkConnPair& bp, const Bulk& bk) { bcT.CalTrans(bp, bk); }

protected:
    /// 批量连接区域的计算
    BulkConnTrans     bcT;
};

/*!
 * \class OCPFlux01
 * \brief 用于等温达西流量计算的类
 *
 * OCPFlux01 类继承自 OCPFlux 类，专门用于计算等温达西流量。
 */
class OCPFlux01 : public OCPFlux
{
public:
    /*!
     * \brief 默认构造函数
     */
    OCPFlux01() = default;

    /*!
     * \brief 构造函数，初始化类属性
     * \param npin 相的数量
     * \param ncin 组件的数量
     */
    OCPFlux01(const USI& npin, const USI& ncin) {
        Setup(npin, ncin);
        dFdXpB.resize((nc + 1) * (nc + 1));
        dFdXpE.resize((nc + 1) * (nc + 1));
        dFdXsB.resize((nc + 1) * (nc + 1) * np);
        dFdXsE.resize((nc + 1) * (nc + 1) * np);
    }

    // 覆盖基类的纯虚函数
    void CalFlux(const BulkConnPair& bp, const Bulk& bk) override;
    void AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override;
    void AssembleMatAIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override;
    void AssembleMatIMPEC(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override;
};

/*!
 * \class OCPFlux02
 * \brief 用于等温达西流量和标准重力排水模型计算的类
 *
 * OCPFlux02 类继承自 OCPFlux 类，除了计算等温达西流量，还包括标准重力排水模型的计算。
 */
class OCPFlux02 : public OCPFlux
{
public:
    /*!
     * \brief 默认构造函数
     */
    OCPFlux02() = default;

    /*!
     * \brief 构造函数，初始化类属性
     * \param npin 相的数量
     * \param ncin 组件的数量
     */
    OCPFlux02(const USI& npin, const USI& ncin) {
        Setup(npin, ncin);
        dFdXpB.resize((nc + 1) * (nc + 1));
        dFdXpE.resize((nc + 1) * (nc + 1));
        dFdXsB.resize((nc + 1) * (nc + 1) * np);
        dFdXsE.resize((nc + 1) * (nc + 1) * np);
    }

    // 覆盖基类的纯虚函数
    void CalFlux(const BulkConnPair& bp, const Bulk& bk) override;
    void AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override;
    void AssembleMatAIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override { OCP_ABORT("NOT USED!"); }
    void AssembleMatIMPEC(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override { OCP_ABORT("NOT USED!"); }
};

/*!
 * \class OCPFluxT01
 * \brief 用于热达西流量计算的类
 *
 * OCPFluxT01 类继承自 OCPFlux 类，专门用于计算热达西流量。
 */
class OCPFluxT01 : public OCPFlux
{
public:
    /*!
     * \brief 默认构造函数
     */
    OCPFluxT01() = default;

    /*!
     * \brief 构造函数，初始化类属性
     * \param npin 相的数量
     * \param ncin 组件的数量
     */
    OCPFluxT01(const USI& npin, const USI& ncin) {
        Setup(npin, ncin);
        flux_Hj.resize(np);
        dFdXpB.resize((nc + 2) * (nc + 2));
        dFdXpE.resize((nc + 2) * (nc + 2));
        dFdXsB.resize((nc + 2) * (nc + 1) * np);
        dFdXsE.resize((nc + 2) * (nc + 1) * np);
    }

    // 覆盖基类的纯虚函数
    void CalFlux(const BulkConnPair & bp, const Bulk& bk) override;
    void AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override;
    void AssembleMatAIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override{}
    void AssembleMatIMPEC(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk) override{}
};

#endif // __OCPFLUX_HEADER__


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           May/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/
