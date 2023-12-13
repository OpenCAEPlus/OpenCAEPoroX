```cpp
/*! 
 * \file    OCPMixtureVarSet.hpp
 * \brief   OCPMixtureVarSet类的声明
 * \author  Shizhe Li
 * \date    Oct/02/2023
 *
 * OCPMixtureVarSet类用于表示和初始化流体混合物的各种属性和状态。
 * 它包括了相的类型、相的存在性、各相的物性数据以及它们的导数等。
 * 本文件还定义了相关的枚举类型以表示不同的混合物类型和相类型。
 *
 *-----------------------------------------------------------------------------------
 * Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 * Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPMIXTUREVARSET_HEADER__
#define __OCPMIXTUREVARSET_HEADER__

#include "OCPConst.hpp"
#include <vector>
using namespace std;

/*!
 * \enum    OCPMixtureType
 * \brief   混合物类型
 * 
 * 定义了不同的流体混合物类型，包括单相、黑油模型和组分模型等。
 */
enum class OCPMixtureType : USI {
    SP,         ///< 单相
    BO_OG,      ///< 黑油模型：油气
    BO_OW,      ///< 黑油模型：油水
    BO_GW,      ///< 黑油模型：气水
    BO_OGW,     ///< 黑油模型：油气水
    COMP,       ///< 组分模型，恒温
    COMPT,      ///< 组分模型，变温
    THERMALK_OW ///< 热模型：油水
};

/*!
 * \enum    PhaseType
 * \brief   相类型
 * 
 * 定义了流体相的类型，包括油相、气相和水相。
 */
enum class PhaseType : USI {
    oil,    ///< 油相
    gas,    ///< 气相
    wat     ///< 水相
};

/*!
 * \class   OCPMixtureVarSet
 * \brief   表示流体混合物的类
 *
 * OCPMixtureVarSet类用于存储和处理流体混合物的状态和属性。
 * 它包含了初始化混合物属性的方法、计算总流体体积和相饱和度的方法，
 * 以及获取相组分摩尔分数的方法。
 */
class OCPMixtureVarSet {
public:
    /*!
     * \brief   默认构造函数
     */
    OCPMixtureVarSet() = default;

    /*!
     * \brief   初始化混合物属性
     * \param   mixType 混合物类型
     * \param   numPhase 相的数量
     * \param   numCom 组件的数量
     */
    void Init(const OCPMixtureType& mixType, const USI& numPhase, const USI& numCom);

public:
    /*!
     * \brief   计算总流体体积和相饱和度
     */
    void CalVfS();

public:
    /*!
     * \brief   获取相j的组分摩尔分数
     * \param   j 相的索引
     * \return  返回指向相j组分摩尔分数数组的指针
     */
    const OCP_DBL* GetXj(const USI& j) const;

public:
    OCPMixtureType          mixtureType; ///< 混合物类型
    INT                     o, g, w;     ///< 油、气、水相的索引
    vector<INT>             l;           ///< 液相列表
    USI                     np, nc;      ///< 相和组分的数量
    OCP_DBL                 P, T;        ///< 压力、温度
    vector<OCP_DBL>         Pj;          ///< 相压力
    OCP_DBL                 Pb;          ///< 饱和压力
    OCP_DBL                 Nt;          ///< 组件总摩尔数
    OCP_DBL                 Vf;          ///< 组件总体积
    vector<OCP_DBL>         Ni;          ///< 组件摩尔数（某些条件下为质量）
    USI                     phaseNum;    ///< 存在的相数量
    vector<OCP_BOOL>        phaseExist;  ///< 相的存在性
    vector<OCP_DBL>         S;           ///< 相饱和度
    vector<OCP_DBL>         nj;          ///< 相摩尔数
    vector<OCP_DBL>         vj;          ///< 相体积
    vector<OCP_DBL>         x;           ///< 组件i在相j中的摩尔分数
    vector<OCP_DBL>         rho;         ///< 相的质量密度
    vector<OCP_DBL>         xi;          ///< 相的摩尔密度（某些条件下为质量密度）
    vector<OCP_DBL>         mu;          ///< 相的粘度
    // 导数（全导数）
    OCP_DBL                 vfP;         ///< dVf / dP
    OCP_DBL                 vfT;         ///< dVf / dT
    vector<OCP_DBL>         vfi;         ///< dVf / dNi
    vector<OCP_DBL>         vjP;         ///< dVj / dP
    vector<vector<OCP_DBL>> vji;         ///< dVj / dNi
    // 导数（偏导数）
    vector<OCP_DBL>         rhoP;        ///< drho / dP
    vector<OCP_DBL>         rhoT;        ///< drho / dT
    vector<OCP_DBL>         rhox;        ///< drho / dx
    vector<OCP_DBL>         xiP;         ///< dxi / dP
    vector<OCP_DBL>         xiT;         ///< dxi / dT
    vector<OCP_DBL>         xix;         ///< dxi / dx
    vector<OCP_DBL>         muP;         ///< dmu / dP
    vector<OCP_DBL>         muT;         ///< dmu / dT
    vector<OCP_DBL>         mux;         ///< dmu / dx
    OCP_DBL                 Uf;          ///< 单位体积流体的内能
    OCP_DBL                 UfP;         ///< dUf / dP
    OCP_DBL                 UfT;         ///< dUf / dT
    vector<OCP_DBL>         Ufi;         ///< dUf / dNi
    vector<OCP_DBL>         H;           ///< 焓
    vector<OCP_DBL>         HT;          ///< d Hj / d T
    vector<OCP_DBL>         Hx;          ///< d Hj / d xij
    vector<OCP_DBL>         dXsdXp;      ///< d(Sj, xij) / d(P,Ni,(T))
};

#endif /* end if __OCPMIXTURE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/02/2023      Create file                          */
/*----------------------------------------------------------------------------*/
```

