/*! 
 *  \file    OCPSurfacTension.hpp
 *  \brief   本文件包含OCPSurfacTension类的声明。
 *  \author  Shizhe Li
 *  \date    Jul/03/2023
 *
 *  本文件定义了与表面张力计算相关的类和方法。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPSURFACETENSION_HEADER__
#define __OCPSURFACETENSION_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPMixtureVarSet.hpp"
#include <vector>

using namespace std;

/*! 
 *  \class   SurTenVarSet
 *  \brief   表面张力变量集合
 *  \details 用于存储和设置表面张力计算过程中的相关变量。
 */
class SurTenVarSet
{
public:
    /*! 
     *  \brief   设置相的数量。
     *  \param   nbin 相的数量。
     *  \details 此方法用于设置系统中存在的相的数量，这是表面张力计算的前提条件。
     */
    void SetNb(const OCP_USI& nbin) { nb = nbin; }

public:
    OCP_USI         nb;          //!< 相的数量，类型为无符号短整型。
    vector<OCP_DBL> surTen;      //!< 表面张力值的向量，存储各相之间的表面张力值，类型为双精度浮点数向量。
};

/*! 
 *  \class   SurTenMethod
 *  \brief   表面张力计算方法的基类
 *  \details 提供了计算表面张力的接口。
 */
class SurTenMethod
{
public:
    /// 默认构造函数
    SurTenMethod() = default;

    /*! 
     *  \brief   根据给定参数计算表面张力。
     *  \param   bId 相的标识。
     *  \param   stvs 表面张力变量集合。
     *  \param   mvs 混合物变量集合。
     *  \details 此方法定义了表面张力计算的接口，具体实现将由派生类完成。
     */
    virtual void CalSurfaceTension(const OCP_USI& bId, SurTenVarSet& stvs, const OCPMixtureVarSet& mvs) const = 0;
};

/*! 
 *  \class   SurTenMethod01
 *  \brief   Macleod - Sugden 相关法
 *  \details 使用Macleod-Sugden公式来计算油气表面张力。
 */
class SurTenMethod01 : public SurTenMethod
{
public:
    /*! 
     *  \brief   构造函数。
     *  \param   parachorin 液滴的参数。
     *  \param   stvs 表面张力变量集合。
     *  \details 此构造函数初始化Macleod-Sugden方法所需的参数，包括液滴参数和表面张力变量集合。
     */
    SurTenMethod01(const vector<OCP_DBL>& parachorin, SurTenVarSet& stvs) {
        parachor = parachorin;
        NC       = parachor.size();
        stvs.surTen.resize(stvs.nb);
    }

    void CalSurfaceTension(const OCP_USI& bId, SurTenVarSet& stvs, const OCPMixtureVarSet& mvs) const override;

protected:
    vector<OCP_DBL> parachor; //!< Macleod-Sugden相关法中使用的液滴参数，类型为双精度浮点数向量。
    USI             NC;       //!< 组件数量，类型为无符号短整型。
};

/*!
 *  \class   SurfaceTension
 *  \brief   表面张力类
 *  \details 用于管理和计算表面张力相关的方法和变量。
 */
class SurfaceTension
{
public:
    /// 默认构造函数
    SurfaceTension() = default;

    /*! 
     *  \brief   设置表面张力计算方法。
     *  \param   rs_param 储层参数。
     *  \param   i 方法索引。
     *  \param   nb 相的数量。
     *  \return  返回无符号短整型，表示设置结果。
     *  \details 不同的混合物使用不同的设置方法，此函数负责选择并设置适当的表面张力计算方法。
     */
    USI Setup(const ParamReservoir& rs_param, const USI& i, const OCP_USI& nb);

    /*! 
     *  \brief   使用指定的方法计算表面张力。
     *  \param   bId 相的标识。
     *  \param   mIndex 方法索引。
     *  \param   mvs 混合物变量集合。
     *  \details 如果设置了使用表面张力计算，则调用指定方法的计算接口进行计算。
     */
    void CalSurfaceTension(const OCP_USI& bId, const USI& mIndex, const OCPMixtureVarSet& mvs) {
        if (ifUse) {
            stMethod[mIndex]->CalSurfaceTension(bId, vs, mvs);
        }
    }

    const auto& GetVS()const { return vs; } //!< 获取表面张力变量集合的引用。
    const auto IfUse() const { return ifUse; } //!< 返回布尔值，表示是否使用表面张力计算。
    const auto GetSurfaceTension(const OCP_USI& bId) const {
        OCP_ASSERT(ifUse, "Surface Tension is not available!");
        return vs.surTen[bId]; //!< 获取指定相的表面张力值。
    }

    void ResetTolastTimeStep() { } //!< 重置至上一个时间步，待实现。
    void UpdateLastTimeStep()  { } //!< 更新至上一个时间步，待实现。

protected:
    OCP_BOOL              ifUse{ OCP_FALSE }; //!< 布尔变量，用于标识是否计算表面张力。
    SurTenVarSet          vs;                //!< 表面张力变量集合。
    vector<SurTenMethod*> stMethod;          //!< 表面张力计算方法的向量。
};

#endif /* end if __OCPSURFACETENSION_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/03/2023      Create file                          */
/*----------------------------------------------------------------------------*/
