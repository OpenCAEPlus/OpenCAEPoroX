/*! \file    OptionalModules.hpp
 *  \brief   OptionalModules 类声明文件
 *  \author  Shizhe Li
 *  \date    Dec/25/2022
 *
 *  本文件包含 OptionalModules 类的声明，此类负责管理模拟软件中的可选模块。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OPTIONALMODULE_HEADER__
#define __OPTIONALMODULE_HEADER__

#include "AcceleratePEC.hpp"
#include "OCPMiscible.hpp"
#include "OCPScalePcow.hpp"
#include "OCPBoundary.hpp"
#include "HeatConduct.hpp"

/*!
 * \class OptionalModules
 * \brief 管理模拟软件中的可选模块。
 *
 * OptionalModules 类负责初始化和更新模拟过程中的可选模块，如边界条件处理和热传导计算。
 */
class OptionalModules
{
    friend class MixtureUnit;
    friend class OCPMixtureComp;
    friend class FlowUnit;
    friend class FlowUnit_OW;
    friend class Reservoir;
    // For Output
    friend class Out4RPT;

public:
    /*!
     * \brief 重置到上一个时间步的状态。
     */
    void ResetToLastTimeStep()
    {
        skipPSA.ResetToLastTimeStep();
        surTen.ResetTolastTimeStep();
        misFac.ResetTolastTimeStep();
        misCur.ResetTolastTimeStep();
        scalePcow.ResetTolastTimeStep();
        boundary.ResetToLastTimeStep();
        heatConduct.ResetToLastTimeStep();
    }

    /*!
     * \brief 更新到上一个时间步的状态。
     */
    void UpdateLastTimeStep()
    {
        skipPSA.UpdateLastTimeStep();
        surTen.UpdateLastTimeStep();
        misFac.UpdateLastTimeStep();
        misCur.UpdateLastTimeStep();
        scalePcow.UpdateLastTimeStep();
        boundary.UpdateLastTimeStep();
        heatConduct.UpdateLastTimeStep();
    }

    /////////////////////////////////////////////////////////////////////
    // Common vars
    /////////////////////////////////////////////////////////////////////

protected:
    /// \brief 存储流体块数量的变量
    OCP_USI nb;

///////////////////////////////////////////////////////////////////////
// Attachment Module
/////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////
    // Accelerate PVT
    /////////////////////////////////////////////////////////////////////

protected:
    /// \brief 用于跳过稳定性分析的模块
    SkipPSA        skipPSA;

    /////////////////////////////////////////////////////////////////////
    // Phase Permeability Curve
    /////////////////////////////////////////////////////////////////////

protected:
    /// \brief 表面张力计算模块
    SurfaceTension surTen;
    /// \brief 可混合因子计算模块
    MiscibleFactor misFac;
    /// \brief 可混合曲线校正模块
    MiscibleCurve  misCur;
    /// \brief 调整水油毛细压力的模块
    ScalePcow      scalePcow;

///////////////////////////////////////////////////////////////////////
// Independent Module
/////////////////////////////////////////////////////////////////////

public:
    /*!
     * \brief 设置独立模块。
     *
     * \param rs_param 油藏参数
     * \param bvs 流体块变量集
     */
    void SetupIndependentModule(const ParamReservoir& rs_param, const BulkVarSet& bvs)
    {
        boundary.Setup(rs_param, bvs);
        heatConduct.Setup(rs_param, bvs);
    }

public:
    /// \brief 边界条件处理模块
    OCPBoundary    boundary;
    /// \brief 热传导计算模块
    HeatConduct    heatConduct;
};

#endif /* end if __OptionalFeatures_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/25/2022      Create file                          */
/*----------------------------------------------------------------------------*/
