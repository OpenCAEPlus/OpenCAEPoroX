/*! \file    OCPScalePcow.hpp 
 *  \brief   OCPScalePcow类声明
 *  \author  Shizhe Li 
 *  \date    Jul/01/2023 
 *
 *----------------------------------------------------------------------------------- 
 *  版权所有 (C) 2021--至今 OpenCAEPoroX团队。保留所有权利。
 *  根据GNU Lesser General Public License 3.0或更高版本的条款发布。
 *----------------------------------------------------------------------------------- */

#ifndef __OCPSCALEPCOW__HEADER__
#define __OCPSCALEPCOW__HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPFlow.hpp"
#include <vector>

using namespace std;

///////////////////////////////////////////////////////////////////////
// 缩放水油毛细压力曲线（来自SWATINIT）
/////////////////////////////////////////////////////////////////////

/**
 * \brief 缩放Pcow的变量集合类
 */
class ScalePcowVarSet{
public:
    /**
     * \brief 设置bulks的数量
     * \param nbin bulks的数量
     */
    void SetNb(const OCP_USI& nbin) { nb = nbin; }

public:
    /// bulks的数量
    OCP_USI         nb;
    /// 初始水分布
    vector<OCP_DBL> swatInit;
    /// Pcow的缩放值，将从swatInit计算得到
    vector<OCP_DBL> scaleVal;
};

/**
 * \brief 用于缩放水渗透率曲线的方法类
 */
class ScalePcowMethod{
public:
    /// 默认构造函数
    ScalePcowMethod() = default;

    /**
     * \brief 设置缩放系数
     * \param bId bulks的ID
     * \param spvs 缩放Pcow的变量集合
     * \param Swinout 输入输出饱和度
     * \param Pcowin 输入的Pcow
     */
    virtual void SetScaleVal(const OCP_USI& bId, ScalePcowVarSet& spvs, OCP_DBL& Swinout, const OCP_DBL& Pcowin) const = 0;

    /**
     * \brief 缩放Pcow和dPcowdSw
     * \param bId bulks的ID
     * \param spvs 缩放Pcow的变量集合
     */
    virtual void ScaleDer(const OCP_USI& bId, const ScalePcowVarSet& spvs) const = 0;

    /**
     * \brief 缩放Pcow
     * \param bId bulks的ID
     * \param spvs 缩放Pcow的变量集合
     */
    virtual void Scale(const OCP_USI& bId, const ScalePcowVarSet& spvs) const = 0;
};

/**
 * \brief 用于缩放水渗透率曲线的方法类
 */
class ScalePcowMethod01 : public ScalePcowMethod{
public:
    /**
     * \brief 构造函数，设置Pmax和Pmin
     * \param flowin 流动模型
     * \param spvs 缩放Pcow的变量集合
     */
    ScalePcowMethod01(OCPFlow* flowin, ScalePcowVarSet& spvs);

    /**
     * \brief 设置缩放系数
     * \param bId bulks的ID
     * \param spvs 缩放Pcow的变量集合
     * \param Swinout 输入输出饱和度
     * \param Pcowin 输入的Pcow
     */
    void SetScaleVal(const OCP_USI& bId, ScalePcowVarSet& spvs, OCP_DBL& Swinout, const OCP_DBL& Pcowin) const override;

    /**
     * \brief 缩放Pcow和dPcowdSw
     * \param bId bulks的ID
     * \param spvs 缩放Pcow的变量集合
     */
    void ScaleDer(const OCP_USI& bId, const ScalePcowVarSet& spvs) const override;

    /**
     * \brief 缩放Pcow
     * \param bId bulks的ID
     * \param spvs 缩放Pcow的变量集合
     */
    void Scale(const OCP_USI& bId, const ScalePcowVarSet& spvs) const override;

protected:
    /// 固结水饱和度
    OCP_DBL        Swco;
    /// 最大毛细压力：Po - Pw
    OCP_DBL        maxPcow;
    /// 最小毛细压力：Po - Pw
    OCP_DBL        minPcow;

    // 依赖模块
    /// 流动模型
    OCPFlow*       flow;
};

/**
 * \brief 缩放Pcow的类
 */
class ScalePcow{
    // 输出类
    friend class Out4RPT;

public:
    /**
     * \brief 设置缩放Pcow项
     * \param flowin 流动模型
     * \return 错误码
     */
    USI Setup(OCPFlow* flowin);

    /**
     * \brief 设置缩放系数
     * \param bId bulks的ID
     * \param mIndex 模型的索引
     * \param Swinout 输入输出饱和度
     * \param Pcowin 输入的Pcow
     */
    void SetScaleVal(const OCP_USI& bId, const USI& mIndex, OCP_DBL& Swinout, const OCP_DBL& Pcowin);

    /**
     * \brief 缩放Pcow和dPcowdSw
     * \param bId bulks的ID
     * \param mIndex 模型的索引
     */
    void ScaleDer(const OCP_USI& bId, const USI& mIndex) const;

    /**
     * \brief 缩放Pcow
     * \param bId bulks的ID
     * \param mIndex 模型的索引
     */
    void Scale(const OCP_USI& bId, const USI& mIndex) const;

    /**
     * \brief 获取swatinit
     * \return swatinit
     */
    auto& GetSwatInit() { return vs.swatInit; }

    /**
     * \brief 重置到上一个时间步长
     */
    void ResetTolastTimeStep() { }

    /**
     * \brief 更新上一个时间步长
     */
    void UpdateLastTimeStep() { }

protected:
    /// 是否缩放
    OCP_BOOL ifScale{ OCP_FALSE };

    /// 缩放Pcow的变量集合
    ScalePcowVarSet vs;

    /// 缩放Pcow的方法
    vector<ScalePcowMethod*> scalePcowMethod;
};

#endif /* end if __OCPSCALEPCOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/01/2023      Create file                          */
/*----------------------------------------------------------------------------*/