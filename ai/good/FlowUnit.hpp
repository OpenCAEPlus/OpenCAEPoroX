/*! \file    FlowUnit.hpp
 *  \brief   流体单元类的声明
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *  本文件包含了FlowUnit类的声明，该类设计用于处理与饱和度表相关的事宜，
 *  包括相对渗透率和毛细压力的计算。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __FLOWUNIT_HEADER__
#define __FLOWUNIT_HEADER__

#include <math.h>
// OpenCAEPoroX头文件
#include "ParamReservoir.hpp"
#include "OptionalModules.hpp"
#include "OCPFlow.hpp"

/**
 *  \class  FlowUnit
 *  \brief  用于处理与饱和度表相关的计算，如相对渗透率和毛细压力。
 */
class FlowUnit
{
public:
    /**
     *  \brief  构造函数
     *  \param  rs_param  油藏参数
     *  \param  i         单元索引
     *  \param  opts      可选模块
     */
    FlowUnit(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts);
    
    /**
     *  \brief  设置比例尺
     *  \param  bId       块ID
     *  \param  Swinout   饱和度输出
     *  \param  Pcowin    毛细压力输入
     */
    void SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin) const;

    /**
     *  \brief  计算相对渗透率和毛细压力
     *  \param  bId       块ID
     *  \param  S         饱和度数组
     */
    void CalKrPc(const OCP_USI& bId, const OCP_DBL* S) const;

    /**
     *  \brief  使用FIM方法计算相对渗透率和毛细压力
     *  \param  bId       块ID
     *  \param  S         饱和度数组
     */
    void CalKrPcFIM(const OCP_USI& bId, const OCP_DBL* S) const;

    /**
     *  \brief  获取残余水饱和度
     *  \return 残余水饱和度值
     */
    OCP_DBL GetSwco() const;

    /**
     *  \brief  通过水饱和度计算毛细压力
     *  \param  Sw        水饱和度
     *  \return 毛细压力值
     */
    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const;

    /**
     *  \brief  通过毛细压力计算水饱和度
     *  \param  Pcow      毛细压力
     *  \return 水饱和度值
     */
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const;

    /**
     *  \brief  通过气饱和度计算毛细压力
     *  \param  Sg        气饱和度
     *  \return 毛细压力值
     */
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const;

    /**
     *  \brief  通过毛细压力计算气饱和度
     *  \param  Pcgo      毛细压力
     *  \return 气饱和度值
     */
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const;

    /**
     *  \brief  通过毛细压力计算水饱和度
     *  \param  Pcgw      毛细压力
     *  \return 水饱和度值
     */
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const;

    /**
     *  \brief  获取相对渗透率数组
     *  \return 相对渗透率数组的常量引用
     */
    const vector<OCP_DBL>& GetKr() const;

    /**
     *  \brief  获取毛细压力数组
     *  \return 毛细压力数组的常量引用
     */
    const vector<OCP_DBL>& GetPc() const;

    /**
     *  \brief  获取相对渗透率导数数组
     *  \return 相对渗透率导数数组的常量引用
     */
    const vector<OCP_DBL>& GetdKrdS() const;

    /**
     *  \brief  获取毛细压力导数数组
     *  \return 毛细压力导数数组的常量引用
     */
    const vector<OCP_DBL>& GetdPcdS() const;

protected:
    OCPFlow*        flow;           ///< 用于水油毛细压力曲线的缩放

    /// 缩放Pcow
    ScalePcow*      scalePcow;      ///< 用于缩放水油毛细压力曲线的类实例
    USI             spMethodIndex;  ///< 缩放方法的索引

    // 对于可混溶
    MiscibleCurve*  misCurve;       ///< 可混溶曲线校正的类实例
    USI             mcMethodIndex;  ///< 可混溶校正方法的索引
};

#endif /* end if __FLOWUNIT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/05/2022      Format file                          */
/*----------------------------------------------------------------------------*/