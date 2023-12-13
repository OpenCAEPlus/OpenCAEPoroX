/*! \file    OCPRock.hpp
 *  \brief   本文件包含OCPRock类及其派生类的声明，用于模拟岩石的物理特性。
 *  \author  Shizhe Li
 *  \date    Nov/15/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __ROCK_HEADER__
#define __ROCK_HEADER__

#include <math.h>
// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "BulkVarSet.hpp"

/*!
 *  \class  OCPRock
 *  \brief  OCPRock类是一个抽象类，用于定义岩石物理特性的基本接口。
 */
class OCPRock
{
public:
    OCPRock() = default;

    /*!
     *  \brief  计算岩石的孔隙率以及孔隙率对压力和温度的导数。
     *  \param  P 压力，类型为OCP_DBL。
     *  \param  T 温度，类型为OCP_DBL。
     *  \param  poroInit 初始孔隙率，类型为OCP_DBL。
     *  \param  bcType 体积内容类型，类型为BulkContent。
     */
    virtual void CalPoro(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& poroInit, const BulkContent& bcType) = 0;

    // 获取当前孔隙率
    OCP_DBL GetPoro() const { return poro; }
    // 获取孔隙率对压力的导数
    OCP_DBL GetdPorodP() const { return dPorodP; }
    // 获取孔隙率对温度的导数
    OCP_DBL GetdPorodT() const { return dPorodT; }
    // 获取非孔隙部分的比例
    OCP_DBL Get_Poro() const { return _poro; }
    // 获取非孔隙部分对压力的导数
    OCP_DBL Get_dPorodP() const { return _dPorodP; }
    // 获取非孔隙部分对温度的导数
    OCP_DBL Get_dPorodT() const { return _dPorodT; }
    // 获取单位体积岩石的焓
    OCP_DBL GetHr() const { return Hr; }
    // 获取焓对温度的导数
    OCP_DBL GetdHrdT() const { return dHrdT; }

protected:
    OCP_DBL   poro;      ///< 岩石孔隙率
    OCP_DBL   dPorodP;   ///< 孔隙率对压力的导数
    OCP_DBL   dPorodT;   ///< 孔隙率对温度的导数
    OCP_DBL   _poro;     ///< 非孔隙部分的比例
    OCP_DBL   _dPorodP;  ///< 非孔隙部分对压力的导数
    OCP_DBL   _dPorodT;  ///< 非孔隙部分对温度的导数
    OCP_DBL   Hr;        ///< 单位体积岩石的焓
    OCP_DBL   dHrdT;     ///< 焓对温度的导数
};

/*!
 *  \class  OCPRockIsoT
 *  \brief  OCPRockIsoT类是OCPRock的派生类，用于等温条件下的岩石物理特性模拟。
 */
class OCPRockIsoT : public OCPRock
{
public:
    OCPRockIsoT() = default;

    /*!
     *  \brief  根据岩石参数设置岩石的属性。
     *  \param  param 岩石参数，类型为RockParam。
     */
    void Assign(const RockParam& param) {
        Pref = param.Pref;
        cp1 = param.cp1;
        cp2 = param.cp2;
    }

protected:
    OCP_DBL Pref;   ///< 参考压力
    OCP_DBL cp1;    ///< 第一个系数
    OCP_DBL cp2;    ///< 第二个系数
};

/*!
 *  \class  OCPRockIsoT_Linear
 *  \brief  OCPRockIsoT_Linear类是OCPRockIsoT的派生类，实现了线性孔隙率模型。
 */
class OCPRockIsoT_Linear : public OCPRockIsoT
{
    // poro = poroInit * (1 + phi)
    // poro = poroInit * (1 + cp1 * (P - Pref) + cp2 / 2 * (P - Pref) * (P - Pref))
public:
    OCPRockIsoT_Linear() = default;
    OCPRockIsoT_Linear(const RockParam& param) { Assign(param); };

    void CalPoro(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& poroInit, const BulkContent& bcType) override;
};

/*!
 *  \class  OCPRockT
 *  \brief  OCPRockT类是OCPRock的派生类，用于考虑温度变化对岩石物理特性影响的模型。
 */
class OCPRockT : public OCPRock
{
public:
    OCPRockT() = default;

    /*!
     *  \brief  根据岩石参数设置岩石的属性。
     *  \param  param 岩石参数，类型为RockParam。
     */
    void Assign(const RockParam& param)
    {
        Pref      = param.Pref;
        Tref      = param.Tref;
        cp        = param.cp1;
        ct        = param.ct;
        cpt       = param.cpt;
        ConstRock = param.ConstRock;
        hcp1      = param.HCP1;
        hcp2      = param.HCP2;
    }

    /*!
     *  \brief  计算岩石的焓值。
     *  \param  T 温度，类型为OCP_DBL。
     */
    void CalRockHT(const OCP_DBL& T);

protected:
    OCP_DBL  Pref;      ///< 初始孔隙率时的参考压力
    OCP_DBL  Tref;      ///< 初始孔隙率时的参考温度
    OCP_DBL  cp;        ///< 岩石压缩系数
    OCP_DBL  ct;        ///< 岩石膨胀系数，仅在考虑热效应时使用
    OCP_DBL  cpt;       ///< 交叉项，仅在考虑热效应时使用
    OCP_BOOL ConstRock; ///< 如果为true，岩石体积保持不变，否则，整体体积保持不变
    OCP_DBL hcp1;       ///< 岩石焓值公式的系数，单位Btu/ft^3 - F
    OCP_DBL hcp2;       ///< 岩石焓值公式的系数，单位Btu/ft^3 - F
};

/*!
 *  \class  OCPRockT_Linear
 *  \brief  OCPRockT_Linear类是OCPRockT的派生类，实现了考虑温度变化的线性孔隙率模型。
 */
class OCPRockT_Linear : public OCPRockT
{
    // poro = poroInit * (1 + phi)
    // poro = poroInit * (1 + cp*(P-Pref) - ct*(T-Tref) + cpt*(P-Pref)*(T-Tref))
public:
    OCPRockT_Linear() = default;
    OCPRockT_Linear(const RockParam& param) { Assign(param); };

    void CalPoro(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& poroInit, const BulkContent& bcType) override;
};

/*!
 *  \class  OCPRockT_Exp
 *  \brief  OCPRockT_Exp类是OCPRockT的派生类，实现了考虑温度变化的指数孔隙率模型。
 */
class OCPRockT_Exp : public OCPRockT
{
    // poro = poroInit * (1 + exp(phi))
    // poro = poroInit * exp(cp*(P-Pref) - ct*(T-Tref) + cpt*(P-Pref)*(T-Tref))
public:
    OCPRockT_Exp() = default;
    OCPRockT_Exp(const RockParam& param) { Assign(param); };

    void CalPoro(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& poroInit, const BulkContent& bcType) override;
};

#endif /* end if __ROCK_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/15/2022      Create file                          */
/*----------------------------------------------------------------------------*/
