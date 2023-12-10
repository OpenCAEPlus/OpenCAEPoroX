```cpp
/*! \file    OCPMixtureMethodK.hpp
 *  \brief   本文件包含OCPMixtureMethodK类的声明及其派生类
 *  \author  Shizhe Li
 *  \date    Jul/31/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPMIXTUREMETHODK_HEADER__
#define __OCPMIXTUREMETHODK_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPFuncPVT.hpp"
#include "OCPMixtureVarSet.hpp"
#include "OCPPhaseEquilibrium.hpp"
#include "BulkVarSet.hpp"
#include <vector>

using namespace std;

/**
 * \class OCPMixtureMethodK
 * \brief OCPMixtureMethodK是一个在非EoS模型中使用的基础类。
 */
class OCPMixtureMethodK
{
public:
    OCPMixtureMethodK() = default;
    
    /**
     * \brief 设置变量集
     * \param bId 块ID
     * \param bvs BulkVarSet对象的引用
     * \param mvs OCPMixtureVarSet对象的引用
     */
    virtual void SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const = 0;
    
    /**
     * \brief 设置变量集
     * \param P 压力
     * \param T 温度
     * \param Ni 组分摩尔数数组
     * \param mvs OCPMixtureVarSet对象的引用
     */
    virtual void SetVarSet(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, OCPMixtureVarSet& mvs) const = 0;
    
    /**
     * \brief 执行闪蒸计算并计算VfP,Vfi
     * \param Vp 体积百分比
     * \param vs OCPMixtureVarSet对象的引用
     */
    virtual void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    
    /**
     * \brief 执行闪蒸计算并计算VfP,Vfi
     * \param vs OCPMixtureVarSet对象的引用
     */
    virtual void Flash(OCPMixtureVarSet& vs) = 0;
    
    /**
     * \brief 执行闪蒸计算并计算VfP,Vfi,dXsdXp
     * \param Vp 体积百分比
     * \param vs OCPMixtureVarSet对象的引用
     */
    virtual void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) = 0;
    
    /**
     * \brief 执行闪蒸计算并计算VfP,Vfi,dXsdXp
     * \param vs OCPMixtureVarSet对象的引用
     */
    virtual void FlashDer(OCPMixtureVarSet& vs) = 0;
    
    /**
     * \brief 在标准条件下进行闪蒸计算
     * \param vs OCPMixtureVarSet对象的引用
     */
    virtual void CalVStd(OCPMixtureVarSet& vs) = 0;
    
    /**
     * \brief 计算目标相的摩尔密度
     * \param P 压力
     * \param Pb 泡点压力
     * \param T 温度
     * \param z 组分摩尔分数数组
     * \param pt 相类型
     * \return 摩尔密度
     */
    virtual OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;
    
    /**
     * \brief 在标准条件下计算目标相的摩尔密度
     * \param P 压力
     * \param Pb 泡点压力
     * \param T 温度
     * \param z 组分摩尔分数数组
     * \param pt 相类型
     * \return 摩尔密度
     */
    virtual OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;
    
    /**
     * \brief 计算目标相的质量密度
     * \param P 压力
     * \param Pb 泡点压力
     * \param T 温度
     * \param z 组分摩尔分数数组
     * \param pt 相类型
     * \return 质量密度
     */
    virtual OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;
    
    /**
     * \brief 输出迭代次数
     */
    virtual void OutIters() const { }
    
    /**
     * \brief 判断当前方法是否适用于井
     * \return 是否适用于井的布尔值
     */
    virtual OCP_BOOL IfWellFriend() const = 0;
    
    /**
     * \brief 计算焓值
     * \param T 温度
     * \param zi 组分摩尔分数数组
     * \return 焓值
     */
    virtual OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const = 0;
};


/// 使用PVDO表格和PVTW表格计算油水两相混合性质
class OCPMixtureMethodK_OW01 : public OCPMixtureMethodK
{
public:
    /// 构造函数
    OCPMixtureMethodK_OW01(const ParamReservoir& rs_param, const USI& i, OCPMixtureVarSet& vs);
    /// 设置变量集
    void SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const override;
    /// 设置变量集
    void SetVarSet(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, OCPMixtureVarSet& mvs) const override;
    /// 用于初始化的Flash计算
    void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    /// 模拟过程中的Flash计算
    void Flash(OCPMixtureVarSet& vs) override;
    /// 用于初始化的Flash计算，计算大量导数
    void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    /// 模拟过程中的Flash计算，计算大量导数
    void FlashDer(OCPMixtureVarSet& vs) override;
    /// 计算标态下的各相体积
    void CalVStd(OCPMixtureVarSet& vs) override;
    /// 计算指定相的摩尔密度
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    /// 计算指定相的质量密度
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    /// 计算指定相标态下的摩尔体积
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    /// 是否对井友好
    OCP_BOOL IfWellFriend() const override { return OCP_TRUE; }
    /// 计算混合物的焓，此处不可用
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override { OCP_ABORT("Not Used!"); }

protected:
    /// 计算油相摩尔浓度
    OCP_DBL CalXiO(const OCP_DBL& P) { return PVDO->CalXiO(P); }
    /// 计算水相摩尔浓度
    OCP_DBL CalXiW(const OCP_DBL& P) { return PVTW.CalXiW(P); }
    /// 计算油相密度
    OCP_DBL CalRhoO(const OCP_DBL& P) { return PVDO->CalRhoO(P); }
    /// 计算水相密度
    OCP_DBL CalRhoW(const OCP_DBL& P) { return PVTW.CalRhoW(P); }


protected:

    OCP_PVDO*     PVDO;          ///< PVDO表格
    OCP_PVTW      PVTW;          ///< PVTW表格
    const OCP_DBL stdVo{ 1 };    ///< 标态下油相的摩尔体积
    const OCP_DBL stdVw{ 1 };    ///< 标态下水相的摩尔体积
};


/// 计算油水两相热流混合性质
class OCPMixtureMethodK_OW01T : public OCPMixtureMethodK
{
public:
    /// 构造函数
    OCPMixtureMethodK_OW01T(const ComponentParam& param, const USI& tarId, OCPMixtureVarSet& vs);
    /// 设置变量集
    void SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const override;
    /// 设置变量集
    void SetVarSet(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, OCPMixtureVarSet& mvs) const override;
    /// 用于初始化的Flash计算
    void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    /// 模拟过程中的Flash计算
    void Flash(OCPMixtureVarSet& vs) override;
    /// 用于初始化的Flash计算，计算大量导数
    void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    /// 模拟过程中的Flash计算，计算大量导数
    void FlashDer(OCPMixtureVarSet& vs) override;
    /// 计算标态下的各相体积
    void CalVStd(OCPMixtureVarSet& vs) override;
    /// 计算指定相的摩尔密度
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    /// 计算指定相的质量密度
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    /// 计算指定相标态下的摩尔体积
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override { return 1 / CalXi(P, 0, T, 0, pt); }
    /// 是否对井友好
    OCP_BOOL IfWellFriend() const override { return OCP_FALSE; }
    /// 计算混合物的焓
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override { return eC.CalEnthalpy(T + CONV4, zi); }

protected:
    /// 计算油相摩尔浓度
    OCP_DBL CalXiO(const OCP_DBL& P, const OCP_DBL& T);
    /// 计算水相摩尔浓度
    OCP_DBL CalXiW(const OCP_DBL& P, const OCP_DBL& T);
    /// 计算油相密度
    OCP_DBL CalRhoO(const OCP_DBL& P, const OCP_DBL& T);
    /// 计算水相密度
    OCP_DBL CalRhoW(const OCP_DBL& P, const OCP_DBL& T);

protected:
    EnthalpyCalculation  eC;  ///< 焓计算模块
    ViscosityCalculation vC;  ///< 粘度计算模块

protected:
    OCP_DBL Pref;                 ///< 参考温度
    OCP_DBL Tref;                 ///< 参考压力
    vector<OCP_DBL> xi_ref;       ///< 参考摩尔浓度
    vector<OCP_DBL> MWc;          ///< 组分分子质量
    vector<OCP_DBL> MWp;          ///< 相分子质量
                                  
    vector<OCP_DBL> cp;           ///< 组分压缩系数
    vector<OCP_DBL> ct1;          ///< 第一热膨胀系数
    vector<OCP_DBL> ct2;          ///< 第二热膨胀系数
    vector<OCP_DBL> cpt;          ///< 压力温度依赖的密度系数
                                  
    vector<OCP_DBL> avg;          ///< 粘度计算变量
    vector<OCP_DBL> bvg;          ///< 粘度计算变量
};


/////////////////////////////////////////////////////
// OCPMixtureMethodK_OGW01
/////////////////////////////////////////////////////


/// 使用PVDO和PVTW来计算油气水三相混合物的性质
class OCPMixtureMethodK_OGW01 : public OCPMixtureMethodK
{
public:
    /// 构造函数
    OCPMixtureMethodK_OGW01(const ParamReservoir& rs_param, const USI& i, OCPMixtureVarSet& vs);
    /// 设置变量集
    void SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const override;
    /// 设置变量集
    void SetVarSet(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, OCPMixtureVarSet& mvs) const override;
    /// 用于初始化的Flash计算
    void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    /// 模拟过程中的Flash计算
    void Flash(OCPMixtureVarSet& vs) override;
    /// 用于初始化的Flash计算，计算大量导数
    void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    /// 模拟过程中的Flash计算，计算大量导数
    void FlashDer(OCPMixtureVarSet& vs) override;
    /// 计算标态下的各相体积
    void CalVStd(OCPMixtureVarSet& vs) override;
    /// 计算指定相的摩尔密度
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    /// 计算指定相的质量密度
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    /// 计算指定相标态下的摩尔体积
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    /// 是否对井友好
    OCP_BOOL IfWellFriend() const override { return OCP_TRUE; }
    /// 计算混合物的焓
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override { OCP_ABORT("Not Used!"); }

protected:
    /// 计算油相摩尔浓度
    OCP_DBL CalXiO(const OCP_DBL& P, const OCP_DBL& Pb) { return PVCO.CalXiO(P, Pb); }
    /// 计算气相摩尔浓度
    OCP_DBL CalXiG(const OCP_DBL& P) { return PVDG.CalXiG(P); }
    /// 计算水相摩尔浓度
    OCP_DBL CalXiW(const OCP_DBL& P) { return PVTW.CalXiW(P); }
    /// 计算油相密度
    OCP_DBL CalRhoO(const OCP_DBL& P, const OCP_DBL& Pb) { return PVCO.CalRhoO(P, Pb); }
    /// 计算气相密度
    OCP_DBL CalRhoG(const OCP_DBL& P) { return PVDG.CalRhoG(P); }
    /// 计算水相密度
    OCP_DBL CalRhoW(const OCP_DBL& P) { return PVTW.CalRhoW(P); }
    /// 计算组分数
    void CalNi(const OCP_DBL& Vp, OCPMixtureVarSet& vs);

protected:
    OCP_PVCO        PVCO;           ///< PVCO表格
    OCP_PVDG        PVDG;           ///< PVDG表格
    OCP_PVTW        PVTW;           ///< PVTW表格
    OCP_DBL         x;              ///< 1摩尔油组分能吸收的气组分摩尔数
    const OCP_DBL   stdVo{ 1 };     ///< 标态下油相摩尔体积
    const OCP_DBL   stdVg{ 1 };     ///< 标态下气相摩尔体积
    const OCP_DBL   stdVw{ 1 };     ///< 标态下水相摩尔体积
};


/////////////////////////////////////////////////////
// OCPMixtureMethodK_GW01
/////////////////////////////////////////////////////


/// 利用PVTCO2和PVTH2O计算气水两相的混合物性质
class OCPMixtureMethodK_GW01 : public OCPMixtureMethodK
{
public:
    /// 构造函数
    OCPMixtureMethodK_GW01(const ParamReservoir& rs_param, const USI& i, OCPMixtureVarSet& vs);
    /// 设置变量集
    void SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const override;
    /// 设置变量集
    void SetVarSet(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, OCPMixtureVarSet& mvs) const override;
    /// 用于初始化的Flash计算
    void InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    /// 模拟过程中的Flash计算
    void Flash(OCPMixtureVarSet& vs) override;
    /// 用于初始化的Flash计算，计算大量导数
    void InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs) override;
    /// 模拟过程中的Flash计算，计算大量导数
    void FlashDer(OCPMixtureVarSet& vs) override;
    /// 计算标态下的各相体积
    void CalVStd(OCPMixtureVarSet& vs) override;
    /// 计算指定相的摩尔密度
    OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    /// 计算指定相的质量密度
    OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    /// 计算指定相标态下的摩尔体积
    OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) override;
    /// 是否对井友好
    OCP_BOOL IfWellFriend() const override { return OCP_FALSE; }
    /// 计算混合物的焓，此处不可用
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override { OCP_ABORT("Not Used!"); }

protected:
    /// 计算组分数
    void CalNi(const OCP_DBL& Vp, OCPMixtureVarSet& vs);
    /// 计算气相摩尔浓度
    OCP_DBL CalXiG(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z) const { return CalRhoG(P, T, z); }
    /// 计算水相摩尔浓度
    OCP_DBL CalXiW(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z) const { return CalRhoW(P, T, z); }
    /// 计算气相密度
    OCP_DBL CalRhoG(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z) const { return PVTCO2.CalRho(P, T); }
    /// 计算水相密度
    OCP_DBL CalRhoW(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z) const;

protected:
    OCP_PVTCO2    PVTCO2;   ///< PVTCO2表格
    OCP_PVTH2O    PVTH2O;   ///< PVTH2O表格
    Garciaw       garciaw;  ///< 是否使用Garciaw模型校正水相密度
};


#endif /* end if __OCPMIXTUREKMETHOD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/02/2023      Create file                          */
/*----------------------------------------------------------------------------*/

#endif /* end if __OCPMIXTUREKMETHOD_HEADER__ */
```

