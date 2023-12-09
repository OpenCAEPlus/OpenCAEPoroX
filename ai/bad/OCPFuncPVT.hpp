/*! \file    OCPFuncPVT.hpp
*   \brief   Functions for PVT in OCP
*   \author  Shizhe Li
*   \date    Jun/18/2023
*
*-----------------------------------------------------------------------------------
*  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
*  Released under the terms of the GNU Lesser General Public License 3.0 or later.  
*-----------------------------------------------------------------------------------
*/

#ifndef __OCPFUNCPVT_HEADER__ 
#define __OCPFUNCPVT_HEADER__  

// OpenCAEPoroX header files
#include "OCPFuncTable.hpp"  
#include "UtilMath.hpp"
using namespace std;

/** @defgroup PVT 函数
* @brief 用于计算流体的PVT性质
* @details 包含了众多类型的PVT函数，他们负责计算流体的PVT性质如密度，粘度等
* @{
*/

/** @brief 水相性质表格函数(PVTW)
* @details 计算水相的密度，粘度等性质，他们是关于水相压力的函数
*/
class OCP_PVTW : public OCPFuncTable
{
    public:
        //! 默认构造函数
        OCP_PVTW() = default; 

        //! 建立PVT表格
        void Setup(const vector<vector<OCP_DBL>>& src, const OCP_DBL& stdRhoWin, const OCP_DBL& stdVwin);
        
        //! 计算水相的摩尔浓度
        OCP_DBL CalXiW(const OCP_DBL& P) const;  

        //! 计算水相的密度
        OCP_DBL CalRhoW(const OCP_DBL& P) const;

        //! 计算水相的质量密度，摩尔浓度，粘度以及相关的导数性质
        void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu, 
                           OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const;

    protected:
        //! 计算水相的体积变化系数
        OCP_DBL CalBw(const OCP_DBL& P) const;

        //! 计算水相的体积变化系数，粘度等导数
        void CalBwMuwDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu, OCP_DBL& bP, OCP_DBL& muP) const;

    protected:
        //! 标态下水相的密度
        OCP_DBL stdRhoW;
        
        //! 标态下水相的摩尔体积
        OCP_DBL stdVw;  
};

/** @brief 可压缩的活油的PVT性质 
* @details 计算油相的密度，粘度等性质，他们是关于油相压力的函数
*/
class OCP_PVCO : public OCPFuncTable  
{
    public:
        //! 默认构造函数
        OCP_PVCO() = default;

        //! 建立PVCO表格
        void Setup(const vector<vector<OCP_DBL>>& src, const OCP_DBL& stdRhoOin, 
                   const OCP_DBL& stdRhoGin, const OCP_DBL& stdVoin, const OCP_DBL& stdVgin);
        
        //! 计算油相的密度
        OCP_DBL CalRhoO(const OCP_DBL& P, const OCP_DBL& Pb) const;

        //! 计算油相的摩尔浓度
        OCP_DBL CalXiO(const OCP_DBL& P, const OCP_DBL& Pb) const;

        //! 计算饱和油相的密度，摩尔浓度，粘度，气油比及相应的导数
        void CalRhoXiMuRsDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu, 
                             OCP_DBL& rs, OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP, OCP_DBL& rsP) const;
                             
        //! 计算不饱和油相的密度，摩尔浓度，粘度及相应的导数
        void CalRhoXiMuDer(const OCP_DBL& rs, const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, 
                           OCP_DBL& mu, OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP, 
                           OCP_DBL& rhoRs, OCP_DBL& xiRs, OCP_DBL& muRs) const;
                           
        //! 计算饱和油相的气油比
        OCP_DBL CalRs(const OCP_DBL& P) const;

    protected:
        //! 计算饱和油相的气油比，体积变化因子，粘度，及相应的导数 
        void CalRsBoMuoDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& rs, OCP_DBL& mu,
                           OCP_DBL& bP, OCP_DBL& rsP, OCP_DBL& muP) const;
                           
        //! 计算不饱和油相的体积变化因子，粘度，及相应的导数
        void CalBoMuoDer(const OCP_DBL& rs, const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu,  
                         OCP_DBL& bP, OCP_DBL& muP, OCP_DBL& bRs, OCP_DBL& muRs) const;
                         
    protected:
        //! 标态下油相密度
        OCP_DBL stdRhoO;   
        
        //! 标态下气相密度
        OCP_DBL stdRhoG;
        
        //! 标态下油相摩尔体积
        OCP_DBL stdVo;
        
        //! 标态下气相摩尔体积
        OCP_DBL stdVg;              
};

/** @brief 干燥气相的PVT性质
* @details 计算干燥气相的密度，粘度等性质，他们是关于油相压力的函数
*/
class OCP_PVDG : public OCPFuncTable
{
    public:
        //! 默认构造函数
        OCP_PVDG() = default;

        //! 建立PVDG表格函数
        void Setup(const vector<vector<OCP_DBL>>& src, 
                   const OCP_DBL& stdRhoGin, const OCP_DBL& stdVgin);
                   
        //! 计算气相的摩尔浓度
        OCP_DBL CalXiG(const OCP_DBL& P) const;

        //! 计算气相的密度
        OCP_DBL CalRhoG(const OCP_DBL& P) const;

        //! 计算气相的密度，摩尔浓度，粘度，及相应的导数 
        void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, OCP_DBL& mu,
                           OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const;
                           
    protected:
        //! 计算气相的体积变化因子 
        OCP_DBL CalBg(const OCP_DBL& P) const;

        //! 计算气相的体积变化因子，粘度和相应的导数
        void CalBgMugDer(const OCP_DBL& P, OCP_DBL& b, OCP_DBL& mu, 
                         OCP_DBL& bP, OCP_DBL& muP) const;
                         
    protected:
        //! 标态下气相的密度
        OCP_DBL stdRhoG;
        
        //! 标态下气相的摩尔体积
        OCP_DBL stdVg;
};

/** @brief 死油的PVT性质
* @details 计算死油相的密度，粘度等性质，他们是关于油相压力的函数
*/  
class OCP_PVDO : public OCPFuncTable
{
    public:
        //! 默认构造函数
        OCP_PVDO() = default;

        //! 建立PVDO表格函数
        virtual void Setup(const vector<vector<OCP_DBL>>& src, 
                           const OCP_DBL& stdRhoOin, const OCP_DBL& stdVoin);
                           
        //! 计算油相的摩尔浓度
        virtual OCP_DBL CalXiO(const OCP_DBL& P) const;

        //! 计算油相的密度
        virtual OCP_DBL CalRhoO(const OCP_DBL& P) const;

        //! 计算油相的密度，摩尔浓度，粘度，及相应的导数
        virtual void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, 
                                   OCP_DBL& mu, OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const;
                                   
    protected:
        //! 计算油相的体积变化因子  
        virtual OCP_DBL CalBo(const OCP_DBL& P) const;

        //! 计算不饱和油相的体积变化因子，粘度，及相应的导数
        virtual void CalBoMuoDer(const OCP_DBL& P, OCP_DBL& bo, OCP_DBL& muo, 
                                 OCP_DBL& dBodP, OCP_DBL& dMuodP) const;
                                 
    protected:
        //! 标态下油相的密度
        OCP_DBL stdRhoO;    
        
        //! 标态下油相的摩尔体积
        OCP_DBL stdVo;     
};

/** @brief 死油的PVT性质(可压性系数为常数)
* @details 计算死油相的密度，粘度等性质，他们是关于油相压力的函数，可压缩系数为常数
*/
class OCP_PVCDO : public OCP_PVDO
{
    public:
        //! 默认构造函数
        OCP_PVCDO() = default;

        //! 建立PVCDO表格函数
        void Setup(const vector<vector<OCP_DBL>>& src, 
                   const OCP_DBL& stdRhoOin, const OCP_DBL& stdVoin) override;
                   
        //! 计算油相的摩尔浓度
        OCP_DBL CalXiO(const OCP_DBL& P) const override;

        //! 计算油相的密度
        OCP_DBL CalRhoO(const OCP_DBL& P) const override;

        //! 计算油相的密度，摩尔浓度，粘度，及相应的导数
        void CalRhoXiMuDer(const OCP_DBL& P, OCP_DBL& rho, OCP_DBL& xi, 
                           OCP_DBL& mu, OCP_DBL& rhoP, OCP_DBL& xiP, OCP_DBL& muP) const override;
                           
    protected:
        //! 计算油相的体积变化因子
        OCP_DBL CalBo(const OCP_DBL& P) const;

        //! 计算不饱和油相的体积变化因子，粘度，及相应的导数 
        void CalBoMuoDer(const OCP_DBL& P, OCP_DBL& bo, OCP_DBL& muo, 
                         OCP_DBL& dBodP, OCP_DBL& dMuodP) const;
                         
    protected:
        //! 参考压力
        OCP_DBL Pref;
        
        //! 参考压力下的油相体积变化因子
        OCP_DBL Bref;
        
        //! 参考压力下的油相可压缩性系数
        OCP_DBL Cb;
        
        //! 参考压力下的油相粘度
        OCP_DBL muref;
        
        //! 参考压力下的油相粘度系数  
        OCP_DBL Cmu;
};

/** @brief 一般3维PVT表格，变量性质随温度和压力变化 
* @details 计算流体的密度，粘度等性质，他们是关于流体压力和温度的函数
*/
class OCP_PVT2 : public OCPFuncTable2
{
    public:
        //! 计算流体密度
        OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& T) const;

        //! 计算密度，粘度和溶解度
        void CalRhoMuSol(const OCP_DBL& P, const OCP_DBL& T, 
                         OCP_DBL& rho, OCP_DBL& mu, OCP_DBL& sol) const;
                         
        //! 计算密度，粘度和溶解度和相应的导数
        void CalRhoMuSolDer(const OCP_DBL& P, const OCP_DBL& T, 
                            OCP_DBL& rho, OCP_DBL& mu, OCP_DBL& sol,
                            OCP_DBL& rhoP, OCP_DBL& muP, OCP_DBL& solP) const;
};

// Typedefs for specialized PVT tables
typedef OCP_PVT2 OCP_PVTCO2; 
typedef OCP_PVT2 OCP_PVTH2O;

/** @brief Garciaw 模型 
* @details 计算溶解了CO2的水相密度
*/
class Garciaw
{
    public:
        //! 建立Garciaw模型 
        void Setup(const OCP_BOOL& flag);

        //! 是否使用该模型
        auto IfUse() const;

        //! 计算水相密度
        void CalRho(const OCP_DBL& T, const OCP_DBL& xGw, OCP_DBL& rhow) const;

        //! 计算水相密度和相关的导数 
        void CalRhoDer(const OCP_DBL& T, const OCP_DBL& xGw, const OCP_DBL& xGwP, 
                       OCP_DBL& rhow, OCP_DBL& rhowP, OCP_DBL& drhow_dxGw) const;
                       
        //! 计算水相密度和相关的导数 
        void CalRhoDer(const OCP_DBL& T, const OCP_DBL& xGw, const OCP_DBL& xGwP,  
                       const OCP_DBL& xGwT, OCP_DBL& rhow, OCP_DBL& rhowP, 
                       OCP_DBL& rhowT, OCP_DBL& drhow_dxGw) const;
                       
    protected:
        OCP_BOOL ifUse;      ///< 是否使用该模型       
        const OCP_DBL MWCO2; ///< 二氧化碳分子质量 
};

/** @brief 粘度计算参数集  
*/
class ViscosityParams 
{
    public:
        //! 构造函数
        ViscosityParams(const OCP_DBL* Pin, const OCP_DBL* Tin, const OCP_DBL* xin);

        //! 构造函数 
        ViscosityParams(const OCP_DBL* Pin, const OCP_DBL* Tin, 
                        const OCP_DBL* xin, const OCP_DBL* xiin);
                        
        //! 构造函数
        ViscosityParams(const OCP_DBL* Pin, const OCP_DBL* Tin,
                        const OCP_DBL* xin, const OCP_DBL* xiin,  
                        const OCP_DBL* xiPin, const OCP_DBL* xiTin,
                        const OCP_DBL* xixin);
                        
    public:
        const OCP_DBL* P;     ///< 压力
        const OCP_DBL* T;     ///< 温度
        const OCP_DBL* x;     ///< 组分摩尔分数
        const OCP_DBL* xi;    ///< 摩尔浓度
        const OCP_DBL* xiP;   ///< 摩尔浓度对压力的导数
        const OCP_DBL* xiT;   ///< 摩尔浓度对温度的导数
        const OCP_DBL* xix;   ///< 摩尔浓度对摩尔分数的导数
};

/** @brief 粘度计算方法 
*/
class ViscosityMethod 
{
    public:
        //! 默认构造函数
        ViscosityMethod() = default;

        //! 计算粘度
        virtual OCP_DBL CalViscosity(const ViscosityParams& vp) = 0;

        //! 计算粘度和相应的导数 
        virtual OCP_DBL CalViscosity(const ViscosityParams& vp, 
                                     OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) = 0;
                                     
    protected:
        USI nc;                ///< 组分数
        vector<OCP_DBL> muc;   ///< 组分的粘度
        vector<OCP_DBL> mucP;  ///< 组分的粘度对压力的导数
        vector<OCP_DBL> mucT;  ///< 组分的粘度对温度的导数
};

/** @brief 利用组分的压力温度依赖的粘度表计算粘度，线性混合规则
*/
class ViscosityMethod01 : public ViscosityMethod
{
    public:
        //! 构造函数 
        ViscosityMethod01(const Table2& tab);

        //! 计算粘度
        OCP_DBL CalViscosity(const ViscosityParams& vp) override;

        //! 计算粘度和相应的导数
        OCP_DBL CalViscosity(const ViscosityParams& vp,  
                             OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) override;
                             
    protected:
        //! 组分的压力温度依赖的粘度表
        OCPTable2 viscTab;
};

/** @brief 利用组分的粘度关联参数计算粘度，线性混合规则
*/ 
class ViscosityMethod02 : public ViscosityMethod
{
    public:
        //! 构造函数
        ViscosityMethod02(const vector<OCP_DBL>& av, const vector<OCP_DBL>& bv);

        //! 计算粘度 
        OCP_DBL CalViscosity(const ViscosityParams& vp) override;

        //! 计算粘度和相应的导数
        OCP_DBL CalViscosity(const ViscosityParams& vp,  
                             OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) override;
                             
    protected:
        //! 组分的粘度关联参数 
        vector<OCP_DBL> avisc;
        
        //! 组分的粘度关联参数
        vector<OCP_DBL> bvisc;
};

/** @brief 利用Lohrenz-Bray-Clark 粘度计算公式
*/
/// Lohrenz-Bray-Clark formula 
class ViscosityMethod03 : public ViscosityMethod
{
public:
    ViscosityMethod03(const ComponentParam& param, const USI& tarId);
    /// 计算粘度
    OCP_DBL CalViscosity(const ViscosityParams& vp) override;
    /// 计算粘度和相应的导数
    OCP_DBL CalViscosity(const ViscosityParams& vp, OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) override;

protected:
    USI             nc;       ///< 组分数
    vector<OCP_DBL> coef;     ///< LBC系数
    vector<OCP_DBL> Tc;       ///< 组分临界温度
    vector<OCP_DBL> Pc;       ///< 组分临界压力
    vector<OCP_DBL> Vcvis;    ///< 组分粘度临界体积
    vector<OCP_DBL> MWC;      ///< 组分分子质量
    OCP_DBL MW;               ///< 相分子质量
    vector<OCP_DBL> sqrtMWC;  ///< 组分分子质量的平方根
    OCP_DBL xPc, xTc, xVc;    ///< 辅助变量
    vector<OCP_DBL> auxA;     ///< 辅助变量
    vector<OCP_DBL> auxB;     ///< 辅助变量
};


/** @brief 粘度计算接口类
*/
class ViscosityCalculation
{
public:
    /// 默认构造函数
    ViscosityCalculation() = default;
    /// 组装
    void Setup(const ComponentParam& param, const USI& tarId);
    /// 计算粘度
    OCP_DBL CalViscosity(const ViscosityParams& vp) {
        return vM->CalViscosity(vp);
    }
    /// 计算粘度和相应的导数
    OCP_DBL CalViscosity(const ViscosityParams& vp, OCP_DBL& muP, OCP_DBL& muT, OCP_DBL* mux) {
        return vM->CalViscosity(vp, muP, muT, mux);
    }
protected:
    ViscosityMethod* vM;  ///< 粘度计算方法集
};


/** @brief 焓计算方法
*/
class EnthalpyMethod
{
public:
    /// 默认构造函数
    EnthalpyMethod() = default;
    /// 计算焓
    virtual OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const = 0;
    /// 计算焓和相应的导数
    virtual OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& HT, OCP_DBL* Hz) const = 0;
};


/** @brief 液相焓和简单焓计算方法
*/
class EnthalpyMethod01 : public EnthalpyMethod
{
public:
    /// 构造函数
    EnthalpyMethod01(const OCP_DBL& Trefin, const vector<OCP_DBL>& cpl1in, const vector<OCP_DBL>& cpl2in,
        const vector<OCP_DBL>& cpl3in, const vector<OCP_DBL>& cpl4in);
    /// 计算焓
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override;
    /// 计算焓和相应的导数
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& HT, OCP_DBL* Hz) const override;

protected:
    OCP_DBL         Tref;  ///< 参考温度
    USI             nc;    ///< 组分数
    vector<OCP_DBL> cpl1;  ///< 液相焓计算系数
    vector<OCP_DBL> cpl2;  ///< 液相焓计算系数
    vector<OCP_DBL> cpl3;  ///< 液相焓计算系数
    vector<OCP_DBL> cpl4;  ///< 液相焓计算系数
};


/** @brief 气相焓计算方法
*/
class EnthalpyMethod02 : public EnthalpyMethod
{
public:
    /// 构造函数
    EnthalpyMethod02(const OCP_DBL& Trefin, const vector<OCP_DBL>& Tcritin,
        const vector<OCP_DBL>& cpg1in, const vector<OCP_DBL>& cpg2in,
        const vector<OCP_DBL>& cpg3in, const vector<OCP_DBL>& cpg4in,
        const vector<OCP_DBL>& hvaprin, const vector<OCP_DBL>& hvrin,
        const vector<OCP_DBL>& evin);
    /// 计算焓
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const override;
    /// 计算焓和相应的导数
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& HT, OCP_DBL* Hz) const override;

protected:
    OCP_DBL         Tref;    ///< 参考温度
    vector<OCP_DBL> Tcrit;   ///< 组分临界温度
    USI             nc;      ///< 组分数
    vector<OCP_DBL> cpg1;    ///< 气相焓计算系数
    vector<OCP_DBL> cpg2;    ///< 气相焓计算系数
    vector<OCP_DBL> cpg3;    ///< 气相焓计算系数
    vector<OCP_DBL> cpg4;    ///< 气相焓计算系数
    vector<OCP_DBL> hvapr;   ///< 蒸发焓计算系数
    vector<OCP_DBL> hvr;     ///< 蒸发焓计算系数
    vector<OCP_DBL> ev;      ///< 蒸发焓计算系数
};


/** @brief 焓计算接口
*/
class EnthalpyCalculation
{
public:
    /// 组装计算
    void Setup(const ComponentParam& param, const USI& tarId);
    /// 计算焓
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const { return eM->CalEnthalpy(T, zi); }
    /// 计算焓和相应的导数
    OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi, OCP_DBL& HT, OCP_DBL* Hz) const { return eM->CalEnthalpy(T, zi, HT, Hz); }
protected:
    EnthalpyMethod* eM; ///< 焓计算方法集
};