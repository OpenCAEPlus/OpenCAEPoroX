/*! \file    OCPMixture.hpp 
 *  \brief   OCPMixture类声明
 *  \author  Shizhe Li 
 *  \date    Jul/12/2023 
 * 
 *-----------------------------------------------------------------------------------
 *  版权所有 (C) 2021--现在 OpenCAEPoroX团队. 保留所有权利.
 *  根据GNU较小通用公共许可证3.0或更高版本的条款发布.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPMIXTURE_HEADER__
#define __OCPMIXTURE_HEADER__

#include "OCPMixtureMethodComp.hpp"
#include "OCPMixtureMethodK.hpp"
#include "OptionalModules.hpp"

using namespace std;

/*! \class OCPMixture
 *  \brief OCPMixture类
 * 
 *  OCPMixture类是一个基类，包含了一些虚函数以及一些基本的成员变量和函数
 */
class OCPMixture{
public:    
    /// 默认构造函数
    OCPMixture() = default;    

    /// 返回混合类型
    auto MixtureType() const { return vs.mixtureType; }    

    /// 获取变量集
    const OCPMixtureVarSet& GetVarSet() const { return vs; }    

    /// 虚函数Flash
    virtual void Flash(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) = 0;    

    /// 虚函数InitFlash
    virtual void InitFlash(const OCP_USI& bId, const BulkVarSet& bvs) = 0;    

    /// 虚函数Flash
    virtual void Flash(const OCP_USI& bId, const BulkVarSet& bvs) = 0;    

    /// 虚函数InitFlashDer
    virtual void InitFlashDer(const OCP_USI& bId, const BulkVarSet& bvs) = 0;    

    /// 虚函数FlashDer
    virtual void FlashDer(const OCP_USI& bId, const BulkVarSet& bvs) = 0;    

    /// 虚函数CalVStd
    virtual void CalVStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni) = 0;    

    /// 虚函数CalVmStd
    virtual OCP_DBL CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;    

    /// 虚函数CalXi
    virtual OCP_DBL CalXi(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;    

    /// 虚函数CalRho
    virtual OCP_DBL CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt) = 0;    

    /// 虚函数OutputIters
    virtual void OutputIters() const = 0;    

    /// 虚函数CalEnthalpy
    virtual OCP_DBL CalEnthalpy(const OCP_DBL& T, const OCP_DBL* zi) const = 0;    

    /// 虚函数IfWellFriend
    virtual OCP_BOOL IfWellFriend() const = 0;

public:    
    /// 获取油指数
    auto OilIndex() const { return vs.o; }    

    /// 获取气指数
    auto GasIndex() const { return vs.g; }    

    /// 获取水指数
    auto WatIndex() const { return vs.w; }    

    /// 获取液体指数
    auto LiquidIndex() const { return vs.l; }

protected:    
    /// 混合变量集
    OCPMixtureVarSet vs;
};

#endif /* end if __OCPMIXTURE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/12/2023      Create file                          */
/*----------------------------------------------------------------------------*/