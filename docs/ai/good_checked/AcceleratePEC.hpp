/*! \file    AcceleratePEC.hpp *  
\brief   AcceleratePEC类声明（AcceleratePEC class declaration）
\author  Shizhe Li *  
\date    Dec/25/2022 * 

-----------------------------------------------------------------------------------
版权所有 (C) 2021--现在 OpenCAEPoroX团队. 版权所有.
根据GNU Lesser General Public License 3.0或更高版本的条款发布.
-----------------------------------------------------------------------------------
*/

#ifndef __ACCELERATEPEC_HEADER__
#define __ACCELERATEPEC_HEADER__

#include "OCPMixtureMethodComp.hpp"
#include <vector>

using namespace std;

// SkipPSAVarset类
class SkipPSAVarset{
public:      
    // 设置nb的值
    void SetNb(const USI& nbin) { nb = nbin; }    
    // 重置SkipPSA变量到上一时间步
    void ResetToLastTimeStep();    
    // 更新SkipPSA变量在上一时间步
    void UpdateLastTimeStep();

public:    
    OCP_BOOL         ifSetup{ OCP_FALSE };  ///< 只需要一个设置（Only one setup is needed）
    OCP_USI          nb;                    ///< 块数（Num of bulk num）
    USI              np;                    ///< 平衡计算中使用的相数（Num of phase used in phase equilibrium calculation）
    USI              nc;                    ///< 平衡计算中使用的组分数（Num of components used in phase equilibrium calculation）
    vector<OCP_BOOL> flag;                  ///< 如果为真，将测试跳过（If true, skip will be test）
    vector<OCP_DBL>  minEigen;              ///< 用于测试跳过的最小特征值（minimum eigenvalue used for testing skipping）
    vector<OCP_DBL>  P;                     ///< 上一步的压力（Pressure at last step）
    vector<OCP_DBL>  T;                     ///< 上一步的温度（Temperature at last step）
    vector<OCP_DBL>  zi;                    ///< 上一步的组分摩尔分数（Mole fraction of components(for test) at last step）
    vector<OCP_BOOL> lflag;                 ///< 上一次的标志（Last flag）
    vector<OCP_DBL>  lminEigen;             ///< 上一次的最小特征值（Last min eigenvalue）
    vector<OCP_DBL>  lP;                    ///< 上一次的压力（Last P）
    vector<OCP_DBL>  lT;                    ///< 上一次的温度（Last T）
    vector<OCP_DBL>  lzi;                   ///< 上一次的组分摩尔分数（Last zi）
};

// SkipPSAMethod类
class SkipPSAMethod{
public:    
    SkipPSAMethod() = default;    
    // 没有预测饱和度计算ftype
    virtual USI CalFtype01(const OCP_USI& bId, const SkipPSAVarset& svs, const OCPMixtureVarSet& mvs) = 0;    
    // 有预测饱和度计算ftype
    virtual USI CalFtype02(const OCP_USI& bId, const SkipPSAVarset& svs, const OCPMixtureVarSet& mvs) = 0;    
    // 计算下一步的跳过信息
    virtual void CalSkipForNextStep(const OCP_USI& bId, SkipPSAVarset& svs, const OCPMixtureVarSet& mvs) = 0;
};

// SkipPSAMethod01类
class SkipPSAMethod01 : public SkipPSAMethod{
public:    
    SkipPSAMethod01(SkipPSAVarset& svs, const OCPMixtureMethodComp* compMin);    
    // 没有预测饱和度计算ftype
    USI CalFtype01(const OCP_USI& bId, const SkipPSAVarset& svs, const OCPMixtureVarSet& mvs) override;    
    // 有预测饱和度计算ftype
    USI CalFtype02(const OCP_USI& bId, const SkipPSAVarset& svs, const OCPMixtureVarSet& mvs) override;    
    // 计算下一步的指示器
    void CalSkipForNextStep(const OCP_USI& bId, SkipPSAVarset& svs, const OCPMixtureVarSet& mvs) override;

protected:    
    OCP_BOOL IfSkip(const OCP_USI& bId, const SkipPSAVarset& svs, const OCPMixtureVarSet& mvs) const;

protected:      
    /// d ln phi[i][j] / d n[k][j]
    vector<OCP_DBL>             lnphiN;    
    /// 用于跳过稳定性分析的矩阵
    vector<OCP_SIN>             skipMatSTA;    
    /// 跳过稳定性分析的矩阵的特征值。只使用最小特征值
    vector<OCP_SIN>             eigenSkip;    
    /// 用于计算特征值的工作空间
    vector<OCP_SIN>             eigenWork;    
    /// 支持模块
    const OCPMixtureMethodComp* compM;
};

// SkipPSA类
class SkipPSA{
public:    
    // 设置SkipPSA
    USI Setup(const OCP_USI& nb, const OCPMixtureMethodComp* compsin);    
    // 设置ifUse为真或假
    void SetUseSkip(const OCP_BOOL& flag) { ifUse = flag; }    
    // 返回ifUse
    OCP_BOOL IfUseSkip() const { return ifUse; }    
    // 没有预测饱和度计算ftype
    USI CalFtype01(const OCP_USI& bId, const USI& mIndex, const OCPMixtureVarSet& mvs)    
    {        
        if (ifUse)  return sm[mIndex]->CalFtype01(bId, vs, mvs);        
        else        return 0;    
    }    
    // 有预测饱和度计算ftype
    USI CalFtype02(const OCP_USI& bId, const USI& mIndex, const OCPMixtureVarSet& mvs)    
    {        
        if (ifUse)  return sm[mIndex]->CalFtype02(bId, vs, mvs);        
        else        return 0;    
    }    
    // 计算下一步的跳过信息
    void CalSkipForNextStep(const OCP_USI& bId, const USI& mIndex, const OCPMixtureVarSet& mvs)    
    {        
        if (ifUse)             sm[mIndex]->CalSkipForNextStep(bId, vs, mvs);    
    }    
    // 重置SkipPSA变量到上一时间步
    void ResetToLastTimeStep() { if(ifUse) vs.ResetToLastTimeStep(); }    
    // 更新SkipPSA变量在上一时间步
    void UpdateLastTimeStep() { if (ifUse) vs.UpdateLastTimeStep(); }

protected:    
    /// 如果使用跳过PSA选项
    OCP_BOOL                 ifUse{ OCP_TRUE };    
    /// 变量集
    SkipPSAVarset            vs;    
    /// 跳过方法
    vector<SkipPSAMethod*>   sm;
};

#endif /* end if __ACCELERATEPEC_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/05/2023      Create file                          */
/*----------------------------------------------------------------------------*/