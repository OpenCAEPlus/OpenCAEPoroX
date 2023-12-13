/*! \file    OCPPhaseEquilibrium.hpp
 *  \brief   OCPPhaseEquilibrium 类的声明
 *  \author  Shizhe Li
 *  \date    Jul/28/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPPHASEEQUILIBRIUM_HEADER__
#define __OCPPHASEEQUILIBRIUM_HEADER__

#include "OCPConst.hpp"
#include "OCPEoS.hpp"
#include <vector>

using namespace std;

/*! \class PEIterTol
 *  \brief 迭代容忍度参数类
 *
 *  存储最大迭代次数、容忍度、残差和收敛标志等信息。
 */
class PEIterTol
{
public:
    /// 最大迭代次数
    USI      maxIt;
    /// 容忍度
    OCP_DBL  tol;
    /// 容忍度的平方
    OCP_DBL  tol2;
    /// 残差
    OCP_DBL  res;
    /// 收敛标志，如果收敛，conflag = OCP_TRUE
    OCP_BOOL conflag;
    /// 当前迭代次数
    USI curIt;
};

/*! \class SSMparamSTA
 *  \brief 相稳定性分析中SSM（Successive Substitution Method）的参数类
 *
 *  继承自 PEIterTol 类，添加了 K 容忍度、dY 容忍度等参数。
 */
class SSMparamSTA : public PEIterTol
{
public:
    OCP_DBL  Ktol{ 1E-4 }; ///< K 容忍度的平方
    OCP_DBL  dYtol{ 1E-6 };
    OCP_DBL  eYt{ 1E-8 }; ///< 如果 Yt > 1 + eYt，则单相不稳定
    OCP_DBL curSk;
};

/*! \class NRparamSTA
 *  \brief 相稳定性分析中NR（Newton-Raphson）的参数类
 *
 *  继承自 PEIterTol 类，目前没有额外参数。
 */
class NRparamSTA : public PEIterTol { };

/*! \class SSMparamSP
 *  \brief 相分裂计算中SSM的参数类
 *
 *  继承自 PEIterTol 类，目前没有额外参数。
 */
class SSMparamSP : public PEIterTol { };

/*! \class NRparamSP
 *  \brief 相分裂计算中NR的参数类
 *
 *  继承自 PEIterTol 类，目前没有额外参数。
 */
class NRparamSP : public PEIterTol { };

/*! \class RRparam
 *  \brief 解Rachford-Rice方程的参数类
 *
 *  继承自 PEIterTol 类，目前没有额外参数。
 */
class RRparam : public PEIterTol { };

/*! \class FlashCtrl
 *  \brief 闪蒸控制类
 *
 *  存储相稳定性分析和相分裂计算的参数。
 */
class FlashCtrl
{
public:
    SSMparamSTA SSMsta;
    NRparamSTA  NRsta;
    SSMparamSP  SSMsp;
    NRparamSP   NRsp;
    RRparam     RR;
};

/*! \class OCPPhaseEquilibrium
 *  \brief 相平衡计算类
 *
 *  实现相平衡计算的主要功能，包括设置相平衡参数、计算相平衡以及获取相平衡结果。
 */
class OCPPhaseEquilibrium
{
public:
    OCPPhaseEquilibrium() = default;
    /*! \brief 设置相平衡
     *  \param param 组件参数
     *  \param tarId 目标ID
     *  \param eosin 方程状态计算指针
     */
    void Setup(const ComponentParam& param, const USI& tarId, EoSCalculation* eosin);
    /*! \brief 相平衡API
     *  \param Pin 输入压力
     *  \param Tin 输入温度
     *  \param ziin 组件摩尔分数数组
     *  \param ftypein 闪蒸类型
     *  \param lNPin 上一次闪蒸的相数
     *  \param lxin 上一次闪蒸的组分摩尔分数
     *  \param nc 组件数量
     */
    void PhaseEquilibrium(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* ziin,
                          const USI& ftypein, const OCP_DBL& lNPin, const OCP_DBL* lxin,
                           const USI& nc);
    /*! \brief 相平衡API
     *  \param Pin 输入压力
     *  \param Tin 输入温度
     *  \param ziin 组件摩尔分数数组
     */
    void PhaseEquilibrium(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* ziin);
    /*! \brief 获取当前相数
     *  \return 当前相数
     */
    const auto& GetNP() const { return NP; }
    /*! \brief 获取j相的组分摩尔分数
     *  \param j 相索引
     *  \return j相的组分摩尔分数
     */
    const auto& GetX(const USI& j) const { return x[j]; }
    /*! \brief 获取j相的摩尔分数
     *  \param j 相索引
     *  \return j相的摩尔分数
     */
    const auto& GetNu(const USI& j) const { return nu[j]; }
    /*! \brief 返回闪蒸类型
     *  \return 闪蒸类型
     */
    const auto& GetFtype() const { return ftype; }
protected:
    /*! \brief 设置初始值
     *  \param Pin 输入压力
     *  \param Tin 输入温度
     *  \param Niin 组件摩尔数数组
     *  \param ftypein 闪蒸类型
     *  \param lNPin 上一次闪蒸的相数
     *  \param lxin 上一次闪蒸的组分摩尔分数
     *  \param nc 组件数量
     */
    void SetInitalValue(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin,
                         const USI& ftypein, const OCP_DBL& lNPin, const OCP_DBL* lxin,
                        const USI& nc);

///////////////////////////////////////////////// 基础组分属性 ///////////////////////////////////////////////
protected:
    /// 烃类组分数量
    USI             NC;
    /// 烃类组分名称
    vector<string>  Cname;
    /// 烃类组分临界温度
    vector<OCP_DBL> Tc;
    /// 烃类组分临界压力
    vector<OCP_DBL> Pc;
    /// 烃类组分临界体积
    vector<OCP_DBL> Vc;
    /// 烃类组分分子量
    vector<OCP_DBL> MWC;
    /// 烃类组分偏心因子
    vector<OCP_DBL> Acf;

///////////////////////////////////////////////// 状态方程 ///////////////////////////////////////////////
protected:
    // 状态方程计算
    EoSCalculation*  eos;

///////////////////////////////////////////////// 初始变量 ///////////////////////////////////////////////
protected:
     /// 当前压力
    OCP_DBL         P;
    /// 当前温度
    OCP_DBL         T;
    /// 烃类组分摩尔分数
    vector<OCP_DBL> zi;

///////////////////////////////////////////////// 基础变量 ///////////////////////////////////////////////
protected:
    /// 为基础变量分配内存
    void AllocateBasicVars();
    /// x 转换到 n
    void x2n();
    /// 计算分子量
    void CalMW();
    /// 找到最重的相
    USI  FindMWmax();

protected:
    /// 上一次NR步骤中外部迭代的相数
    USI                      lNP;
    /// 当前相数
    USI                      NP;
    /// 相的摩尔分数
    vector<OCP_DBL>          nu;
    /// j相中i组分的摩尔分数
    vector<vector<OCP_DBL>>  x;
    /// j相中i组分的摩尔数
    vector<vector<OCP_DBL>>  n;
    /// 相分裂计算中NR迭代的上一次摩尔数
    vector<vector<OCP_DBL>>  ln;
    /// 闪蒸前的吉布斯能量（非真实值）
    OCP_DBL                  GibbsEnergyB;
    /// 闪蒸后的吉布斯能量（非真实值）
    OCP_DBL                  GibbsEnergyE;
    /// 相的分子量
    vector<OCP_DBL>          MW;
    /// j相中i组分的逸度系数
    vector<vector<OCP_DBL>>  phi;
    /// j相中i组分的逸度
    vector<vector<OCP_DBL>>  fug;

///////////////////////////////////////////////// 方法 ///////////////////////////////////////////////
protected:
    /// 为相平衡分配内存
    void     AllocateMethodVars();

protected:
    // 闪蒸方法控制
    /// 允许的最大烃类相数
    USI             NPmax;
    /// 解决相平衡问题的方法参数
    FlashCtrl       flashCtrl;

protected:
    void     CalKwilson();
    ///////////////////////////////////////////////
    // 相稳定性分析
    ///////////////////////////////////////////////
protected:
    /// 相稳定性分析API
    OCP_BOOL PhaseStable();
    /// 连续替代方法
    OCP_BOOL StableSSM(const USI& Id);
    /// 放松的连续替代方法
    OCP_BOOL StableSSM01(const USI& Id);
    /// NR 方法
    OCP_BOOL StableNR(const USI& Id);
    /// 为StableNR组装雅可比矩阵
    void     AssembleJmatSTA();

protected:
    // 方法变量
    /// 稳定性分析中测试相的索引
    USI                     testPId;
    /// Whilson平衡常数
    vector<vector<OCP_DBL>> Kw;
    /// SSM中平衡常数的近似值
    vector<vector<OCP_DBL>> Ks;
    /// 上一次的Ks
    vector<OCP_DBL>         lKs;
    /// 相稳定性分析中使用的逸度系数
    vector<OCP_DBL>         phiSta;
    /// 相稳定性分析中使用的逸度
    vector<OCP_DBL>         fugSta;
    // SSM稳定性分析
    /// Y[i] / Yt
    vector<OCP_DBL>         Y;
    /// Y的总和
    OCP_DBL                 Yt;
    /// 测试相中的phi[id] * x[id]
    vector<OCP_DBL>         di;
    // NR稳定性分析
    /// TPD方程的残差
    vector<OCP_DBL>         resSTA;
    ///< TPD方程的雅可比矩阵
    vector<OCP_DBL>         JmatSTA;
    /// d ln fij / d xkj, 在每个子向量中，按k排序。
    vector<vector<OCP_DBL>> lnfugX;
    ///////////////////////////////////////////////
    // 相分裂计算
    ///////////////////////////////////////////////
    /// 相分裂计算API
    void     PhaseSplit();
    /// 连续替代方法API
    void     SplitSSM(const OCP_BOOL& flag);
    /// 两相的连续替代方法
    void     SplitSSM2(const OCP_BOOL& flag);
    /// 大于等于3相的连续替代方法
    void     SplitSSM3(const OCP_BOOL& flag);
    /// 当NP=2时，使用NR解RR方程
    void     RachfordRice2();
    /// 当NP=2时，使用改进但似乎不太稳定的NR解RR方程
    void     RachfordRice2P();
    /// 当NP>=3时，使用NR解RR方程
    void     RachfordRice3();
    /// 使用更新的nu更新x
    void     UpdateXRR();
    /// 使用BFGS解逸度平衡方程
    void     SplitBFGS();
    /// 使用NR解逸度平衡方程
    void     SplitNR();
    /// 计算逸度平衡方程的残差
    void     CalResSP();
    /// 组装逸度平衡方程的雅可比矩阵
    void     AssembleJmatSP();
    /// 计算NR法的步长
    OCP_DBL  CalStepNRsp();
    /// 检查相分裂计算的正确性
    OCP_BOOL CheckSplit();
    // SSM相分裂
    /// Rachford-Rice方程的残差。
    vector<OCP_DBL> resRR;
    // NR相分裂
    /// 逸度平衡方程的残差。
    vector<OCP_DBL>         resSP;
    /// 逸度平衡方程的雅可比矩阵。
    vector<OCP_DBL>         JmatSP;
    /// 上一次的resSP，BFGS中使用（待完成）
    vector<OCP_DBL>         lresSP;
    /// d ln fij / d nkj, 在每个子向量中，按k排序。
    vector<vector<OCP_DBL>> lnfugN;
    // 使用lapack进行线性求解
    /// 在lapack的dgesv_中使用
    vector<OCP_INT>         pivot;
    /// JmatSTA和SP的工作空间
    vector<OCP_DBL>         JmatWork;
    /// JmatWork的长度
    OCP_INT                 lenJmatWork;
    /// 使用dsysv_时存储在上三角中
    char                    uplo{ 'U' };

///////////////////////////////////////////////// 统计 ///////////////////////////////////////////////
public:
    /// 打印PE中总共使用的迭代次数
    void OutMixtureIters() const;

protected:
    /// SSMSTA的总迭代次数
    OCP_ULL itersSSMSTA{ 0 };
    /// NRSTA的总迭代次数
    OCP_ULL itersNRSTA{ 0 };
    /// SSMSP的总迭代次数
    OCP_ULL itersSSMSP{ 0 };
    /// NRSP的总迭代次数
    OCP_ULL itersNRSP{ 0 };
    /// RR的总迭代次数
    OCP_ULL itersRR{ 0 };
    /// 实施SSMSTA的总次数
    OCP_ULL countsSSMSTA{ 0 };
    /// 实施NRSTA的总次数
    OCP_ULL countsNRSTA{ 0 };
    /// 实施SSMSP的总次数
    OCP_ULL countsSSMSP{ 0 };
    /// 实施NRSP的总次数
    OCP_ULL countsNRSP{ 0 };

///////////////////////////////////////////////// 其他变量 ///////////////////////////////////////////////
protected:
    /// 当前闪蒸和下一次闪蒸的起点
    // 输入:
    // ftype == 0: 从单相闪蒸
    // ftype == 1: 单相（跳过单相稳定性分析）
    // ftype == 2: 两相（跳过单相稳定性分析）
    // 输出:
    // ftype == 0: 尝试在下一次闪蒸中跳过当前批次的单相稳定性分析，
    //             并重新计算判断范围
    // ftype == 1: 尝试在下一次闪蒸中跳过当前批次的单相稳定性分析，
    //             但不重新计算判断范围
    // ftype == 2: 在下一次闪蒸中不跳过当前批次的单相稳定性分析。
    USI         ftype{ 0 };
};

#endif /* end if __OCPPHASEEQUILIBRIUM_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/28/2023      Create file                          */
/*----------------------------------------------------------------------------*/
