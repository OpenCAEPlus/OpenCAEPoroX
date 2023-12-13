/*! \file    HeatConduct.hpp 
 *  \brief   HeatConduct类声明
 *  \author  Shizhe Li 
 *  \date    Aug/26/2023 
 *
 *-----------------------------------------------------------------------------------
 *  版权所有 (C) 2021--present OpenCAEPoroX团队。保留所有权利。
 *  根据GNU Lesser General Public License 3.0或更高版本的条款发布。
 *----------------------------------------------------------------------------------- */

#ifndef __HEATCONDUCT_HEADER__
#define __HEATCONDUCT_HEADER__

// OpenCAEPoroX头文件
#include "ParamReservoir.hpp"
#include "BulkVarSet.hpp"
#include <vector>

using namespace std;

/*! \class HeatConductVarSet
 *  \brief HeatConduct的变量集合
 */
class HeatConductVarSet{
    friend class HeatConduct;
    friend class HeatConductMethod01;

public:
    /*! \brief 设置变量集合
     *  \param nbin bulks的数量
     *  \param npin 相的数量
     *  \param oIndex 油相的索引
     *  \param gIndex 气相的索引
     *  \param wIndex 水相的索引
     */
    void SetUp(const OCP_USI& nbin, const USI& npin, const INT& oIndex, const INT& gIndex, const INT& wIndex) {
        nb = nbin;
        np = npin;
        o  = oIndex;
        g  = gIndex;
        w  = wIndex;
    }

    /*! \brief 恢复到上一个时间步
     */
    void ResetToLastTimeStep() {
        kt  = lkt;
        ktP = lktP;
        ktT = lktT;
        ktS = lktS;
    }

    /*! \brief 更新上一个时间步
     */
    void UpdateLastTimeStep() {
        lkt  = kt;
        lktP = ktP;
        lktT = ktT;
        lktS = ktS;
    }

public:
    OCP_USI         nb;     ///< bulks的数量
    USI             np;     ///< 相的数量
    INT             o, g, w;///< 油、气、水的索引
    vector<OCP_DBL> kt;    ///< 热导率
    vector<OCP_DBL> ktP;   ///< dkt / dP
    vector<OCP_DBL> ktT;   ///< dkt / dT
    vector<OCP_DBL> ktS;   ///< dkt / dS
    vector<OCP_DBL> lkt;   ///< 上一个时间步的热导率
    vector<OCP_DBL> lktP;  ///< 上一个时间步的dkt / dP
    vector<OCP_DBL> lktT;  ///< 上一个时间步的dkt / dT
    vector<OCP_DBL> lktS;  ///< 上一个时间步的dkt / dS
};

/*! \class HeatConductMethod
 *  \brief 热传导方法的基类
 */
class HeatConductMethod{
public:
    /*! \brief 默认构造函数
     */
    HeatConductMethod() = default;

    /*! \brief 计算热传导
     *  \param bId bulks的索引
     *  \param hlvs HeatConductVarSet的引用
     *  \param bvs BulkVarSet的引用
     */
    virtual void CalHeatConduct(const OCP_USI& bId, HeatConductVarSet& hlvs, const BulkVarSet& bvs) const = 0;
};

/*! \class HeatConductMethod01
 *  \brief 热传导方法01的类
 */
class HeatConductMethod01 : public HeatConductMethod{
public:
    /*! \brief 构造函数
     *  \param rs_param Reservoir参数的引用
     *  \param hlvs HeatConductVarSet的引用
     */
    HeatConductMethod01(const ParamReservoir& rs_param, HeatConductVarSet& hlvs);

    /*! \brief 计算热传导
     *  \param bId bulks的索引
     *  \param hlvs HeatConductVarSet的引用
     *  \param bvs BulkVarSet的引用
     */
    void CalHeatConduct(const OCP_USI& bId, HeatConductVarSet& hlvs, const BulkVarSet& bvs) const override;

protected:
    vector<OCP_DBL> thconp; ///< 相的热导率
    OCP_DBL         thconr; ///< 岩石的热导率
};

/*! \class HeatConduct
 *  \brief 热传导的类
 */
class HeatConduct{
public:
    /*! \brief 默认构造函数
     */
    HeatConduct() = default;

    /*! \brief 设置热传导
     *  \param rs_param Reservoir参数的引用
     *  \param bvs BulkVarSet的引用
     */
    void Setup(const ParamReservoir& rs_param, const BulkVarSet& bvs);

    /*! \brief 计算热传导
     *  \param bvs BulkVarSet的引用
     */
    void CalHeatConduct(const BulkVarSet& bvs);

    /*! \brief 获取变量集合
     *  \return 变量集合的引用
     */
    const auto& GetVarSet() const { return vs; }

    /*! \brief 恢复到上一个时间步
     */
    void ResetToLastTimeStep() { if (ifUse)  vs.ResetToLastTimeStep(); }

    /*! \brief 更新上一个时间步
     */
    void UpdateLastTimeStep() { if (ifUse)  vs.UpdateLastTimeStep(); }

protected:
    OCP_BOOL                   ifUse{ OCP_FALSE };            ///< 是否使用热损失
    HeatConductVarSet          vs;                            ///< 热损失变量集合
    vector<HeatConductMethod*> hcM;                           ///< 热传导计算方法
};

#endif /* end if __HEATCONDUCT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/26/2023      Create file                          */
/*----------------------------------------------------------------------------*/
