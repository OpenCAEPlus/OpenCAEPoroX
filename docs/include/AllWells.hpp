/*! 
\file    AllWells.hpp 
\brief   AllWells类声明
\author  Shizhe Li
\date    Oct/01/2021
*/

#ifndef __WELLGROUP_HEADER__
#define __WELLGROUP_HEADER__

#include "ParamReservoir.hpp"
#include "Well.hpp"
#include "WellPeaceman.hpp"

using namespace std;

/// \class WellGroup
/// \brief 包含一个井组，负责管理一些井的生产和注入目标及其互动等，它会在模拟开始时初始化，如果需要，应该更新，例如，井的类型更改或井被重新分组。
class WellGroup
{
    friend class AllWells;

public:
    /// \brief 默认构造函数
    WellGroup() = default;

    /// \brief 构造函数
    /// \param gname 井组名称
    WellGroup(const string& gname)
        : name(gname){};

private:
    string      name;    ///< 井组名称
    vector<USI> wId;     ///< 井组中的井索引
    vector<USI> wIdINJ;  ///< AllWells中的注入井索引
    vector<USI> wIdPROD; ///< AllWells中的生产井索引
};

/// \class AllWells
/// \brief 包含所有现在的井，用于统一管理储层中的所有井。实际上，你可以把它看作是井和其他模块之间的接口。
class AllWells
{
    friend class Reservoir;
    friend class OCPNRsuite;
    friend class Out4RPT;
    friend class Out4VTK;
    friend class IsoT_FIM;
    friend class IsoT_IMPEC;
    friend class IsoT_AIMc;
    friend class IsoT_FIMddm;
    friend class T_FIM;

public:
    /// \brief 默认构造函数
    AllWells() = default;

    /// \brief 输入井参数
    /// \param paramWell 井参数
    /// \param domain 域
    void InputParam(const ParamWell& paramWell, const Domain& domain);

    /// \brief 设置井
    /// \param bk Bulk对象
    void Setup(const Bulk& bk);

    /// \brief 设置井组信息
    /// \param bk Bulk对象
    void SetupWellGroup(const Bulk& bk);

    /// \brief 设置当前阶段被井穿透的Bulk
    /// \param bk Bulk对象
    void SetupWellBulk(Bulk& bk) const;

    /// \brief 应用第i个关键时间的操作模式
    /// \param i 时间索引
    void ApplyControl(const USI& i);

    /// \brief 设置初始井压力
    /// \param bk Bulk对象
    void InitBHP(const Bulk& bk);

    /// \brief 计算每个时间步开始时的井属性
    /// \param bk Bulk对象
    void PrepareWell(const Bulk& bk);

    /// \brief 计算每个射孔的体积流速和物质流速
    /// \param bk Bulk对象
    void CalFlux(const Bulk& bk);

    /// \brief 计算注入率、总注入、生产率、总生产
    /// \param bk Bulk对象
    /// \param dt 时间步长
    void CalIPRT(const Bulk& bk, OCP_DBL dt);

    /// \brief 检查是否发生了不合理的井压力或射孔压力
    /// \param bk Bulk对象
    /// \return 储层状态
    ReservoirState CheckP(const Bulk& bk);

    /// \brief 返回井的数量
    /// \return 井的数量
    USI GetWellNum() const { return numWell; }

    /// \brief 返回指定井的名称
    /// \param i 井索引
    /// \return 井的名称
    string GetWellName(const USI& i) const { return wells[i]->name; }

    /// \brief 返回指定井的索引
    /// \param name 井的名称
    /// \return 井的索引
    USI GetIndex(const string& name) const;

    /// \brief 返回井i的射孔数量
    /// \param i 井索引
    /// \return 射孔数量
    USI GetWellPerfNum(const USI& i) const { return wells[i]->numPerf; }

    /// \brief 返回所有井的射孔数量
    /// \return 射孔数量
    USI GetWellPerfNum() const;

    /// \brief 计算所有井的最大射孔数量
    /// \return 最大射孔数量
    USI     GetMaxWellPerNum() const;

    /// \brief 返回wth井的BHP
    /// \param w 井索引
    /// \return BHP
    OCP_DBL GetWBHP(const USI& w) const
    {
        if (wells[w]->IsOpen())  return wells[w]->bhp;
        else                     return 0;
    }

    /// \brief 显示井状态
    /// \param bk Bulk对象
    void    ShowWellStatus(const Bulk& bk)
    {
        for (USI w = 0; w < numWell; w++) wells[w]->ShowPerfStatus(bk);
    }

    /// \brief 获取井操作是否改变
    /// \return 井操作是否改变
    OCP_BOOL    GetWellOptChange() const { return wellOptChange; }

    /// \brief 获取打开的井的数量
    /// \return 打开的井的数量
    USI GetNumOpenWell()const;

protected:
    USI                     numWell;   ///< 井的数量
    vector<Well*>           wells;     ///< 井集合
    USI                     numGroup;  ///< 组的数量
    vector<WellGroup>       wellGroup; ///< 井组集合
    OCP_BOOL           wellOptChange; ///< 如果井发生变化，则为OCP_TRUE
    vector<SolventINJ> solvents;   ///< 溶剂集合
    OCP_DBL            dPmax{0};   ///< 最大BHP改变
    OCP_DBL          Psurf;    ///< 井参考压力
    OCP_DBL          Tsurf;    ///< 井参考温度
};
#endif /* end if __WELLGROUP_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Shizhe Li           Feb/08/2022      Rename to AllWells                   */
/*----------------------------------------------------------------------------*/