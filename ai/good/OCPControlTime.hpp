/*! \file    OCPControlTime.hpp
 *  \brief   OCPControlTime类的声明
 *  \author  Shizhe Li
 *  \date    Dec/05/2023
 *
 *  本文件包含了OCPControlTime类的声明，用于控制油田模拟软件中的时间步长选择和控制。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPCONTROLTIME_HEADER__
#define __OCPCONTROLTIME_HEADER__

// Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "ParamControl.hpp"
#include "OCPControlFast.hpp"
#include "OCPNRsuite.hpp"

using namespace std;

/**
 *  \class ControlTimeParam
 *  \brief 选择时间步长的参数类
 *
 *  该类包含了进行时间步长选择时需要的参数。
 */
class ControlTimeParam
{
    friend class ControlTime;

public:
    /// 默认构造函数
    ControlTimeParam() = default;

    /**
     *  \brief 构造函数，只设置控制参数
     *  \param src 时间步长控制参数源
     *  \param Tstep 时间步长序列
     *  \param i 时间步长序列的索引
     */
    ControlTimeParam(const TuningPair& src, const vector<OCP_DBL>& Tstep, const USI& i);

    /**
     *  \brief 设置快速控制参数
     *  \param fCtrl 快速控制对象
     */
    void SetFastControl(const FastControl& fCtrl);

protected:
    OCP_DBL timeInit; ///< 下一个TSTEP开始的第一个时间步长
    OCP_DBL timeMax;  ///< 运行过程中的最大时间步长
    OCP_DBL timeMin;  ///< 运行过程中的最小时间步长
    OCP_DBL maxIncreFac; ///< 最大时间步长增加因子
    OCP_DBL minChopFac;  ///< 最小时间步长削减因子
    OCP_DBL cutFacNR;    ///< 收敛失败后的削减因子
    // 时间步长预测参数，即下一个时间步长的变化限制
    OCP_DBL dPlim;      ///< 理想的最大压力变化
    OCP_DBL dTlim;      ///< 理想的最大温度变化
    OCP_DBL dSlim;      ///< 理想的最大饱和度变化
    OCP_DBL dNlim;      ///< 理想的最大相对Ni（组分摩尔数）变化
    OCP_DBL eVlim;      ///< 理想的最大相对Verr（孔隙-流体）变化
    OCP_DBL begin_time; ///< TSTEP间隔的开始
    OCP_DBL end_time;   ///< TSTEP间隔的结束
};

/**
 *  \class ControlTime
 *  \brief OCP时间控制器类
 *
 *  该类负责控制时间步长，包括时间步长的计算、调整和控制。
 */
class ControlTime
{
public:
    /**
     *  \brief 设置通信器
     *  \param domain 模拟域对象
     */
    void SetupComm(const Domain& domain);

    /**
     *  \brief 设置控制参数
     *  \param src 时间步长控制参数源
     *  \param Tstep 时间步长序列
     *  \param i 时间步长序列的索引
     */
    void SetCtrlParam(const TuningPair& src, const vector<OCP_DBL>& Tstep, const USI& i);

    /**
     *  \brief 设置快速控制参数
     *  \param fCtrl 快速控制对象
     */
    void SetFastControl(const FastControl& fCtrl);

    /**
     *  \brief 削减时间步长
     *  \param fac 削减因子，默认为-1
     */
    void CutDt(const OCP_DBL& fac = -1);

    /**
     *  \brief 根据NR套件削减时间步长
     *  \param NRs NR套件对象
     */
    void CutDt(const OCPNRsuite& NRs);

    /**
     *  \brief 设置下一个TSTEP的参数
     *  \param i 时间步长序列的索引
     *  \param wellOptChange_loc 井优化更改的位置标记
     */
    void SetNextTSTEP(const USI& i, const OCP_BOOL& wellOptChange_loc);

    /**
     *  \brief 计算下一个时间步长
     *  \param NRs NR套件对象
     *  \param il 初始化列表，包含字符串
     */
    void CalNextTimeStep(const OCPNRsuite& NRs, const initializer_list<string>& il);

    /**
     *  \brief 获取总模拟时间
     *  \return 总模拟时间
     */
    auto GetTotalTime() const { return ps.back().end_time; }

    /**
     *  \brief 获取TSTEP间隔的数量
     *  \return TSTEP间隔的数量
     */
    auto GetNumTstepInterval() const { return ps.size(); }

protected:
    MPI_Comm myComm; ///< MPI通信器
    OCP_INT numproc; ///< 进程数量
    OCP_INT myrank;  ///< 当前进程的排名

public:
    /**
     *  \brief 获取当前时间
     *  \return 当前时间
     */
    auto GetCurrentTime() const { return current_time; }

    /**
     *  \brief 获取当前时间步长
     *  \return 当前时间步长
     */
    auto GetCurrentDt() const { return current_dt; }

    /**
     *  \brief 获取上一个时间步长
     *  \return 上一个时间步长
     */
    auto GetLastDt() const { return last_dt; }

    /**
     *  \brief 判断是否达到关键时间点
     *  \return 是否结束的标记
     */
    auto IfEnd() { return ((wp->end_time - current_time) < TINY); }

protected:
    vector<ControlTimeParam> ps; ///< 控制参数集合
    const ControlTimeParam* wp;  ///< 工作参数
    OCP_DBL predict_dt;          ///< 预测的下一个TSTEP的时间步长
    OCP_DBL current_dt;          ///< 当前时间步长
    OCP_DBL last_dt{ 0 };        ///< 上一个时间步长
    OCP_DBL current_time{ 0 };   ///< 当前时间
};

#endif /* end if __OCPControlTime_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/05/2023      Create file                          */
/*----------------------------------------------------------------------------*/
