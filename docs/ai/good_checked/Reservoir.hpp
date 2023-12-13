/*! \file    Reservoir.hpp
 *  \brief   油藏类声明
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __RESERVOIR_HEADER__
#define __RESERVOIR_HEADER__

// OpenCAEPoroX header files
#include "AllWells.hpp"
#include "Bulk.hpp"
#include "BulkConn.hpp"
#include "ParamRead.hpp"
#include "PreProcess.hpp"

/*!
 * \class VarInfo
 * \brief 用于存储变量信息的模板类
 *
 * \tparam T 数据类型
 */
template <typename T>
class VarInfo {
public:
    /*!
     * \brief 构造函数，用于初始化变量信息
     * \param n 变量名称
     * \param src 源数据指针
     * \param dst 目标数据指针
     */
    VarInfo(const string n, const T* src, T* dst) :name(n), src_ptr(src), dst_ptr(dst) {};

    /*!
     * \brief 构造函数，用于初始化变量信息，目标指针设为源指针
     * \param n 变量名称
     * \param src 源数据指针
     */
    VarInfo(const string n, const T* src) :name(n), src_ptr(src) { dst_ptr = const_cast<T*>(src_ptr); };

    const string   name;      ///< 变量名称
    const T*       src_ptr;   ///< 源数据指针
    T*             dst_ptr;   ///< 目标数据指针
};

/*!
 * \class VarInfoGrid
 * \brief 用于存储网格变量信息的类
 */
class VarInfoGrid{
public:
    /*!
     * \brief 默认构造函数
     */
    VarInfoGrid() = default;

    /*!
     * \brief 构造函数，用于初始化网格变量信息
     * \param dbl 双精度类型变量信息列表
     * \param usi 无符号短整型变量信息列表
     */
    VarInfoGrid(const vector<VarInfo<vector<OCP_DBL>>> dbl,
                const vector<VarInfo<vector<OCP_USI>>> usi):
        var_dbl(dbl), var_usi(usi) { CalNumByte(); };

    /*!
     * \brief 计算变量占用的字节数
     */
    void CalNumByte() {
        numvar_dbl    = var_dbl.size();
        numvar_usi    = var_usi.size();
        numByte_dbl   = sizeof(OCP_DBL);
        numByte_usi   = sizeof(OCP_USI);
        numByte_total = numvar_dbl * numByte_dbl + numvar_usi * numByte_usi;
    }

    vector<VarInfo<vector<OCP_DBL>>>  var_dbl;    ///< 双精度类型变量信息列表
    vector<VarInfo<vector<OCP_USI>>>  var_usi;    ///< 无符号短整型变量信息列表
    OCP_INT  numvar_dbl;                          ///< 双精度类型变量数量
    OCP_INT  numvar_usi;                          ///< 无符号短整型变量数量
    OCP_INT  numByte_dbl;                         ///< 双精度类型变量占用字节数
    OCP_INT  numByte_usi;                         ///< 无符号短整型变量占用字节数
    OCP_INT  numByte_total;                       ///< 总占用字节数
};

/*!
 * \class Reservoir
 * \brief 油藏核心组件类，包含所有油藏信息以及相关操作。
 *
 * 油藏有四个核心组件：
 * Grids 包含作为油藏数据库的所有网格的基本信息。
 * Bulk 仅存储活动网格，定义了用于计算的区域。
 * AllWells 包含井信息，用于管理与井相关的操作。
 * BulkConn 包含活动网格之间的连接。
 */
class Reservoir{
    friend class OCPControl;
    friend class ControlTime;
    friend class Summary;
    friend class CriticalInfo;
    friend class Out4RPT;
    friend class Out4VTK;
    friend class OCPNRsuite;
    friend class IsothermalMethod;
    friend class IsoT_IMPEC;
    friend class IsoT_FIM;
    friend class IsoT_AIMc;
    friend class IsoT_FIMddm;
    friend class T_FIM;
    friend class Solver;

    /////////////////////////////////////////////////////////////////////
    // Param Distribute
    /////////////////////////////////////////////////////////////////////
public:
    /*!
     * \brief 输入参数
     * \param prepro 预处理类
     * \param rsparam 参数读取类
     */
    void InputParam(PreProcess& prepro, ParamRead& rsparam);

    /*!
     * \brief 设置域
     * \param myDomain 域类
     */
    void SetupDomain(Domain& myDomain) { swap(domain, myDomain); }

    /*!
     * \brief 输入基于网格的分布参数
     * \param prepro 预处理网格井类
     */
    void InputDistParamGrid(PreParamGridWell& prepro);

    /*!
     * \brief 输入基于井的分布参数
     * \param param 参数读取类
     */
    void InputDistParamOthers(const ParamRead& param);

    /*!
     * \brief 获取域信息
     * \return 域类的常量引用
     */
    const Domain& GetDomain() const { return domain; }

protected:
    /*!
     * \brief 获取发送变量的名称和地址
     * \param varInfo 变量信息类
     * \param send_var 发送变量信息类
     * \return 整型，表示操作状态
     */
    OCP_INT GetSendVarInfo(const VarInfoGrid& varInfo, VarInfoGrid& send_var);

    /////////////////////////////////////////////////////////////////////
    // General
    /////////////////////////////////////////////////////////////////////
public:
    /*!
     * \brief 设置油藏的静态信息
     */
    void Setup();

    /*!
     * \brief 应用第i个临界时间点的控制
     * \param i 临界时间点索引
     */
    void ApplyControl(const USI& i);

    /*!
     * \brief 计算注入和生产的数量
     * \param dt 时间步长
     */
    void CalIPRT(const OCP_DBL& dt);

    /*!
     * \brief 获取活动网格数量
     * \return 活动网格数量
     */
    OCP_USI GetBulkNum() const { return bulk.GetBulkNum(); }

    /*!
     * \brief 获取内部活动网格数量
     * \return 内部活动网格数量
     */
    OCP_USI GetInteriorBulkNum() const { return bulk.GetInteriorBulkNum(); }

    /*!
     * \brief 获取井数量
     * \return 井数量
     */
    USI GetWellNum() const { return allWells.GetWellNum(); }

    /*!
     * \brief 获取组分数量
     * \return 组分数量
     */
    USI GetComNum() const { return bulk.GetComNum(); }

    /*!
     * \brief 获取开放井数量
     * \return 开放井数量
     */
    USI GetNumOpenWell() const { return allWells.GetNumOpenWell(); }

    /*!
     * \brief 判断是否存在油
     * \return 布尔值，存在油则为真
     */
    OCP_BOOL IfOilExist() const { return bulk.vs.o >= 0; }

    /*!
     * \brief 判断是否存在气
     * \return 布尔值，存在气则为真
     */
    OCP_BOOL IfGasExist() const { return bulk.vs.g >= 0; }

    /*!
     * \brief 判断是否存在水
     * \return 布尔值，存在水则为真
     */
    OCP_BOOL IfWatExist() const { return bulk.vs.w >= 0; }

protected:
    Bulk             bulk;        ///< 活动网格信息
    AllWells         allWells;    ///< 井类信息
    BulkConn         conn;        ///< 活动网格的连接信息
    Domain           domain;      ///< 域分解信息

public:
    /*!
     * \brief 打印FIM解决方案
     * \param outfile 输出文件名
     */
    void    PrintSolFIM(const string& outfile) const;

    /*!
     * \brief 输出最终信息
     */
    void    OutInfoFinal() const { if (domain.numproc == 1) bulk.OutMixtureIters(); }
};

#endif /* end if __RESERVOIR_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/
