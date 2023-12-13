```cpp
/*! \file    OCPOutput.hpp
 *  \brief   本文件包含OCPOutput类的声明，用于管理OpenCAEPoroX软件的输出信息。
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCP_OUTPUT_HEADER__
#define __OCP_OUTPUT_HEADER__

// 标准头文件
#include <iomanip>
#include <iostream>
#include <cstdio>

// OpenCAEPoroX头文件
#include "OCPControl.hpp"
#include "Output4Vtk.hpp"
#include "ParamOutput.hpp"
#include "Reservoir.hpp"
#include "UtilOutput.hpp"
#include "UtilTiming.hpp"

using namespace std;

/// \class OCPIJK
/// \brief 3D坐标表示，用于OpenCAEPoroX软件。
class OCPIJK
{
public:
    OCPIJK() = default;
    /*! \brief 构造函数
     *
     *  \param i 沿I方向的索引
     *  \param j 沿J方向的索引
     *  \param k 沿K方向的索引
     */
    OCPIJK(const USI& i, const USI& j, const USI& k)
        : I(i), J(j), K(k){};
    
    /*! \brief 拷贝构造函数
     *
     *  \param src 源OCPIJK对象
     */
    OCPIJK(const COOIJK& src)
    {
        I = src.I;
        J = src.J;
        K = src.K;
    };

    /*! \brief 赋值运算符重载
     *
     *  \param src 源COOIJK对象
     *  \return 返回当前对象的引用
     */
    OCPIJK& operator=(const COOIJK& src)
    {
        I = src.I;
        J = src.J;
        K = src.K;
        return *this;
    }

    USI I; ///< 沿I方向的索引
    USI J; ///< 沿J方向的索引
    USI K; ///< 沿K方向的索引
};

/// \class OCPType_Sum
/// \brief 用于记录和汇总特定类型数据的类模板。
template <typename T>
class OCPType_Sum
{
public:
    /*! \brief 赋值运算符重载
     *
     *  \param src 源Type_A_o对象
     *  \return 返回当前对象的引用
     */
    OCPType_Sum& operator=(const Type_A_o& src)
    {
        activity = src.activity;
        obj      = src.obj;
        return *this;
    }

    /*! \brief 赋值运算符重载
     *
     *  \param src 源Type_B_o对象
     *  \return 返回当前对象的引用
     */
    OCPType_Sum& operator=(const Type_B_o& src)
    {
        activity = src.activity;
        obj.assign(src.obj.begin(), src.obj.end());
        return *this;
    }

    OCP_BOOL  activity{OCP_FALSE}; ///< 活动状态标志，默认为不活动
    vector<T> obj;                 ///< 存储T类型对象的向量
    vector<USI> index;             ///< 记录将要打印属性的批量或井的索引
};

/// \class ItersInfo
/// \brief 记录迭代信息的类。
class ItersInfo
{
public:
    /*! \brief 更新所有迭代信息
     *
     *  \param NRs Newton-Raphson迭代信息
     */
    void Update(const OCPNRsuite& NRs);

    /*! \brief 获取时间步数
     *
     *  \return 时间步数
     */
    auto GetNumTimeStep() const { return numTstep; }

    /*! \brief 获取总的Newton迭代次数
     *
     *  \return 总的Newton迭代次数
     */
    auto GetNRt() const { return NRt; }

    /*! \brief 获取浪费的Newton迭代次数
     *
     *  \return 浪费的Newton迭代次数
     */
    auto GetNRwt() const { return NRwt; }

    /*! \brief 获取总的线性迭代次数
     *
     *  \return 总的线性迭代次数
     */
    auto GetLSt() const { return LSt; }

    /*! \brief 获取浪费的线性迭代次数
     *
     *  \return 浪费的线性迭代次数
     */
    auto GetLSwt() const { return LSwt; }

protected:
    USI numTstep{ 0 }; ///< 时间步数
    USI NRt{ 0 };      ///< 总的Newton迭代次数
    USI NRwt{ 0 };     ///< 浪费的Newton迭代次数
    USI LSt{ 0 };      ///< 总的线性迭代次数
    USI LSwt{ 0 };     ///< 浪费的线性迭代次数
};

/// \class SumItem
/// \brief 汇总数据输出的辅助结构类。
class SumItem
{
public:
    SumItem() = default;

    /*! \brief 构造函数
     *
     *  \param item 项目名称
     *  \param obj 对象名称
     */
    SumItem(const string& item, const string& obj)
        : Item(item), Obj(obj) {};

    /*! \brief 构造函数
     *
     *  \param item 项目名称
     *  \param n 值的数量
     */
    SumItem(const string& item, const USI& n)
        : Item(item) { val.reserve(n); };

    /*! \brief 构造函数
     *
     *  \param item 项目名称
     *  \param obj 对象名称
     *  \param unit 单位
     *  \param type 类型
     *  \param n 值的数量
     */
    SumItem(const string& item,
            const string& obj,
            const string& unit,
            const string& type,
            const OCP_USI& n)
        : Item(item)
        , Obj(obj)
        , Unit(unit)
        , Type(type) { val.reserve(n); };

    /*! \brief 等于运算符重载
     *
     *  \param other 另一个SumItem对象
     *  \return 如果项目名称和对象名称都相同，则返回真，否则返回假
     */
    OCP_BOOL operator==(const SumItem& other)
    {
        if (this->Item == other.Item && this->Obj == other.Obj)
            return OCP_TRUE;
        else
            return OCP_FALSE;
    }

    string          Item; ///< 项目名称
    string          Obj;  ///< 对象名称
    string          Unit; ///< 单位
    string          Type; ///< 类型
    vector<OCP_DBL> val;  ///< 存储数值的向量
};

/// \class SUMMARY文件管理
/// \brief 负责管理要在SUMMARY文件中输出的内容
class Summary
{
public:
    /// 读入要输出的变量名
    void InputParam(const OutputSummary& summary_param);

    /// 组装输出操作
    void Setup(const Reservoir& reservoir);

    /// 设置要输出的变量的值
    void SetVal(const Reservoir& reservoir, const OCPControl& ctrl, const ItersInfo& iters);

    /// 将内容打印到文件
    void PrintInfo(const string& dir, const string& filename, const OCP_INT& rank) const;

    /// 并行运行下主进程合并各进程的输出文件
    void PostProcess(const string& dir, const string& filename, const OCP_INT& numproc) const;

protected:
    vector<SumItem> Sumdata;    ///< 所有要打印的信息

    OCP_BOOL FPR{ OCP_FALSE };  ///< 油田平均压力
    OCP_BOOL FTR{ OCP_FALSE };  ///< 油田平均温度
    OCP_BOOL FOPR{ OCP_FALSE }; ///< 油田产油速率
    OCP_BOOL FOPT{ OCP_FALSE }; ///< 油田产油总量
    OCP_BOOL FGPR{ OCP_FALSE }; ///< 油田产气速率
    OCP_BOOL FGPt{ OCP_FALSE }; ///< 油田产气总量
    OCP_BOOL FWPR{ OCP_FALSE }; ///< 油田产水速率
    OCP_BOOL FWPT{ OCP_FALSE }; ///< 油田产水总量
    OCP_BOOL FGIR{ OCP_FALSE }; ///< 油田注气速率
    OCP_BOOL FGIT{ OCP_FALSE }; ///< 油田注气总量
    OCP_BOOL FWIR{ OCP_FALSE }; ///< 油田注水速率
    OCP_BOOL FWIT{ OCP_FALSE }; ///< 油田注水总量

    OCPType_Sum<string> WOPR; ///< 井产油速率
    OCPType_Sum<string> WOPT; ///< 井产油总量
    OCPType_Sum<string> WGPR; ///< 井产气速率
    OCPType_Sum<string> WGPT; ///< 井产气总量
    OCPType_Sum<string> WWPR; ///< 井产水速率
    OCPType_Sum<string> WWPT; ///< 井产水总量
    OCPType_Sum<string> WGIR; ///< 井注气速率
    OCPType_Sum<string> WGIT; ///< 井注气总量
    OCPType_Sum<string> WWIR; ///< 井注水速率
    OCPType_Sum<string> WWIT; ///< 井注水总量
    OCPType_Sum<string> WBHP; ///< 井底压力
    OCPType_Sum<string> DG;   ///< 井射孔压差

    OCPType_Sum<OCPIJK> BPR;  ///< 网格块压力
    OCPType_Sum<OCPIJK> SOIL; ///< 网格块油相饱和度
    OCPType_Sum<OCPIJK> SGAS; ///< 网格块气相饱和度
    OCPType_Sum<OCPIJK> SWAT; ///< 网格块水相饱和度
};

/// \class 关键信息输出文件管理
/// \brief 负责管理输出运行过程中的关键信息
class CriticalInfo
{
public:
    /// 组装输出
    void Setup();

    /// 输出变量赋值
    void SetVal(const OCPControl& ctrl, const OCPNRsuite& NR);

    /// 打印到文件
    void PrintFastReview(const string& dir, const string& filename, const OCP_INT& rank, const ItersInfo& iters) const;

    /// 并行运行下主进程合并其余进程输出文件
    void PostProcess(const string& dir, const string& filename, const OCP_INT& numproc, const ItersInfo& iters) const;

private:
    vector<SumItem> Sumdata; ///< 输出变量信息
};

/// \class 基于网格的可输出基本信息
/// \brief 用于管理基于网格的可输出基本信息
class BasicGridProperty
{
    friend class Out4RPT;
    friend class Out4VTK;

public:
    void SetBasicGridProperty(const BasicGridPropertyParam& param);
    void Check(const Reservoir& rs);

private:
    OCP_BOOL PRE{ OCP_FALSE };    ///< 网格压力
    OCP_BOOL PGAS{ OCP_FALSE };   ///< 气相压力
    OCP_BOOL PWAT{ OCP_FALSE };   ///< 水相压力
    OCP_BOOL SOIL{ OCP_FALSE };   ///< 油相饱和度
    OCP_BOOL SGAS{ OCP_FALSE };   ///< 气相饱和度
    OCP_BOOL SWAT{ OCP_FALSE };   ///< 水相饱和度
    OCP_BOOL DENO{ OCP_FALSE };   ///< 油相密度
    OCP_BOOL DENG{ OCP_FALSE };   ///< 气相密度
    OCP_BOOL DENW{ OCP_FALSE };   ///< 水相密度
    OCP_BOOL KRO{ OCP_FALSE };    ///< 油相相对渗透率
    OCP_BOOL KRG{ OCP_FALSE };    ///< 气相相对渗透率
    OCP_BOOL KRW{ OCP_FALSE };    ///< 水相相对渗透率
    OCP_BOOL BOIL{ OCP_FALSE };   ///< 油相体积转化因子
    OCP_BOOL BGAS{ OCP_FALSE };   ///< 气相体积转化因子
    OCP_BOOL BWAT{ OCP_FALSE };   ///< 水相体积转化因子
    OCP_BOOL VOIL{ OCP_FALSE };   ///< 油相粘度
    OCP_BOOL VGAS{ OCP_FALSE };   ///< 气相粘度
    OCP_BOOL VWAT{ OCP_FALSE };   ///< 水相粘度
    OCP_BOOL XMF{ OCP_FALSE };    ///< 油相各组分摩尔占比
    OCP_BOOL YMF{ OCP_FALSE };    ///< 气相各组分摩尔占比
    OCP_BOOL PCW{ OCP_FALSE };    ///< 油水毛管力
    OCP_BOOL CO2{ OCP_FALSE };    ///< 水相中CO2浓度
    OCP_BOOL SATNUM{ OCP_FALSE }; ///< 饱和度区域分布

    USI             bgpnum;       ///< 输出变量数量
};

/// \class 打印每个网格指定变量的详细信息
/// \brief 打印每个网格指定变量的详细信息，仅适用结构网格
class Out4RPT
{
public:
    /// 读入需打印的变量信息
    void InputParam(const OutputRPTParam& RPTparam);
    /// 组装打印操作
    void Setup(const string& dir, const Reservoir& reservoir);
    /// 打印到文件
    void PrintRPT(const string& dir, const Reservoir& rs, const OCP_DBL& days) const;
    /// 获取被打印网格的坐标
    void GetIJKGrid(USI& i, USI& j, USI& k, const OCP_ULL& n) const;

private:
    OCP_BOOL          useRPT{ OCP_FALSE };        ///< 是否使用
    OCP_ULL           numGrid;                    ///< 网格总数
    OCP_USI           nx;                         ///< x方向网格数量
    OCP_USI           ny;                         ///< y方向网格数量
    USI               IJKspace;                   ///< 输出间隙
    BasicGridProperty bgp;                        ///< 基本网格信息输出控制
};


/// \class VTK输出
/// \brief 以VTK文件的格式打印每个网格指定变量的详细信息，仅适用结构网格
class Out4VTK
{
public:
    /// 输入要打印的变量名
    void InputParam(const OutputVTKParam& VTKParam);
    /// 组装打印操作
    void Setup(const string& dir, const Reservoir& rs);
    /// 打印到文件
    void PrintVTK(const Reservoir& rs) const;
    /// 输出文件后处理
    void PostProcess(const string& dir, const string& filename, const OCP_INT& numproc) const;
    /// 多进程运行时整合各进程输出并后处理
    void PostProcessP(const string& dir, const string& filename, const OCP_INT& numproc) const;
    /// 单进程运行时后处理
    void PostProcessS(const string& dir, const string& filename) const;

private:
    OCP_BOOL          useVTK{ OCP_FALSE }; ///< 是否使用
    BasicGridProperty bgp;                 ///< 基本网格信息输出控制
    Output4Vtk        out4vtk;             ///< VTK输出模块
    string            myFile;              ///< 输出文件名
    mutable OCP_ULL   numGrid;             ///< 网格数量
};

/// \class 输出管理模块
/// \brief 管理各类输出
class OCPOutput
{
    friend class OpenCAEPoroX;

public:
    /// 读入要输出变量名
    void InputParam(const ParamOutput& paramOutput);
    /// 组装各类别输出
    void Setup(const Reservoir& reservoir, const OCPControl& ctrl, const Domain& domain);
    /// 对输出变量赋值
    void SetVal(const Reservoir& reservoir, const OCPControl& ctrl, const OCPNRsuite& NR);
    /// 打印基于时间步的输出
    void PrintInfo() const;
    /// 打印基于TSTEP的输出
    void PrintInfoSched(const Reservoir& rs, const OCPControl& ctrl, const OCP_DBL& time) const;
    /// 后处理输出文件
    void PostProcess() const;
    /// 打印当前时间步迭代步信息
    void PrintCurrentTimeIter(const OCPControl& ctrl) const;

protected:
    MPI_Comm  myComm;              ///< 通信器
    OCP_INT   numproc, myrank;     ///< 进程数量，进程编号


protected:
    string       workDir;          ///< 工作文件夹
    string       fileName;         ///< 文件名
    ItersInfo    iters;            ///< 迭代步信息
    Summary      summary;          ///< SUMMARY输出
    CriticalInfo crtInfo;          ///< 运行时信息输出
    Out4RPT      out4RPT;          ///< RPT输出
    Out4VTK      out4VTK;          ///< VTK输出
};

#endif /* end if __OCPOUTPUT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/08/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/
