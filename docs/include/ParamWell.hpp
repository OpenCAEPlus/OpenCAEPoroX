```cpp
/*! 
 * \file    ParamWell.hpp
 * \brief   ParamWell类声明文件
 * \author  Shizhe Li
 * \date    Oct/01/2021
 *
 * 本文件包含ParamWell类及其相关类的声明，这些类用于存储和处理油井的参数信息。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARAMWELL_HEADER__
#define __PARAMWELL_HEADER__

// 标准头文件
#include <cassert>
#include <fstream>
#include <vector>

// OpenCAEPoroX头文件
#include "OCPConst.hpp"
#include "UtilInput.hpp"
#include "UtilOutput.hpp"

using namespace std;

/*!
 * \class   WellOptParam
 * \brief   包含油井操作所需的所有参数
 *
 * WellOptParam类存储了一个油井在不同时间可能改变的所有操作参数。
 */
class WellOptParam
{
public:
    /*!
     * \brief   构造函数，初始化油井操作参数
     * \param   intype 操作类型，注入或生产
     * \param   vbuf   参数值的字符串向量
     */
    WellOptParam(string intype, vector<string>& vbuf);
    
    string type;      ///< 油井类型，注入还是生产？
    string fluidType; ///< 注入井的注入流体类型（仅限注入井）
    string state;     ///< 油井状态，开启或关闭？
    string mode;      ///< 油井模式，按率控还是按底压控？
    OCP_DBL maxRate;  ///< 允许的最大流入/流出油井的流速
    OCP_DBL maxBHP;   ///< 注入井允许的最大井底压力
    OCP_DBL minBHP;   ///< 生产井允许的最小井底压力
    OCP_DBL injTemp;  ///< 注入流体的温度
};

/*!
 * \class   WellOptPair
 * \brief   包含油井操作模式及其应用的开始时间
 *
 * WellOptPair类包含两部分信息：一部分是油井的操作模式，另一部分是操作开始应用的时间，
 * 该时间通过在关键时间点的索引来表示。
 */
class WellOptPair
{
public:
    /*!
     * \brief   构造函数，初始化油井操作模式和开始时间
     * \param   i     关键时间点的索引
     * \param   type  操作类型
     * \param   vbuf  参数值的字符串向量
     */
    WellOptPair(USI i, string type, vector<string>& vbuf);
    
    USI          d;    ///< 操作开始应用的时间点的索引
    WellOptParam opt;  ///< 油井的操作参数
};

/*!
 * \class   WellParam
 * \brief   包含油井静态和动态参数信息
 *
 * WellParam类用于存储和处理油井的静态信息（如井的位置和尺寸）和动态信息（如操作模式）。
 */
class WellParam
{
public:
    /// 输入结构网格下的井参数
    WellParam(vector<string>& info);
    /// 输入非结构网格下的井参数
    WellParam(vector<string>& info, const string& unstructured);
    /// 输入COMPDAT关键字
    void InputCOMPDAT(vector<string>& vbuf);
    /// 输入结构网格下的COMPDAT关键字
    void InputCOMPDATS(vector<string>& vbuf);
    /// 输入非结构网格下的COMPDAT关键字
    void InputCOMPDATUS(vector<string>& vbuf);
    /// 获取射孔数量
    USI GetPerfNum() const { return max(I_perf.size(), X_perf.size()); }

    GridType gridType;             ///< 网格类型
    string   name;                 ///< 井名字
    string   group{ "FEILD" };     ///< 井所在井组
    USI      I;                    ///< 结构网格下X方向索引
    USI      J;                    ///< 结构网格下Y方向索引
    OCP_DBL  depth{ -1.0 };        ///< 井参考深度
    OCP_DBL  X;                    ///< 非结构网格下X坐标
    OCP_DBL  Y;                    ///< 非结构网格下Y坐标
    OCP_DBL  Z;                    ///< 非结构网格下Z坐标

    vector<USI>     I_perf;        ///< 结构网格下射孔X方向索引
    vector<USI>     J_perf;        ///< 结构网格下射孔Y方向索引
    vector<USI>     K_perf;        ///< 结构网格下射孔Z方向索引
    vector<OCP_DBL> X_perf;        ///< 非结构网格下射孔X坐标
    vector<OCP_DBL> Y_perf;        ///< 非结构网格下射孔Y坐标
    vector<OCP_DBL> Z_perf;        ///< 非结构网格下射孔Z坐标
    vector<OCP_DBL> WI;            ///< 井射孔指数
    vector<OCP_DBL> diameter;      ///< 井射孔直径
    vector<OCP_DBL> kh;            ///< 井射孔kh系数
    vector<OCP_DBL> skinFactor;    ///< 井射孔因子
    vector<string>  direction;     ///< 井射孔方向
    OCP_BOOL        ifUseUnweight{ OCP_FALSE }; ///< 是否使用非加权因子

    // dynamic infomation
    vector<WellOptPair> optParam;
};

/// 注入物数据结构
class Solvent
{
public:
    /// 默认构造函数
    Solvent() = default;
    /// 构造函数
    Solvent(const vector<string>& vbuf);
    /// 井名字
    string          name;
    /// 组成比例
    vector<OCP_DBL> comRatio;
};

/// 井参数
class ParamWell
{
public:
    OCP_BOOL          thermal{ OCP_FALSE };   ///< 是否使用热流
    vector<WellParam> well;                   ///< 井参数
    vector<OCP_DBL>   criticalTime;           ///< TSETP时间
    vector<Solvent>   solSet;                 ///< 注入物序列
    OCP_DBL           Psurf;                  ///< 表面压力
    OCP_DBL           Tsurf;                  ///< 表面温度

    /// 初始化井参数
    void Init();
    /// 初始化TSTEP
    void InitTime() { criticalTime.push_back(0); };
    /// 输入WELSPECS
    void InputWELSPECS(ifstream& ifs);
    /// 输入COMPDAT
    void InputCOMPDAT(ifstream& ifs);
    /// 输入WCONINJE
    void InputWCONINJE(ifstream& ifs);
    /// 输入WCONPROD
    void InputWCONPROD(ifstream& ifs);
    /// 输入TSTEP
    void InputTSTEP(ifstream& ifs);
    /// 输入WELTARG
    void InputWELTARG(ifstream& ifs);
    /// 输入WTEMP
    void InputWTEMP(ifstream& ifs);
    /// 输入UNWEIGHT
    void InputUNWEIGHT(ifstream& ifs);
    /// 输入WELLSTRE
    void InputWELLSTRE(ifstream& ifs);
    /// 输入PSURF
    void InputPSURF(ifstream& ifs);
    /// 输入TSURF
    void InputTSURF(ifstream& ifs);
    /// 检查输入参数
    void CheckParam() const;
    /// 检查输入射孔
    void CheckPerf() const;
    /// 检查注入流体
    void CheckINJFluid() const;
};

#endif /* end if __PARAMWELL_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/
```

