/*! \file    ParamOutput.hpp
 *  \brief   本文件包含了ParamOutput类的声明与相关辅助结构体的定义。
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *  本文件定义了用于输出模拟结果的结构体和类，包括输出的参数设置、输出内容的类型等。
 *  这些类和结构体将在模拟过程中用于控制和指示输出结果的详细程度。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARAMOUTPUT_HEADER__
#define __PARAMOUTPUT_HEADER__

// 标准头文件
#include <fstream>
#include <vector>

// OpenCAEPoroX头文件
#include "UtilInput.hpp"

/// \class COOIJK
/// \brief 三维坐标结构体。
class COOIJK{
public:
    COOIJK() = default;
    COOIJK(const USI& i, const USI& j, const USI& k)
        : I(i), J(j), K(k){};
    USI I; ///< I方向坐标。
    USI J; ///< J方向坐标。
    USI K; ///< K方向坐标。
};

/// \class Type_B_o
/// \brief 用于存储坐标形式内容的关键字。
class Type_B_o{
public:
    OCP_BOOL       activity{OCP_FALSE}; ///< 是否活跃。
    vector<COOIJK> obj;                 ///< 坐标对象的向量。
};

/// \class Type_A_o
/// \brief 用于存储字符串形式内容的关键字。
class Type_A_o{
public:
    OCP_BOOL       activity{OCP_FALSE}; ///< 是否活跃。
    vector<string> obj;                 ///< 字符串对象的向量。
};

/// \class OutputSummary
/// \brief 包含所有SUMMARY信息，这些信息指示了哪些结果将被打印到摘要输出文件中。
class OutputSummary{
public:
    // 下面的注释说明了每个成员变量的基本功能和默认值
    OCP_BOOL FPR{OCP_FALSE};  ///< 场均压力。
    OCP_BOOL FTR{OCP_FALSE};  ///< 场均温度。
    OCP_BOOL FOPR{OCP_FALSE}; ///< 场原油生产速率。
    OCP_BOOL FOPT{OCP_FALSE}; ///< 场累计原油生产量。
    OCP_BOOL FGPR{OCP_FALSE}; ///< 场天然气生产速率。
    OCP_BOOL FGPt{OCP_FALSE}; ///< 场累计天然气生产量。
    OCP_BOOL FWPR{OCP_FALSE}; ///< 场水生产速率。
    OCP_BOOL FWPT{OCP_FALSE}; ///< 场累计水生产量。
    OCP_BOOL FGIR{OCP_FALSE}; ///< 场天然气注入速率。
    OCP_BOOL FGIT{OCP_FALSE}; ///< 场累计天然气注入量。
    OCP_BOOL FWIR{OCP_FALSE}; ///< 场水注入速率。
    OCP_BOOL FWIT{OCP_FALSE}; ///< 场累计水注入量。
    Type_A_o WOPR; ///< 井原油生产速率。
    Type_A_o WOPT; ///< 井累计原油生产量。
    Type_A_o WGPR; ///< 井天然气生产速率。
    Type_A_o WGPT; ///< 井累计天然气生产量。
    Type_A_o WWPR; ///< 井水生产速率。
    Type_A_o WWPT; ///< 井累计水生产量。
    Type_A_o WGIR; ///< 井天然气注入速率。
    Type_A_o WGIT; ///< 井累计天然气注入量。
    Type_A_o WWIR; ///< 井水注入速率。
    Type_A_o WWIT; ///< 井累计水注入量。
    Type_A_o WBHP; ///< 井压力。
    Type_A_o DG;   ///< 井与射孔间的压差。
    Type_B_o BPR;  ///< 体积压力。
    Type_B_o SOIL; ///< 体积油饱和度。
    Type_B_o SGAS; ///< 体积气饱和度。
    Type_B_o SWAT; ///< 体积水饱和度。
};

/// \class BasicGridPropertyParam
/// \brief 基础网格属性参数。
class BasicGridPropertyParam{
public:
    // 下面的注释说明了每个成员变量的基本功能和默认值
    OCP_BOOL PRE{OCP_FALSE};   ///< 网格压力。
    OCP_BOOL PGAS{OCP_FALSE};  ///< 网格天然气压力。
    OCP_BOOL PWAT{OCP_FALSE};  ///< 网格水压力。
    OCP_BOOL SOIL{OCP_FALSE};  ///< 网格油饱和度。
    OCP_BOOL SGAS{OCP_FALSE};  ///< 网格气饱和度。
    OCP_BOOL SWAT{OCP_FALSE};  ///< 网格水饱和度。
    OCP_BOOL DENO{OCP_FALSE};  ///< 网格油密度。
    OCP_BOOL DENG{OCP_FALSE};  ///< 网格气密度。
    OCP_BOOL DENW{OCP_FALSE};  ///< 网格水密度。
    OCP_BOOL KRO{OCP_FALSE};   ///< 网格油相相对渗透率。
    OCP_BOOL KRG{OCP_FALSE};   ///< 网格气相相对渗透率。
    OCP_BOOL KRW{OCP_FALSE};   ///< 网格水相相对渗透率。
    OCP_BOOL BOIL{OCP_FALSE};  ///< 网格油藏摩尔密度。
    OCP_BOOL BGAS{OCP_FALSE};  ///< 网格气藏摩尔密度。
    OCP_BOOL BWAT{OCP_FALSE};  ///< 网格水藏摩尔密度。
    OCP_BOOL VOIL{OCP_FALSE};  ///< 网格油粘度。
    OCP_BOOL VGAS{OCP_FALSE};  ///< 网格气粘度。
    OCP_BOOL VWAT{OCP_FALSE};  ///< 网格水粘度。
    OCP_BOOL XMF{OCP_FALSE};   ///< 液相组分摩尔分数。
    OCP_BOOL YMF{OCP_FALSE};   ///< 气相组分摩尔分数。
    OCP_BOOL PCW{OCP_FALSE};   ///< 毛细管压力差：Po - Pw。
    OCP_BOOL CO2{OCP_FALSE};   ///< CO2浓度。
    OCP_BOOL SATNUM{OCP_FALSE};///< 饱和度区域编号。
};

/// \class OutputRPTParam
/// \brief OutputRPTParam是ParamOutput的一部分，用于控制详细油藏信息的输出。
class OutputRPTParam{
public:
    OCP_BOOL               useRPT{OCP_FALSE}; ///< 是否使用RPT文件输出。
    BasicGridPropertyParam bgp;               ///< 基础网格属性参数。
};

/// \class OutputVTKParam
/// \brief OutputVTKParam用于控制VTK格式输出的参数。
class OutputVTKParam{
public:
    OCP_BOOL               useVTK{OCP_FALSE}; ///< 是否使用VTK文件输出。
    BasicGridPropertyParam bgp;               ///< 基础网格属性参数。
};

/// \class ParamOutput
/// \brief ParamOutput是一个内部结构，用于存储从输入文件中获取的输出信息。
class ParamOutput{
public:
    OutputSummary  summary;     ///< 输出摘要参数。
    OutputRPTParam outRPTParam; ///< 详细报告参数。
    OutputVTKParam outVTKParam; ///< VTK输出参数。

    /// \brief 输入关键字SUMMARY，包含许多子关键字，指示用户感兴趣的结果。
    /// \param ifs 输入文件流。
    void InputSUMMARY(ifstream& ifs);

    /// \brief 输入SUMMARY中的子关键字，这些关键字的内容是字符串形式的。
    /// \param ifs 输入文件流。
    /// \param obj Type_A_o类型的对象。
    void InputType_A(ifstream& ifs, Type_A_o& obj);

    /// \brief 输入SUMMARY中的子关键字，这些关键字的内容是坐标形式的。
    /// \param ifs 输入文件流。
    /// \param obj Type_B_o类型的对象。
    void InputType_B(ifstream& ifs, Type_B_o& obj);

    /// \brief 输入关键字RPTSCHED，它指示了将哪些详细信息输出到RPT文件中。
    /// \param ifs 输入文件流。
    /// \param keyword 关键字。
    void InputRPTSCHED(ifstream& ifs, const string& keyword);
};

#endif /* end if __PARAMOUTPUT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/