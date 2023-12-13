/*! \file    ParamReservoir.hpp
 *  \brief   油藏参数输入管理模块
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARAMRESERVOIR_HEADER__
#define __PARAMRESERVOIR_HEADER__

 // Standard header files
#include <fstream>
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "UtilInput.hpp"
#include "UtilOutput.hpp"

using namespace std;

/// 一维表格型参数数据结构
class TableSet
{
public:
    /// 打印表格到屏幕
    void DisplayTable() const;

public:
    string                          name;       ///< 表格名字
    USI                             colNum;     ///< 名字
    vector<vector<vector<OCP_DBL>>> data;       ///< 名字
};

/// 二维表格型参数数据结构
class Table2
{
public:
    /// 构造函数
    Table2(const USI& n) { data.resize(n); }
    /// 设置表格的列数
    void SetColNum() { colNum = data[0].size(); }
public:

    string                          refName;       ///< 第一维变量名字
    vector<OCP_DBL>                 refData;       ///< 第一维变量数据
    USI                             colNum;        ///< 第二维列数
    vector<vector<vector<OCP_DBL>>> data;          ///< 表格数据
};


/// 二维表格序列
class Table2Set
{
public:
    string           name;   ///< 表格名字
    vector<Table2>   data;   ///< 表格序列
};

/// 热损失参数
class HLoss
{
public:
    OCP_BOOL ifHLoss{ OCP_FALSE }; ///< 是否使用
    OCP_BOOL obUse{ OCP_FALSE };   ///< 是否使用覆层热损失
    OCP_DBL  obC{ -1 };            ///< 覆层热损失参数
    OCP_DBL  obK{ -1 };            ///< 覆层热损失参数 
    OCP_BOOL ubUse{ OCP_FALSE };   ///< 是否使用低层热损失 
    OCP_DBL  ubC{ -1 };            ///< 底层热损失参数 
    OCP_DBL  ubK{ -1 };            ///< 底层热损失参数 
};

/// 岩性参数
class RockParam
{
public:
    string   type{ "LINEAR" };      ///< 压缩类型
    OCP_DBL  Pref{ 14.7 };          ///< 参考压力
    OCP_DBL  Tref{ 60 };            ///< 参考温度
    OCP_DBL  cp1{ 3.406E-6 };       ///< 一阶压缩系数
    OCP_DBL  cp2{ 0 };              ///< 二阶压缩系数
    OCP_DBL  ct{ 0 };               ///< 热膨胀系数
    OCP_DBL  cpt{ 0 };              ///< 热膨胀系数
    OCP_BOOL ConstRock{ OCP_TRUE }; ///< 岩石体积是否为常数
    OCP_DBL HCP1{ 35 };             ///< 岩石焓系数
    OCP_DBL HCP2{ 0 };              ///< 岩石焓系数
};


/// 初始油藏平衡计算
class EQUILParam
{
public:
    vector<OCP_DBL> data; ///< 平衡参数
};


/// Brooks-Corey 相对渗透率毛管力模型
class BrooksCoreyParam
{
public:
    OCP_DBL sw_imm;     ///< 不可移动的润湿相饱和度
    OCP_DBL sn_imm;     ///< 不可移动的非润湿相饱和度
    OCP_DBL Pentry;     ///< 气相进入压力
    OCP_DBL Pcmax;      ///< 最大气水毛管力
    OCP_DBL Cw_kr;      ///< 润湿相相对渗透率指数
    OCP_DBL Cn_kr;      ///< 非润湿相相对渗透率指数
    OCP_DBL C_pc;       ///< 毛管力指数
};


/// 边界条件参数
class BoundaryParam
{
public:
    /// 构造函数
    BoundaryParam(const string& Name) : name(Name) {}

    string   name;                     ///< 边界名字
    OCP_BOOL constP{ OCP_FALSE };      ///< 是否常压控制
    OCP_DBL  P;                        ///< 常压控制的压力值
};


/// 参数结构1
template <typename T>
class Type_A_r
{
public:
    OCP_BOOL  activity{ OCP_FALSE }; ///< 是否使用
    vector<T> data;                  ///< 参数数据
};

/// 组分类参数
class ComponentParam
{
public:
    // Basic params
    /// 初始化参数
    void Init();
    /// 读入组分信息
    void InputCOMPONENTS(ifstream& ifs, const string& keyword);
    /// 组分变量指针查找器1
    Type_A_r<vector<OCP_DBL>>* FindPtr01(const string& varName);
    /// 读入参考压力和参考温度
    void InputRefPR(ifstream& ifs, const string& keyword);
    /// 组分变量指针查找器2
    vector<OCP_DBL>* FindPtr02(const string& varName);
    /// 读入烃组分名字
    void InputCNAMES(ifstream& ifs);
    /// 读入LBC系数
    void InputLBCCOEF(ifstream& ifs);
    /// 读入BIC系数
    void InputBIC(ifstream& ifs);
    /// 读入稳定性分析的SSM参数
    void InputSSMSTA(ifstream& ifs);
    /// 读入稳定性分析的NR参数
    void InputNRSTA(ifstream& ifs);
    /// 读入相分裂计算的SSM参数
    void InputSSMSP(ifstream& ifs);
    /// 读入相分裂计算的NR参数
    void InputNRSP(ifstream& ifs);
    /// 读入RR方程的求解参数
    void InputRR(ifstream& ifs);

public:
    USI            NTPVT;               ///< PVT区域
    USI            numCom{ 0 };         ///< 相平衡计算的组分数
    USI            numPhase{ 2 };       ///< 相平衡计算的最大相数
    vector<string> Cname;               ///< 组分的名字
    Type_A_r<vector<OCP_DBL>> Tc;       ///< 组分临界温度
    Type_A_r<vector<OCP_DBL>> Pc;       ///< 组分临界压力
    Type_A_r<vector<OCP_DBL>> Vc;       ///< 组分临界体积
    Type_A_r<vector<OCP_DBL>> Zc;       ///< 组分临界压缩因子
    Type_A_r<vector<OCP_DBL>> MW;       ///< 组分分子质量
    Type_A_r<vector<OCP_DBL>> Acf;      ///< 组分偏心因子
    Type_A_r<vector<OCP_DBL>> OmegaA;   ///< PR-EoS系数
    Type_A_r<vector<OCP_DBL>> OmegaB;   ///< PR-EoS系数
    Type_A_r<vector<OCP_DBL>> Vshift;   ///< 组分体积变换系数
    Type_A_r<vector<OCP_DBL>> parachor; ///< parachor

    // for viscosity calculation
    Type_A_r<vector<OCP_DBL>> Vcvis;    ///< 粘度计算临界体积
    Type_A_r<vector<OCP_DBL>> Zcvis;    ///< 粘度计算压缩因子
    vector<OCP_DBL>           LBCcoef;  ///< LBC粘度计算系数
    vector<vector<OCP_DBL>>   BIC;      ///< BIC参数

    // Thermal only
    Type_A_r<vector<OCP_DBL>> molden;   ///< 参考温度压力下的摩尔浓度
    Type_A_r<vector<OCP_DBL>> cp;       ///< 组分压缩系数
    Type_A_r<vector<OCP_DBL>> ct1;      ///< 第一热膨胀系数
    Type_A_r<vector<OCP_DBL>> ct2;      ///< 第一热膨胀系数
    Type_A_r<vector<OCP_DBL>> cpt;      ///< 温度压力依赖系数
    Type_A_r<vector<OCP_DBL>> cpl1;     ///< 液相焓计算系数
    Type_A_r<vector<OCP_DBL>> cpl2;     ///< 液相焓计算系数
    Type_A_r<vector<OCP_DBL>> cpl3;     ///< 液相焓计算系数
    Type_A_r<vector<OCP_DBL>> cpl4;     ///< 液相焓计算系数
    Type_A_r<vector<OCP_DBL>> cpg1;     ///< 气相焓计算系数
    Type_A_r<vector<OCP_DBL>> cpg2;     ///< 气相焓计算系数
    Type_A_r<vector<OCP_DBL>> cpg3;     ///< 气相焓计算系数
    Type_A_r<vector<OCP_DBL>> cpg4;     ///< 气相焓计算系数
    Type_A_r<vector<OCP_DBL>> hvapr;    ///< 蒸发焓计算系数
    Type_A_r<vector<OCP_DBL>> hvr;      ///< 蒸发焓计算系数
    Type_A_r<vector<OCP_DBL>> ev;       ///< 蒸发焓计算系数
    Type_A_r<vector<OCP_DBL>> avisc;    ///< 液相组分的粘度关联参数 
    Type_A_r<vector<OCP_DBL>> bvisc;    ///< 液相组分的粘度关联参数 
    Type_A_r<vector<OCP_DBL>> avg;      ///< 气相组分的粘度关联参数 
    Type_A_r<vector<OCP_DBL>> bvg;      ///< 气相组分的粘度关联参数 


    Table2Set viscTab;                  ///< 粘度温度压力依赖表

    vector<OCP_DBL> Pref;               ///< 参考压力
    vector<OCP_DBL> Tref;               ///< 参考温度
    vector<string> SSMparamSTA;         ///< 稳定性分析SSM参数
    vector<string> NRparamSTA;          ///< 稳定性分析NR参数
    vector<string> SSMparamSP;          ///< 相分裂计算SSM参数
    vector<string> NRparamSP;           ///< 相分裂计算NR参数
    vector<string> RRparam;             ///< RR求解参数
};


/// 混溶控制参数
class Miscstr
{
public:
    OCP_BOOL        ifMiscible{ OCP_FALSE };///< 是否使用混溶
    OCP_DBL         surTenRef{ -1 };        ///< 参考表面张量
    OCP_DBL         surTenEpt{ -1 };        ///< 表面张量相关参数
    OCP_DBL         surTenPc{ -1 };         ///< 表面张量相关参数
    OCP_DBL         surTenExp{ 0.25 };      ///< 表面张量相关参数
};


/// 油藏输入参数模块
class ParamReservoir
{

public:

    string                   unitType;    ///< 单位类型
                                          
    OCP_DBL                  rsTemp;      ///< 油藏初始温度
    vector<RockParam>        rockSet;     ///< 岩石参数集
    HLoss                    hLoss;       ///< 热损失参数集
    Miscstr                  miscstr;     ///< 混溶控制参数集
    vector<BrooksCoreyParam> BCparam;     ///< Brooks-Corey模型参数集
    vector<BoundaryParam>    BDparam;     ///< 边界条件参数集

    Type_A_r<OCP_DBL> density;             ///< 标态下各相密度
    Type_A_r<OCP_DBL> gravity;             ///< 标态下各相重力系数
    OCP_BOOL          ifThcon{ OCP_FALSE };///< 是否使用热传导
    OCP_DBL           thcono{ 24 };        ///< 油相热传导系数
    OCP_DBL           thcong{ 24 };        ///< 气相热传导系数
    OCP_DBL           thconw{ 24 };        ///< 水相热传导系数
    OCP_DBL           thconr{ 24 };        ///< 岩石热传导系数

    // mixture Models
    OCP_BOOL blackOil{ OCP_FALSE }; ///< 是否使用黑油模型
    OCP_BOOL comps{ OCP_FALSE };    ///< 是否使用组分模型
    OCP_BOOL thermal{ OCP_FALSE };  ///< 是否使用热流模型
    OCP_BOOL oil{ OCP_FALSE };      ///< 油相是否存在
    OCP_BOOL gas{ OCP_FALSE };      ///< 气相是否存在
    OCP_BOOL water{ OCP_FALSE };    ///< 水相是否存在
    OCP_BOOL disGas{ OCP_FALSE };   ///< 油相溶解气是否存在

    // flow model
    OCP_BOOL GRAVDR{ OCP_FALSE };   ///< 是否要在双孔模型中使用重力迟滞

    ComponentParam comsParam;       ///< 组分变量参数集

    // SAT Region & PVT Region
    USI               NTSFUN{ 1 }; ///< 饱和度区域数
    USI               NTPVT{ 1 };  ///< PVT区域
    USI               NTROOC{ 1 }; ///< 岩石区域

    TableSet SWFN_T;               ///< SWFN表格
    TableSet SWOF_T;               ///< SWOF表格
    TableSet SGFN_T;               ///< SGFN表格
    TableSet SGOF_T;               ///< SGOF表格
    TableSet SOF3_T;               ///< SOF3表格
    TableSet PBVD_T;               ///< PBVD表格
    // initial zi vs depth
    TableSet           ZMFVD_T;  ///< ZMFVD表格
    TableSet           TEMPVD_T; ///< TEMPVD表格
    vector<EQUILParam> EQUIL;    ///< 平衡条件参数

    // PVT properties
    USI numPhase;     ///< 相数
    USI numCom;       ///< 组分数
    TableSet PVCO_T;  ///< PVCO表格
    TableSet PVDO_T;  ///< PVDO表格
    TableSet PVCDO_T; ///< PVCDO表格
    TableSet PVDG_T;  ///< PVDG表格
    TableSet PVTW_T;  ///< PVTW表格


    OCP_BOOL  GARCIAW{ OCP_FALSE };   ///< 是否使用GARCIAW方法
    Table2Set PVTH2O;                 ///< PVTH2O表格
    Table2Set PVTCO2;                 ///< PVTCO2表格
    OCP_DBL   Psurf;                  ///< 地表压力
    OCP_DBL   Tsurf;                  ///< 地表温度

public:

    /// 表格寻址函数1
    TableSet* FindPtrTable(const string& varName);
    /// 表格寻址函数2
    Table2Set* FindPtrTable2(const string& varName);

    /// 部分油藏参数初始化
    void Init();

    /// 油藏表格初始化
    void InitTable();

    /// 输入组分参数
    void InputCOMPS(ifstream& ifs);

    /// 输入油藏初始温度
    void InputRTEMP(ifstream& ifs);

    /// 输入一维表格
    void InputTABLE(ifstream& ifs, const string& tabName);
    /// 输入二维表格
    void InputTABLE2(ifstream& ifs, const string& tabName);

    /// 输入等温岩石参数
    void InputROCK(ifstream& ifs);
    /// 输入非等温岩石参数
    void InputROCKT(ifstream& ifs);
    /// 输入热损失参数
    void InputHLOSS(ifstream& ifs);
    /// 输入Brooks-Corey模型参数
    void InputBrooksCorey(ifstream& ifs);

    /// 输入混溶控制参数
    void InputMISCSTR(ifstream& ifs);

    /// 输入重力参数
    void InputGRAVITY(ifstream& ifs);

    /// 输入密度参数
    void InputDENSITY(ifstream& ifs);

    /// 输入热传导参数
    void InputTHCON(ifstream& ifs, const string& keyword);

    /// 输入平衡初始化参数
    void InputEQUIL(ifstream& ifs);

    /// 输入区域数
    void InputTABDIMS(ifstream& ifs);

    /// 输入组分数
    void InputNCOMPS(ifstream& ifs) {
        vector<string> vbuf;
        ReadLine(ifs, vbuf);
        comsParam.numCom = stoi(vbuf[0]);
        numCom = comsParam.numCom;
    }
    /// 输入组分名字
    void InputCNAMES(ifstream& ifs) { comsParam.InputCNAMES(ifs); };
    /// 输入组分参数
    void InputCOMPONENTS(ifstream& ifs, const string& keyword)
    {
        comsParam.InputCOMPONENTS(ifs, keyword);
    }
    /// 输入LBC粘度参数
    void InputLBCCOEF(ifstream& ifs) { comsParam.InputLBCCOEF(ifs); }
    /// 输入BIC参数
    void InputBIC(ifstream& ifs) { comsParam.InputBIC(ifs); };
    /// 输入组分参考压力
    void InputRefPR(ifstream& ifs, const string& keyword)
    {
        comsParam.InputRefPR(ifs, keyword);
    };

    /// 输入稳定性分析SSM参数
    void InputSSMSTA(ifstream& ifs) { comsParam.InputSSMSTA(ifs); };
    /// 输入稳定性分析NR参数
    void InputNRSTA(ifstream& ifs) { comsParam.InputNRSTA(ifs); };
    /// 输入相分裂计算SSM参数
    void InputSSMSP(ifstream& ifs) { comsParam.InputSSMSP(ifs); };
    /// 输入相分裂计算NR参数
    void InputNRSP(ifstream& ifs) { comsParam.InputNRSP(ifs); };
    /// 输入RR方程求解参数
    void InputRR(ifstream& ifs) { comsParam.InputRR(ifs); };

    /// 输入边界参数
    void InputBoundary(ifstream& ifs);

    /// 检查参数输入正确性
    void CheckParam();

    /// 检查岩石参数输入正确性
    void CheckRock();

    /// 检查CPL参数输入正确性
    void CheckCPL();

    /// 检查CPL参数输入正确性
    void CheckCPG();
};

#endif /* end if __PARAMRESERVOIR_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/09/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/