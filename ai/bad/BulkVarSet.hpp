```cpp
/*! 
 *  \file    BulkVarSet.hpp
 *  \brief   此文件包含BulkVarSet类的声明
 *  \author  Shizhe Li
 *  \date    Aug/19/2023
 *
 *  BulkVarSet类用于存储和管理与油藏模拟中的岩石和流体体积单元相关的各种变量。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BULKVARSET_HEADER__
#define __BULKVARSET_HEADER__

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include <vector>

using namespace std;

/*!
 * \enum BulkContent
 * \brief 描述体积单元的内容类型
 */
enum class BulkContent : USI {
    r,      ///< 都是岩石
    rf      ///< 岩石和流体
};

/*!
 * \class BulkVarSet
 * \brief 体积单元的基本变量集
 *
 * BulkVarSet类存储了用于基础模拟的所有体积单元内的变量
 */
class BulkVarSet {
public:
    OCP_USI nb;      ///< 体积单元的数量
    OCP_USI nbI;     ///< 内部体积单元的数量
    USI np;          ///< 相的数量
    USI nc;          ///< 组分的数量

public:
    vector<OCP_DBL> dx;      ///< x方向上的单元大小
    vector<OCP_DBL> dy;      ///< y方向上的单元大小
    vector<OCP_DBL> dz;      ///< z方向上的单元大小
    vector<OCP_DBL> v;       ///< 体积单元的体积
    vector<OCP_DBL> depth;   ///< 体积单元中心的深度
    vector<OCP_DBL> ntg;     ///< 体积单元的净毛比
    vector<OCP_DBL> poroInit;///< 初始岩石孔隙度(*ntg)
    vector<OCP_DBL> poro;    ///< 岩石孔隙度(*ntg)
    vector<OCP_DBL> rockVp;  ///< 岩石孔隙体积 = Vgrid * poro
    vector<OCP_DBL> rockKx;  ///< x方向上的岩石渗透率
    vector<OCP_DBL> rockKy;  ///< y方向上的岩石渗透率
    vector<OCP_DBL> rockKz;  ///< z方向上的岩石渗透率
    vector<OCP_DBL> sigma;   ///< 双重孔隙度矩阵-裂缝耦合项中使用的sigma因子
    vector<OCP_DBL> dzMtrx;  ///< 矩阵材料块的垂直尺寸
    vector<OCP_DBL> vr;      ///< 岩石的体积
    vector<OCP_DBL> Hr;      ///< 岩石的焓

    vector<OCP_DBL>    lporo;   ///< 上一个时间步的岩石孔隙度    
    vector<OCP_DBL>    lrockVp; ///< 上一个时间步的岩石孔隙体积         
    vector<OCP_DBL>    lvr;     ///< 上一个时间步的岩石体积       
    vector<OCP_DBL>    lHr;     ///< 上一个时间步的岩石的焓

    
    vector<OCP_DBL>    poroP; ///< 岩石孔隙度对压力的导数
    vector<OCP_DBL>    poroT; ///< 岩石孔隙度对温度的导数
    vector<OCP_DBL>    vrP;   ///< 岩石体积对压力的导数
    vector<OCP_DBL>    vrT;   ///< 岩石体积对温度的导数
    vector<OCP_DBL>    HrT;   ///< 岩石焓对温度的导数

    vector<OCP_DBL>    lporoP; ///< 上一个时间步的岩石孔隙度对压力的导数
    vector<OCP_DBL>    lporoT; ///< 上一个时间步的岩石孔隙度对温度的导数
    vector<OCP_DBL>    lvrP;   ///< 上一个时间步的岩石体积对压力的导数
    vector<OCP_DBL>    lvrT;   ///< 上一个时间步的岩石体积对温度的导数
    vector<OCP_DBL>    lHrT;   ///< 上一个时间步的岩石焓对温度的导数


public:

    
    INT                 o, g, w;     ///< 油气水相的索引，负表示不存在
    INT                 r;           ///< 润湿相的索引
    vector<OCP_DBL>     T;           ///< 温度
    vector<OCP_DBL>     P;           ///< 压力
    vector<OCP_DBL>     Pb;          ///< 泡点压力
    vector<OCP_DBL>     vf;          ///< 流体体积
    vector<OCP_DBL>     Nt;          ///< 组分总数 
    vector<OCP_DBL>     Ni;          ///< 各组分的数量
    vector<OCP_BOOL>    phaseExist;  ///< 相的存在性
    vector<OCP_DBL>     S;           ///< 相的饱和度
    vector<OCP_DBL>     vj;          ///< 相的体积
    vector<OCP_DBL>     xij;         ///< 各相中各组分的摩尔/质量占比
    vector<OCP_DBL>     xi;          ///< 相的摩尔/质量浓度
    vector<OCP_DBL>     rho;         ///< 相的密度
    vector<OCP_DBL>     mu;          ///< 相的粘度
    vector<OCP_DBL>     kr;          ///< 相的相对渗透率
    vector<OCP_DBL>     Pc;          ///< 相的毛管力
    vector<OCP_DBL>     Pj;          ///< 各相的压力
    vector<OCP_DBL>     Uf;          ///< 流体的内能
    vector<OCP_DBL>     H;           ///< 各相的焓

    vector<OCP_DBL>     lT;           ///< 上一个时间步的温度
    vector<OCP_DBL>     lP;           ///< 上一个时间步的压力
    vector<OCP_DBL>     lPb;          ///< 上一个时间步的泡点压力
    vector<OCP_DBL>     lvf;          ///< 上一个时间步的流体体积
    vector<OCP_DBL>     lNt;          ///< 上一个时间步的组分总数 
    vector<OCP_DBL>     lNi;          ///< 上一个时间步的各组分的数量
    vector<OCP_BOOL>    lphaseExist;  ///< 上一个时间步的相的存在性
    vector<OCP_DBL>     lS;           ///< 上一个时间步的相的饱和度
    vector<OCP_DBL>     lvj;          ///< 上一个时间步的相的体积
    vector<OCP_DBL>     lxij;         ///< 上一个时间步的各相中各组分的摩尔/质量占比
    vector<OCP_DBL>     lxi;          ///< 上一个时间步的相的摩尔/质量浓度
    vector<OCP_DBL>     lrho;         ///< 上一个时间步的相的密度
    vector<OCP_DBL>     lmu;          ///< 上一个时间步的相的粘度
    vector<OCP_DBL>     lkr;          ///< 上一个时间步的相的相对渗透率
    vector<OCP_DBL>     lPc;          ///< 上一个时间步的相的毛管力
    vector<OCP_DBL>     lPj;          ///< 上一个时间步的各相的压力
    vector<OCP_DBL>     lUf;          ///< 上一个时间步的流体的内能
    vector<OCP_DBL>     lH;           ///< 上一个时间步的各相的焓

    vector<OCP_DBL>     vfP;          ///< 流体体积对压力的导数
    vector<OCP_DBL>     vfT;          ///< 流体体积对温度的导数
    vector<OCP_DBL>     vfi;          ///< 流体体积对各组分数的导数
    vector<OCP_DBL>     xiP;          ///< 流体摩尔浓度对压力的导数
    vector<OCP_DBL>     xiT;          ///< 流体摩尔浓度对温度的导数
    vector<OCP_DBL>     xix;          ///< 流体摩尔浓度对各组分占比的导数
    vector<OCP_DBL>     rhoP;         ///< 流体密度对压力的导数
    vector<OCP_DBL>     rhoT;         ///< 流体密度对温度的导数
    vector<OCP_DBL>     rhox;         ///< 流体密度对各组分占比的导数
    vector<OCP_DBL>     muP;          ///< 流体粘度对压力的导数
    vector<OCP_DBL>     muT;          ///< 流体粘度对温度的导数
    vector<OCP_DBL>     mux;          ///< 流体粘度对各组分占比的导数
    vector<OCP_DBL>     dPcdS;        ///< 相毛管力对饱和度的导数
    vector<OCP_DBL>     dKrdS;        ///< 相相对渗透率对饱和度的导数
    vector<OCP_DBL>     UfP;          ///< 流体内能对压力的导数 
    vector<OCP_DBL>     UfT;          ///< 流体内能对温度的导数
    vector<OCP_DBL>     Ufi;          ///< 流体内能对各组分数的导数
    vector<OCP_DBL>     HT;           ///< 相焓对温度的导数
    vector<OCP_DBL>     Hx;           ///< 相焓对各组分占比的导数

    vector<OCP_DBL>     lvfP;         ///< 上一时间步的流体体积对压力的导数
    vector<OCP_DBL>     lvfT;         ///< 上一时间步的流体体积对温度的导数
    vector<OCP_DBL>     lvfi;         ///< 上一时间步的流体体积对各组分数的导数
    vector<OCP_DBL>     lxiP;         ///< 上一时间步的流体摩尔浓度对压力的导数
    vector<OCP_DBL>     lxiT;         ///< 上一时间步的流体摩尔浓度对温度的导数
    vector<OCP_DBL>     lxix;         ///< 上一时间步的流体摩尔浓度对各组分占比的导数
    vector<OCP_DBL>     lrhoP;        ///< 上一时间步的流体密度对压力的导数
    vector<OCP_DBL>     lrhoT;        ///< 上一时间步的流体密度对温度的导数
    vector<OCP_DBL>     lrhox;        ///< 上一时间步的流体密度对各组分占比的导数
    vector<OCP_DBL>     lmuP;         ///< 上一时间步的流体粘度对压力的导数
    vector<OCP_DBL>     lmuT;         ///< 上一时间步的流体粘度对温度的导数
    vector<OCP_DBL>     lmux;         ///< 上一时间步的流体粘度对各组分占比的导数
    vector<OCP_DBL>     ldPcdS;       ///< 上一时间步的相毛管力对饱和度的导数
    vector<OCP_DBL>     ldKrdS;       ///< 上一时间步的相相对渗透率对饱和度的导数
    vector<OCP_DBL>     lUfP;         ///< 上一时间步的流体内能对压力的导数 
    vector<OCP_DBL>     lUfT;         ///< 上一时间步的流体内能对温度的导数
    vector<OCP_DBL>     lUfi;         ///< 上一时间步的流体内能对各组分数的导数
    vector<OCP_DBL>     lHT;          ///< 上一时间步的相焓对温度的导数
    vector<OCP_DBL>     lHx;          ///< 上一时间步的相焓对各组分占比的导数

    
    USI                 lendSdP;      ///< 变量dSec_dPri的长度
    vector<OCP_DBL>     dSec_dPri;    ///< 次要变量对主要变量的导数
    vector<OCP_DBL>     ldSec_dPri;   ///< 上一时间步的次要变量对主要变量的导数

    vector<BulkContent> cType;        ///< 网格块包含物的类型
    vector<OCP_DBL>     initT;        ///< 网格块初始温度
};

#endif /* end if __BulkVarSet_HEADER__ */

/*----------------------------------------------------------------------------
 *  Brief Change History of This File
 *----------------------------------------------------------------------------
 *
 *  Author              Date             Actions
 *----------------------------------------------------------------------------
 *
 *  Shizhe Li           Aug/19/2023      Create file
 *
 *----------------------------------------------------------------------------
 */
```

