/*! \file    PreParamGridWell.hpp
 *  \brief   预处理参数网格和井的类声明
 *  \author  Shizhe Li
 *  \date    Feb/15/2023
 *
 *  \note    在OpenCAEPoroX中使用的参数与SLB的Eclipse大致兼容，
 *           但它有一些自己的规则以便于使用。它是可扩展的和用户友好的。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARAMGRID_HEADER__
#define __PARAMGRID_HEADER__

// 标准头文件
#include <fstream>
#include <iostream>
#include <string>

// OpenCAEPoroX头文件
#include "OCPConst.hpp"
#include "UtilInput.hpp"
#include "ParamWell.hpp"
#include "CornerGrid.hpp"
#include "GmshGrid.hpp"
#include "Output4Vtk.hpp"

using namespace std;

/**
 * \class ActiveGridCheck
 * \brief 检查网格活动性的类
 *
 * 检查网格是否活动，并根据输入文件中的ACTNUM等信息进行处理。
 */
class ActiveGridCheck
{
public:
    OCP_ULL CheckActivity(const OCPModel& Model, const OCP_DBL& ev, const OCP_DBL& ep,
                          const vector<OCP_DBL>& v, const vector<OCP_DBL>& poro);
    OCP_BOOL IfFluid(const OCP_ULL& n, const OCP_DBL& poro);
    void FreeSomeMemory();

public:
    /// 如果ACTNUM全部为1
    OCP_BOOL        allAct{ OCP_FALSE };
    /// 网格活动性输入文件：numGridLocal：0 = 非活动，1 = 活动。
    vector<USI>     ACTNUM;

public:
    /// 网格数量
    OCP_ULL         numGrid;
    /// 活动网格数
    OCP_ULL         activeGridNum;
    /// 从活动网格到所有网格的索引映射
    vector<OCP_ULL> map_Act2All;
    /// 从网格到活动所有网格的映射
    vector<OCP_SLL> map_All2Act;

protected:
    /// 用于检查的模型
    OCPModel        model;
    /// 体积限制
    OCP_DBL         eV;
    /// 孔隙率限制
    OCP_DBL         eP;

protected:
    void CheckActivityIsoT(const vector<OCP_DBL>& v, const vector<OCP_DBL>& poro);
    void CheckActivityT(const vector<OCP_DBL>& v, const vector<OCP_DBL>& poro);
};

/**
 * \class ConnPair
 * \brief 描述邻接单元之间连接的类
 *
 * 连接对表示与邻近单元的交界面积有效区域。
 */
class ConnPair
{
public:
    ConnPair() = default;
    ConnPair(const OCP_ULL& Id, const OCP_INT& Wgt, const ConnDirect& direct,
             const OCP_DBL& AreaB, const OCP_DBL& AreaE)
        : id(Id), wgt(Wgt), direction(direct), areaB(AreaB), areaE(AreaE) {};
    const auto& ID() const { return id; }
    const auto& WGT() const { return wgt; }
    const auto& Direct() const { return direction; }
    const auto& AreaB() const { return areaB; }
    const auto& AreaE() const { return areaE; }
    void SetTransMult(const OCP_DBL& var) { transMult *= var; }
    const auto& TransMult() const { return transMult; }

protected:
    /// 邻近单元的索引
    OCP_ULL    id;
    /// 边的权重
    OCP_INT    wgt;
    /// 连接方向
    ConnDirect direction;
    /// 该单元与邻近单元之间的有效交界面积，本身
    OCP_DBL    areaB;
    /// 该单元与邻近单元之间的有效交界面积，邻近
    OCP_DBL    areaE;
    /// 透过性乘数
    OCP_DBL    transMult{ 1.0 };
};

/**
 * \class InitialReservoir
 * \brief 描述初始储层条件的类
 *
 * 初始储层条件，例如饱和度，压力等。
 */
class InitialReservoir
{
    friend class PreParamGridWell;
    friend class Reservoir;

public:
    /// 检查并填补缺失数据（如果有必要）
    void CheckData(const OCP_USI& numGrid);

protected:
    /// 初始水饱和度
    vector<OCP_DBL> swat;
    /// 每个网格的初始水饱和度并使用饱和度端点缩放
    vector<OCP_DBL> swatInit;
    /// 是否应当缩放Pcow
    OCP_BOOL        scalePcow{ OCP_FALSE };
};

/**
 * \class PreParamGridWell
 * \brief 输入网格信息和井几何信息的类
 *
 * 从输入文件中输入网格信息和井几何信息，生成活动网格与井之间的连接，输出连接到文件。
 */
class PreParamGridWell
{
    friend class Partition;
    friend class Domain;
    friend class Reservoir;

    /////////////////////////////////////////////////////////////////////
    // 输入文件中的参数
    /////////////////////////////////////////////////////////////////////

public:
    /// 输入参数文件
    void InputFile(const string& myFile, const string& myWorkdir);

protected:
    void Input(const string& myFilename);
    /// 检查输入
    void CheckInput();
    /// 输入参数后处理
    void PostProcessInput();
    /// 输入问题模型：等温或热模型
    void InputMODEL(ifstream& ifs);
    /// 输入DUALPORO
    void InputDUALPORO() { DUALPORO = OCP_TRUE; }
    /// 输入DPGRID
    void InputDPGRID() { DPGRID = OCP_TRUE; }
    /// 输入尺寸
    void InputDIMENS(ifstream& ifs);
    /// 输入EQUALS
    void InputEQUALS(ifstream& ifs);
    /// 输入COPY
    void InputCOPY(ifstream& ifs);
    /// 输入MULTIPLY
    void InputMULTIPLY(ifstream& ifs);
    /// 输入网格属性
    void InputGrid(ifstream& ifs, string& keyword);
    /// 输入WELSPECS
    void InputWELSPECS(ifstream& ifs);
    /// 输入COMPDAT
    void InputCOMPDAT(ifstream& ifs);
    /// 输入INCLUDE
    void InputINCLUDE(ifstream& ifs);

#ifdef WITH_GMSH
    /// 输入GMSH
    void InputGMSH(ifstream& ifs);
    /// 输入物理属性
    void InputGMSHPRO(ifstream& ifs);
#endif

    // 输入工具
    /// 查找指定变量的指针。
    vector<OCP_DBL>* FindPtr(const string& varName, const OCP_DBL&);
    vector<USI>* FindPtr(const string& varName, const USI&);
    /// 在InputEQUALS中使用，批量赋值。
    template <typename T>
    void setVal(vector<T>& obj, const T& val, const vector<USI>& index);
    /// 在InputCOPY中使用，将一个变量的值复制到另一个变量。
    template <typename T>
    void CopyVal(vector<T>& obj, const vector<T>& src, const vector<USI>& index);
    /// 在InputMULTIPLY中使用，将变量的一定范围的值乘以一个系数。
    void MultiplyVal(vector<OCP_DBL>& obj, const OCP_DBL& val, const vector<USI>& index);

protected:
    /// 当前工作目录
    string          workdir;
    /// 正交网格或角点网格
    GridType        gridType{ GridType::orthogonal };
    /// 网格数量（numGridM + numGridF）
    OCP_ULL         numGrid{ 0 };
    /// 矩阵网格数量
    OCP_ULL         numGridM{ 0 };
    /// 裂缝网格数量
    OCP_ULL         numGridF{ 0 };
    /// 热模型或等温模型
    OCPModel        model{ OCPModel::none };
    // 正交网格
    /// x方向上的网格数量。
    USI             nx;
    /// y方向上的网格数量。
    USI             ny;
    /// z方向上的网格数量。
    USI             nz;
    /// 每个网格的x方向尺寸。
    vector<OCP_DBL> dx;
    /// 每个网格的y方向尺寸。
    vector<OCP_DBL> dy;
    /// 每个网格的z方向尺寸。
    vector<OCP_DBL> dz;
    /// 最上层网格顶面的深度。
    vector<OCP_DBL> tops;
    // 角点网格
    /// 角点网格的线。
    vector<OCP_DBL> coord;
    /// 角点网格的Z值。
    vector<OCP_DBL> zcorn;
    // 通用网格
#ifdef WITH_GMSH
    /// GMSH网格
    GMSHGrid        gmshGrid;
#endif
    // 岩石参数
    /// 净毛比
    vector<OCP_DBL> ntg;
    /// 孔隙率
    vector<OCP_DBL> poro;
    /// x方向渗透率
    vector<OCP_DBL> kx;
    /// y方向渗透率
    vector<OCP_DBL> ky;
    /// z方向渗透率
    vector<OCP_DBL> kz;
    /// Z方向透过性乘数
    vector<OCP_DBL> multZ;
    /// 记录每个网格的SAT区域索引。
    vector<USI>     SATNUM;
    /// 记录每个网格的PVT区域索引。
    vector<USI>     PVTNUM;
    /// 记录每个网格的ROCK区域索引。
    vector<USI>     ROCKNUM;
    /// 边界指示器
    vector<USI>     boundIndex;
    /// 边界面积
    vector<OCP_DBL> boundArea;
    // 初始条件
    InitialReservoir initR;
    /// 网格活动性检查
    ActiveGridCheck  actGC;
    // 井
    vector<WellParam> well;

    /////////////////////////////////////////////////////////////////////
    // 生成网格连接
    /////////////////////////////////////////////////////////////////////

public:
    void Setup();

protected:
    /// 设置网格
    void SetupGrid();
    /// 设置正交网格。
    void SetupOrthogonalGrid();
    /// 计算正交网格的深度和体积。
    void CalDepthVOrthogonalGrid();
    /// 设置正交网格的邻接信息。
    void SetupActiveConnOrthogonalGrid();
    /// 设置正交网格单介质的连接。
    void SetupActiveConnOrthogonalGridSM();
    /// 设置正交网格双孔隙的连接。
    void SetupActiveConnOrthogonalGridDP();
    /// 输出正交网格的点。
    void OutputPointsOrthogonalGrid();
    /// 设置角点网格。
    void SetupCornerGrid();
    /// 设置角点网格的dx,dy,dz,depth,v。
    void SetupBasicCornerGrid(const OCP_COORD& CoTmp);
    /// 设置角点网格的邻接信息。
    void SetupActiveConnCornerGrid(const OCP_COORD& CoTmp);
    /// 设置角点网格单介质的连接。
    void SetupActiveConnCornerGridSM(const OCP_COORD& CoTmp);
    /// 设置角点网格双孔隙的连接。
    void SetupActiveConnCornerGridDP(const OCP_COORD& CoTmp);
    /// 输出角点网格的点。
    void OutputPointsCornerGrid(const OCP_COORD& mycord);

#ifdef WITH_GMSH
    /// 设置GMSH网格。
    void SetupGmshGrid();
    /// 设置GMSH网格的深度和体积。
    void SetupBasicGmshGrid();
    /// 设置GMSH网格的邻接信息。
    void SetupActiveConnGmshGrid();
    /// 输出GMSH网格的点。
    void OutputPointsGmshGrid();
#endif

    // 结构网格
    /// 设置网格位置：顶部，底部，侧面或内部
    void SetLocationStructral();

    /// 通用网格
    /// 计算网格单元的活动性
    void CalActiveGrid(const OCP_DBL& ev, const OCP_DBL& ep);
    /// 设置透过性乘数
    void SetupTransMult();

protected:
    /// 单元体积
    vector<OCP_DBL>          v;
    /// 网格单元中心的深度
    vector<OCP_DBL>          depth;
    // 连接
    /// 活动网格的邻接信息。
    vector<vector<ConnPair>> gNeighbor;
    /// 活动网格数量。
    OCP_ULL                  activeGridNum;

    /////////////////////////////////////////////////////////////////////
    // 双孔隙选项
    /////////////////////////////////////////////////////////////////////

protected:
    /// 是否使用双孔隙选项。（矩阵在前面，裂缝网格在后面）
    OCP_BOOL        DUALPORO{ OCP_FALSE };
    /// 如果为true，则裂缝的默认属性将会从矩阵复制（考虑一下）
    OCP_BOOL        DPGRID{ OCP_FALSE };
    /// 双孔隙矩阵-裂缝耦合项中使用的sigma因子
    vector<OCP_DBL> sigma;
    /// 矩阵材料块的垂直尺寸
    vector<OCP_DBL> dzMtrx;

    /////////////////////////////////////////////////////////////////////
    // 生成活动网格和井之间的连接
    /////////////////////////////////////////////////////////////////////

protected:
    /// 设置井和活动网格之间的连接
    void SetupConnWellGrid();
    /// 返回井穿透的块体索引
    OCP_ULL GetPerfLocation(const WellParam& well, const USI& p);

protected:
    /// 井数量
    USI                     numWell;
    /// 井和活动网格之间的连接
    vector<vector<OCP_ULL>> connWellGrid;

    /////////////////////////////////////////////////////////////////////
    // 连接权重
    /////////////////////////////////////////////////////////////////////

protected:
    const OCP_INT WEIGHT_GG = 1;         ///< 网格-网格
    const OCP_INT WEIGHT_GW = 1000000;   ///< 网格-井

    /////////////////////////////////////////////////////////////////////
    // 输出基本网格信息和活动网格连接
    /////////////////////////////////////////////////////////////////////

public:
    void OutputBaiscInfo() const;   ///< 计算并返回网格的基本信息

protected:
    OCP_BOOL ifUseVtk{ OCP_FALSE };

public:
    void FreeMemory();
};

#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Feb/15/2023      Create file                          */
/*----------------------------------------------------------------------------*/