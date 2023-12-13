/*! \file    CornerGrid.hpp
 *  \brief   声明与角点网格相关的类
 *  \author  Shizhe Li
 *  \date    Nov/16/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __CORNERGRID_HEADER__
#define __CORNERGRID_HEADER__

// OpenCAEPoroX头文件
#include "UtilMesh.hpp"
using namespace std;

// 角点网格代码使用的常量
const OCP_DBL TEENY        = 1E-3;  ///< 用于检查中心到面的距离
const USI     MAX_NEIGHBOR = 80;    ///< 允许的最大邻居数量

/// 半连接类
class HalfConn{
public:
    OCP_DBL    Ad_dd;        ///< 连接的导数
    Point3D    d;            ///< 连接方向的单位向量
    OCP_ULL    neigh;        ///< 邻居的索引
    ConnDirect directionType;///< 连接的方向类型
};

/// 所有连接的类
class ConnGrid{
public:
    USI              nConn, maxConn; ///< 现有连接数和最大连接数
    vector<HalfConn> halfConn;      ///< 半连接的向量

    /// 分配内存给halfConn向量
    void Allocate(const USI& max_neighbor);

    /// 添加一个半连接到halfConn向量
    void AddHalfConn(const OCP_ULL& n,
                     const Point3D& area,
                     const Point3D& d,
                     const ConnDirect& direction,
                     const OCP_DBL& flag = 1);
};

/// 通用连接类
class GeneralConnect{
public:
    OCP_ULL    begin, end;       ///< 连接的起始和结束索引
    ConnDirect directionType;    ///< 连接的方向类型
    OCP_DBL    Ad_dd_begin;      ///< 起始点的连接导数
    OCP_DBL    Ad_dd_end;        ///< 结束点的连接导数
};

/// 角点网格类
class OCP_COORD{
    friend class PreParamGridWell; ///< 预处理网格井参数的友元类
public:
    /// 分配内存给网格数据
    void Allocate(const USI& Nx, const USI& Ny, const USI& Nz);

    /// 输入网格坐标数据
    void InputData(const vector<OCP_DBL>& coord, const vector<OCP_DBL>& zcorn);

    /// 输入COORD数据
    OCP_BOOL InputCOORDDATA(const vector<OCP_DBL>& coord);

    /// 输入ZCORN数据
    OCP_BOOL InputZCORNDATA(const vector<OCP_DBL>& zcorn);

    /// 设置角点
    void SetupCornerPoints();

    /// 设置所有标志
    void SetAllFlags(const HexahedronFace& oFace, const HexahedronFace& Face);

    /// 符号函数
    OCP_DBL OCP_SIGN(const OCP_DBL& x) { return x >= 0 ? 1 : -1; }

private:
    USI           nx, ny, nz;         ///< 网格在x, y, z方向上的数量
    OCP_DBL***    COORDDATA;          ///< COORD数据数组
    OCP_DBL****   ZCORNDATA;          ///< ZCORN数据数组
    Hexahedron*** cornerPoints;       ///< 角点数组
    OCP_ULL       numGrid;            ///< 网格数量
    OCP_ULL       numConn;            ///< 连接数量
    OCP_ULL       numConnMax;         ///< 最大连接数量
    vector<OCP_DBL> v;                ///< 体积向量
    vector<OCP_DBL> depth;            ///< 深度向量
    vector<OCP_DBL> dx, dy, dz;       ///< 网格尺寸向量
    vector<Point3D> center;           ///< 中心点向量
    vector<GeneralConnect> connect;   ///< 连接向量

    // 辅助变量
    OCP_BOOL       flagQuad;          ///< 四边形标志
    OCP_BOOL       upNNC, downNNC;    ///< 上下非邻近连接标志
    OCP_BOOL       flagJump;          ///< 跳跃标志
    HexahedronFace interFace;         ///< 界面

    /// 根据确定的坐标轴方向，如果沿y+方向放置，则flagForward等于1.0，否则等于-1.0，这与面法线矢量的计算有关
    OCP_DBL flagForward;
};

#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/16/2021      Create file                          */
/*  Chensong Zhang      Jan/16/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/