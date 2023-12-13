/*! \file    GmshGrid.hpp
 *  \brief   GmshGrid类的声明，用于处理Gmsh生成的网格数据。
 *  \author  Shizhe Li
 *  \date    Sep/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifdef WITH_GMSH
#ifndef __GMSHGRID_HEADER__
#define __GMSHGRID_HEADER__

// 标准头文件
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <set>

// Gmsh头文件
#include <gmsh.h>

// OpenCAEPoroX头文件
#include "UtilInput.hpp"
#include "UtilMesh.hpp"

using namespace std;

/// \class Edge
/// \brief 表示一个边界的类，包含边界的起始节点、标签、物理信息等。
class Edge {
public:
    /// 构造函数，确保bId < eId。
    Edge(const OCP_ULL& b, const OCP_ULL& e) {
        if (b < e) {
            bId = b;
            eId = e;
        }
        else {
            bId = e;
            eId = b;
        }
    }
    /// 构造函数，包含面索引和局部边索引。
    Edge(const OCP_ULL& b, const OCP_ULL& e, const OCP_ULL& tag_in, const OCP_ULL& face, const OCP_ULL& localEdge) : Edge(b, e) {
        tag = tag_in;
        faceIndex.push_back(face);
        faceIndex.push_back(localEdge);
    }
    /// 构造函数，包含物理信息和物理索引。
    Edge(const OCP_ULL& b, const OCP_ULL& e, const OCP_ULL& tag_in, const string& physical_in, const USI& phyIndex_in) : Edge(b, e) {
        tag      = tag_in;
        physical = physical_in;
        phyIndex = pow(2, phyIndex_in);
    }
    /// 比较运算符，用于确保边的唯一性。
    auto operator <(const Edge& e) const {
        if (bId < e.bId || (bId == e.bId && eId < e.eId)) return OCP_TRUE;
        else                                              return OCP_FALSE;
    }
public:
    OCP_ULL                 bId; ///< 起始节点索引，bId < eId。
    OCP_ULL                 eId; ///< 结束节点索引。
    OCP_ULL                 tag; ///< 边的标签（用于调试）。
    string                  physical; ///< 物理信息（用于调试）。
    OCP_USI                 phyIndex{0}; ///< 物理信息的二进制编码。
    mutable vector<OCP_ULL> faceIndex; ///< 连接的面索引和边的局部索引。
    mutable vector<OCP_DBL> area; ///< 从bId开始的有效区域。
};

/// \class Polygon
/// \brief 表示一个多边形的类，包含多边形的顶点、标签、物理信息等。
class Polygon {
public:
    /// 构造函数。
    Polygon(const vector<OCP_ULL>& pIndex, const OCP_ULL& tag_in, const string& phyinfo, const OCP_USI& index);
    /// 计算中心点。
    void CalCenter(const vector<OCP_DBL>& points);
    /// 计算面积。
    void CalArea(const vector<OCP_DBL>& points);
    /// 判断一个点是否在元素内。
    OCP_BOOL IfPointInElement(const Point3D& objP, const vector<OCP_DBL>& points);

public:
    vector<OCP_ULL> p; ///< 顶点索引（按顺序或逆序存储）。
    OCP_ULL         tag; ///< 面的标签（用于调试）。
    string          physical; ///< 物理信息（用于调试）。
    OCP_USI         phyIndex; ///< 物理索引。
    string          boundary; ///< 边界（用于调试）。
    OCP_USI         boundIndex{0}; ///< 边界索引（二进制编码）。
    OCP_DBL         boundArea{0}; ///< 边界的有效面积。
    Point3D         center; ///< 中心点。
    OCP_DBL         area; ///< 面积。
};

/// \class Facies
/// \brief 表示岩相属性的类，包含岩相的名称、孔隙度、渗透率等属性。
class Facies {
public:
    /// 构造函数。
    Facies(const string& faciesname);

public:
    string  name; ///< 岩相名称。
    OCP_DBL poro{-1}; ///< 当前岩相的孔隙度。
    OCP_DBL kx{-1}; ///< 当前岩相的渗透率（x方向）。
    OCP_DBL ky{-1}; ///< 当前岩相的渗透率（y方向）。
    OCP_DBL kz{-1}; ///< 当前岩相的渗透率（z方向）。
};

/// \class GMSHGrid
/// \brief GMSHGrid类用于处理Gmsh生成的网格数据，目前仅支持二维网格。
class GMSHGrid {
public:
    /// 输入网格数据。
    void InputGrid(const string& file);
    /// 输入每个区域的属性。
    void InputProperty(ifstream& ifs);
    /// 获取边界名称。
    INT GetBoundaryName(vector<string>& names) const { 
        if (ifUse) {
            if (dimen == 2)      names = physicalNameSet[1];
            else if (dimen == 3) names = physicalNameSet[2];
            else                 OCP_ABORT("WRONG DIMENSION!");
            return names.size();
        }
        else {
            return 0;
        }
    }
protected:
    void InputGrid2D(const string& file);
    void Setup();
    void CalAreaCenter2D();
    void SetupConnAreaAndBoundary2D();

public:
    OCP_BOOL               ifUse{OCP_FALSE}; ///< 是否使用gmsh网格。
    USI                    dimen; ///< 维度。
    vector<OCP_DBL>        points; ///< 点的坐标。
    set<Edge>              edges; ///< 边（目前仅对2d）。
    vector<Polygon>        elements; ///< 元素（目前仅对2d）。
    vector<Facies>         facies; ///< 岩相。
    vector<vector<string>> physicalNameSet; ///< 物理名称集合。
    OCP_DBL                thickness; ///< 厚度（目前仅对2d）。
};

#endif /* end if __GMSHGRID_HEADER__ */
#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/
