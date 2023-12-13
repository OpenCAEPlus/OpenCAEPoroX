/*! \file    UtilMesh.hpp
 *  \brief   UtilMesh 类的声明
 *  \author  Shizhe Li
 *  \date    Sep/10/2023
 *
 *  本文件包含了用于网格处理的各种几何对象的声明，如二维和三维点、六面体单元和面、以及三维矩阵。
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __UTILMESH_HEADER__
#define __UTILMESH_HEADER__

// 标准头文件
#include <math.h>
#include <stdlib.h>
#include <vector>

// OpenCAEPoroX 头文件
#include "OCPConst.hpp"

using namespace std;

/// \class Point2D
/// \brief 二维点类
class Point2D
{
public:
    OCP_DBL x; ///< x坐标值
    OCP_DBL y; ///< y坐标值

public:
    /// \brief 默认构造函数
    Point2D() = default;

    /// \brief 两参数构造函数
    /// \param x0 初始化x坐标
    /// \param y0 初始化y坐标
    Point2D(OCP_DBL x0, OCP_DBL y0)
        : x(x0)
        , y(y0) {};
};

/// \class Point3D
/// \brief 三维点类
class Point3D
{
public:
    OCP_DBL x; ///< x坐标值
    OCP_DBL y; ///< y坐标值
    OCP_DBL z; ///< z坐标值

public:
    /// \brief 默认构造函数
    Point3D() = default;

    /// \brief 三参数构造函数
    /// \param x0 初始化x坐标
    /// \param y0 初始化y坐标
    /// \param z0 初始化z坐标
    Point3D(OCP_DBL x0, OCP_DBL y0, OCP_DBL z0)
        : x(x0)
        , y(y0)
        , z(z0) {};

    /// \brief 从数组构造点
    /// \param points 包含三个坐标值的数组
    Point3D(const OCP_DBL* points) {
        x = points[0];
        y = points[1];
        z = points[2];
    }

    /// \brief 赋值运算符
    /// \param other 另一个Point3D对象
    /// \return 赋值后的当前对象引用
    Point3D& operator=(const Point3D& other);

    /// \brief 加法运算符
    /// \param other 另一个Point3D对象
    /// \return 相加后的新Point3D对象
    Point3D  operator+(const Point3D& other) const;

    /// \brief 减法运算符
    /// \param other 另一个Point3D对象
    /// \return 相减后的新Point3D对象
    Point3D  operator-(const Point3D& other) const;

    /// \brief 点积运算符
    /// \param other 另一个Point3D对象
    /// \return 点积结果
    OCP_DBL  operator*(const Point3D& other) const;

    /// \brief 加法赋值运算符
    /// \param other 另一个Point3D对象
    /// \return 加法赋值后的当前对象引用
    Point3D& operator+=(const Point3D& other);

    /// \brief 数乘赋值运算符
    /// \param a 乘数
    /// \return 数乘赋值后的当前对象引用
    Point3D& operator*=(const OCP_DBL& a);

    /// \brief 除法赋值运算符
    /// \param a 除数
    /// \return 除法赋值后的当前对象引用
    Point3D& operator/=(const OCP_DBL& a);

    /// \brief 重置点坐标为原点
    void     Reset() {
        x = 0;
        y = 0;
        z = 0;
    };
};

/// \brief 点与数乘
/// \param p Point3D对象
/// \param a 乘数
/// \return 乘法结果新Point3D对象
Point3D operator*(const Point3D& p, const OCP_DBL& a);

/// \brief 数与点乘
/// \param a 乘数
/// \param p Point3D对象
/// \return 乘法结果新Point3D对象
Point3D operator*(const OCP_DBL& a, const Point3D& p);

/// \brief 计算两个Point3D的叉积
/// \param p1 第一个Point3D对象
/// \param p2 第二个Point3D对象
/// \return 叉积结果新Point3D对象
Point3D CrossProduct(const Point3D& p1, const Point3D& p2);

/// A hexahedron cell.
//      p3-------p2
//      / |      /|
//      p0------p1|
//      | |     | |
//      | p7----|-p6
//      |/      |/
//     p4------p5   
//

/// \class Hexahedron
/// \brief 六面体单元类
class Hexahedron
{
public:
    /// \brief 计算六面体体积，仅适用于凸多面体
    /// \return 六面体体积
    OCP_DBL CalVolum() const;

    /// \brief 计算六面体中心
    /// \return 六面体中心点
    Point3D CalCenter() const;

public:
    Point3D p0, p1, p2, p3, p4, p5, p6, p7; ///< 六面体的八个顶点
};

/// \class HexahedronFace
/// \brief 六面体面类
class HexahedronFace
{
public:
    /// \brief 计算面的面积法向量
    /// \return 面的面积法向量
    Point3D CalAreaVector() const;

    /// \brief 计算面的中心
    /// \return 面的中心点
    Point3D CalCenter() const;

public:
    Point3D p0, p1, p2, p3; ///< 面的四个顶点
};

/// \class Matrix3
/// \brief 3x3矩阵类
class Matrix3
{
public:
    OCP_DBL M[3][3]; ///< 矩阵元素

    /// \brief 与Point3D相乘
    /// \param v Point3D对象
    /// \return 相乘结果新Point3D对象
    Point3D operator*(const Point3D& v) const;
};

/// \brief 计算两条线的交点
/// \param Line1 第一条线的两个点
/// \param Line2 第二条线的两个点
/// \return 交点Point2D对象
Point2D CalCrossingPoint(const Point2D Line1[2], const Point2D Line2[2]);

/// \brief 计算两个HexahedronFace的交界面积（交界面非四边形）
/// \param FACE1 第一个HexahedronFace对象
/// \param FACE2 第二个HexahedronFace对象
/// \return 交界面积
OCP_DBL CalAreaNotQuadr(const HexahedronFace& FACE1, const HexahedronFace& FACE2);

#endif /* end if __UTILMESH_HEADER__ */


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/
