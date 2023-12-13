/*! \file    Output4Vtk.hpp
 *  \brief   用于输出VTK格式的油藏信息的类
 *  \author  Shizhe Li
 *  \date    Oct/19/2022
 *
 *  本文件包含Output4Vtk类的声明，该类用于创建VTK文件并输出油藏模拟的网格和属性数据。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OUTPUT4VTK_HEADER__
#define __OUTPUT4VTK_HEADER__

#include <fstream>
#include <string>
#include <vector>
using namespace std;

#include "OCPConst.hpp"

// 基本数据类型定义
typedef OCP_DBL VTK_DBL; ///< VTK文件中使用的双精度浮点数类型
typedef OCP_SIN VTK_SIN; ///< VTK文件中使用的有符号整数类型
typedef OCP_USI VTK_USI; ///< VTK文件中使用的无符号短整数类型
typedef OCP_ULL VTK_ULL; ///< VTK文件中使用的无符号长整数类型

// VTK文件基本关键字
const string  VTK_HEADER            = "# vtk DataFile Version 3.0"; ///< VTK文件头部信息
const string  VTK_ASCII             = "ASCII";                      ///< VTK文件的ASCII标记
const string  VTK_DATASET           = "DATASET";                    ///< VTK文件数据集关键字
const string  VTK_UNSTRUCTURED_GRID = "UNSTRUCTURED_GRID";          ///< VTK非结构化网格关键字
const string  VTK_POINTS            = "POINTS";                     ///< VTK点集关键字
const string  VTK_CELLS             = "CELLS";                      ///< VTK单元格集关键字
const string  VTK_CELL_TYPES        = "CELL_TYPES";                 ///< VTK单元格类型关键字
const string  VTK_CELL_DATA         = "CELL_DATA";                  ///< VTK单元格数据关键字
const string  VTK_POINT_DATA        = "POINT_DATA";                 ///< VTK点数据关键字
const VTK_USI VTK_MAX_TITLE_LENGTH  = 256;                          ///< VTK文件标题最大长度
const string  VTK_LOOKUP_TABLE      = "LOOKUP_TABLE";               ///< VTK查找表关键字
const string  VTK_DEFAULT           = "default";                    ///< VTK默认关键字
const string  VTK_SCALARS           = "SCALARS";                    ///< VTK标量数据关键字

// 基本单元格类型定义
const VTK_USI VTK_POLY_LINE  = 4;  ///< VTK多线段单元格类型
const VTK_USI VTK_TRIANGLE   = 5;  ///< VTK三角形单元格类型
const VTK_USI VTK_QUAD       = 9;  ///< VTK四边形单元格类型
const VTK_USI VTK_HEXAHEDRON = 12; ///< VTK六面体单元格类型

// 数据类型定义
const string VTK_FLOAT        = "float";         ///< VTK文件中使用的浮点数类型关键字
const string VTK_UNSIGNED_INT = "unsigned_int";  ///< VTK文件中使用的无符号整数类型关键字

/**
 * \class Output4Vtk
 * \brief 用于输出VTK格式文件的类
 *
 * 该类提供了创建VTK文件和输出网格信息及属性数据的方法。
 */
class Output4Vtk
{
    friend class Out4VTK; ///< 友元类声明

public:
    /// 创建一个新文件并写入通用信息
    OCP_ULL InitASCII(const string& dir, const string& myFile, const string& shortInfo) const;

    /**
     * \brief 输出CELL_DATA_SCALARS数据到VTK文件
     * \param[out] outVtk 输出流对象
     * \param[in] dataName 数据名称
     * \param[in] dataType 数据类型
     * \param[in] tmpV 数据向量
     * \param[in] bId 基础索引
     * \param[in] nb 数据数量
     * \param[in] digits 小数位数
     */
    template <typename T>
    void OutputCELL_DATA_SCALARS(ofstream&        outVtk,
                                 const string&    dataName,
                                 const string&    dataType,
                                 const vector<T>  tmpV,
                                 const OCP_ULL&   bId,
                                 const OCP_ULL&   nb,
                                 const USI&       digits) const;

public:
    /**
     * \brief 输出网格信息到VTK文件
     * \param[in] dir 目录路径
     * \param[in] nG 网格数
     * \param[in] points_xyz 点坐标集
     * \param[in] cell_points 单元格点索引集
     * \param[in] cell_type 单元格类型集
     */
    static void OutputGridInfo(const string& dir, const OCP_ULL& nG, const vector<OCP_DBL>& points_xyz,
                               const vector<OCP_ULL>& cell_points, const vector<USI>& cell_type);

protected:
    /**
     * \brief 从VTK文件读取网格信息
     * \param[in] dir 目录路径
     * \param[out] nG 网格数
     * \param[out] nP 点数
     * \param[out] points_xyz 点坐标集
     * \param[out] cell_points 单元格点索引集
     * \param[out] cell_type 单元格类型集
     */
    void InputGridInfo(const string& dir, OCP_ULL& nG, OCP_ULL& nP, vector<OCP_DBL>& points_xyz, vector<OCP_ULL>& cell_points, vector<USI>& cell_type) const;

protected:
    static const string tmpFile; ///< 临时文件路径
};

template <typename T>
void Output4Vtk::OutputCELL_DATA_SCALARS(ofstream&        outVtk,
                                         const string&    dataName,
                                         const string&    dataType,
                                         const vector<T>  tmpV,
                                         const OCP_ULL&   bId,
                                         const OCP_ULL&   nb,
                                         const USI&       digits) const
{
    outVtk << "\n" << VTK_SCALARS << " " << dataName << " " << dataType << " " << 1;
    outVtk << "\n" << VTK_LOOKUP_TABLE << " " << VTK_DEFAULT << "\n";
    outVtk << fixed << setprecision(digits);
    for (OCP_ULL n = 0; n < nb; n++) {
        outVtk << tmpV[bId + n] << "\n";
    }
}

#endif /* __OUTPUT4VTK_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/19/2022      Create file                          */
/*  Chensong Zhang      Feb/05/2023      Update output in vtk files           */
/*----------------------------------------------------------------------------*/
