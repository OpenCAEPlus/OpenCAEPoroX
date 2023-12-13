/*! \file    OCPFuncTable.hpp 
 *  \brief   Table Functions in OCP 
 *  \author  Shizhe Li 
 *  \date    Jul/11/2023 
 * 
 *----------------------------------------------------------------------------------- 
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved. 
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later. 
 *----------------------------------------------------------------------------------- 
 */

#ifndef __OCPFUNCTABLE_HEADER__
#define __OCPFUNCTABLE_HEADER__

// OpenCAEPoroX header files
#include "OCPTable.hpp"

using namespace std;

///////////////////////////////////////////////////////
// OCPFuncTable
/////////////////////////////////////////////////////
/**
 * \class OCPFuncTable
 * \brief 提供OCP表功能
 */
class OCPFuncTable
{
public:
    OCPFuncTable() = default;

    /**
     * \brief 设置函数
     * \param src 数据源
     */
    void Setup(const vector<vector<OCP_DBL>>& src) {
        table.Setup(src);
        data.resize(table.GetColNum());
        cdata.resize(table.GetColNum());
    }

    /**
     * \brief 检查表是否为空
     * \return 表是否为空
     */
    OCP_BOOL IsEmpty() const { return table.IsEmpty(); }

protected:
    OCPTable                  table;   ///< OCP表
    mutable vector<OCP_DBL>   data;    ///< 数据
    mutable vector<OCP_DBL>   cdata;   ///< 缓存数据
};

///////////////////////////////////////////////////////
// OCPFuncTable2
/////////////////////////////////////////////////////
/**
 * \class OCPFuncTable2
 * \brief 提供OCP表功能(版本2)
 */
class OCPFuncTable2
{
public:
    OCPFuncTable2() = default;

    /**
     * \brief 设置函数
     * \param tab 表格数据
     */
    void Setup(const Table2& tab) {
        table.Setup(tab);
        data.resize(table.GetColNum());
        cdata1.resize(table.GetColNum());
        cdata2.resize(table.GetColNum());
    }

    /**
     * \brief 检查表是否为空
     * \return 表是否为空
     */
    OCP_BOOL IsEmpty() const { return table.IsEmpty(); }

protected:
    OCPTable2                 table;   ///< OCP表(版本2)
    mutable vector<OCP_DBL>   data;    ///< 数据
    mutable vector<OCP_DBL>   cdata1;  ///< 缓存数据1
    mutable vector<OCP_DBL>   cdata2;  ///< 缓存数据2
};

#endif // __OCPFUNCTABLE_HEADER__

/*----------------------------------------------------------------------------
 *  Brief Change History of This File                                         
 *----------------------------------------------------------------------------
 *  Author              Date             Actions                              
 *----------------------------------------------------------------------------
 *  Shizhe Li           Jul/11/2023      Create file                          
 *----------------------------------------------------------------------------
 */