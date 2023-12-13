/*! \file    BulkConnTrans.hpp
 *  \brief   本文件包含BulkConnTrans类的声明以及相关方法。
 *  \author  Shizhe Li
 *  \date    Aug/24/2023
 *
 *  本文件定义了BulkConnTrans类及其派生类，用于计算体积连接传输。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BULKCONNAREA_HEADER__
#define __BULKCONNAREA_HEADER__

// OpenCAEPoroX header files
#include "BulkConnVarSet.hpp"
#include "Bulk.hpp"
using namespace std;

/*!
 *  \class   BulkConnTransMethod
 *  \brief   体积连接传输方法的基类。
 *
 *  该类定义了计算体积连接传输的基本框架。
 */
class BulkConnTransMethod{
public:
    /// 默认构造函数
    BulkConnTransMethod() = default;
    /// 纯虚函数，计算传输
    virtual void CalTrans(BulkConnPair& bp, const Bulk& bk) = 0;
};

/*!
 *  \class   BulkConnTransMethod01
 *  \brief   BulkConnTransMethod的派生类，实现特定的传输计算方法。
 *
 *  该类实现了CalTrans函数，提供了一种计算体积连接传输的方法。
 */
class BulkConnTransMethod01 : public BulkConnTransMethod{
public:
    /// 默认构造函数
    BulkConnTransMethod01() = default;
    /// 重写的函数，根据BulkConnPair和Bulk计算传输
    void CalTrans(BulkConnPair& bp, const Bulk& bk) override;
};

/*!
 *  \class   BulkConnTrans
 *  \brief   体积连接传输的主要类。
 *
 *  该类负责设置传输计算方法并调用相关函数进行计算。
 */
class BulkConnTrans{
public:
    /// 默认构造函数
    BulkConnTrans() = default;
    /// 设置传输计算方法
    void Setup();
    /*!
     *  \brief   计算体积连接的传输。
     *  \param   bp 体积连接对，包含连接的信息。
     *  \param   bk 体积对象，包含体积的信息。
     */
    void CalTrans(BulkConnPair& bp, const Bulk& bk) { bcaM->CalTrans(bp, bk); }

protected:
    /// 体积连接计算方法的指针，用于指向具体的计算方法。
    BulkConnTransMethod* bcaM;
};

#endif /* end if __BULKCONNTRANS_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/24/2023      Create file                          */
/*----------------------------------------------------------------------------*/
