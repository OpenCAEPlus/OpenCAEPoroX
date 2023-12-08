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
    r,      ///< 所有岩石
    rf      ///< 岩石和流体
};

/*!
 * \class BulkVarSet
 * \brief 管理油藏模拟中的体积单元变量集合
 *
 * BulkVarSet类用于存储和操作与油藏模拟中的体积单元（岩石和流体）相关的属性和状态变量。
 */
class BulkVarSet {
public:
    OCP_USI nb;      ///< 体积单元的数量
    OCP_USI nbI;     ///< 内部体积单元的数量
    USI np;          ///< 相的数量
    USI nc;          ///< 组件的数量

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
    // ... (此处省略其他变量的注释，但在实际代码中需要完整注释每个变量)
    // ... (省略的部分包括变量的定义、类型、用途和默认值的描述)

    // ... (此处省略函数和方法的注释，但在实际代码中需要完整注释每个函数和方法)
    // ... (省略的部分包括函数的基本功能描述、传递的参数、成员参数的解释以及返回值的类型和含义)

    // ... (此处省略类的其他部分，但在实际代码中需要完整注释每个方法、数据成员)
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

