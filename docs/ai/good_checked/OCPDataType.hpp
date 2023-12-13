/*! \file    OCPDataType.hpp 
 *  \brief   OpenCAEPoro中使用的数据类型定义
 *  \author  Shizhe Li
 *  \date    Oct/22/2023
 *
 *  本文件包含OpenCAEPoro项目中使用的基本数据类型定义。
 *  根据编译时的不同选项，这些类型可能会映射到不同的基础数据类型上。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPDATATYPE_HEADER__
#define __OCPDATATYPE_HEADER__

#include <string>

#ifdef OCPINT64
#define OCPINTTYPEWIDTH 64
#else
#define OCPINTTYPEWIDTH 32
#endif // OCPINT64

#if OCPINTTYPEWIDTH == 32
typedef int                OCP_SLL;  ///< 长整型有符号整数
typedef unsigned int       OCP_ULL;  ///< 长整型无符号整数
#define OCPMPI_SLL         MPI_INT
#define OCPMPI_ULL         MPI_UNSIGNED
#elif OCPINTTYPEWIDTH == 64
#include "petscksp.h"
// typedef long long          OCP_SLL;  ///< 长整型有符号整数
typedef PetscInt64        OCP_SLL;  ///< 64位有符号整数
// typedef unsigned long long OCP_ULL;  ///< 长整型无符号整数
typedef PetscInt64        OCP_ULL;  ///< 64位无符号整数
#define OCPMPI_SLL         MPI_LONG_LONG_INT
#define OCPMPI_ULL         MPI_UNSIGNED_LONG_LONG
#endif

#define OCPFLOATTYPEWIDTH 64

#if OCPFLOATTYPEWIDTH == 64
typedef double             OCP_DBL;  ///< 双精度浮点数
#define OCPMPI_DBL         MPI_DOUBLE
#elif  OCPFLOATTYPEWIDTH == 128
#ifdef BOOST_MATH_USE_FLOAT128
#include<boost/cstdfloat.hpp>
typedef _Quad              OCP_DBL;  ///< 128位浮点数
#define OCPMPI_DBL         MPI_LONG_DOUBLE
#endif
#endif // OCPFLOATTYPEWIDTH

// 内置数据类型
typedef unsigned int       USI;      ///< 通用无符号整数
typedef unsigned int       OCP_USI;  ///< OpenCAEPoro中的无符号整数
typedef int                INT;      ///< 通用有符号整数
typedef int                OCP_INT;  ///< OpenCAEPoro中的有符号整数
typedef float              OCP_SIN;  ///< 单精度浮点数
typedef unsigned int       OCP_BOOL; ///< OpenCAEPoro中的布尔类型
typedef char               OCP_CHAR; ///< 字符类型

#define OCPMPI_INT         MPI_INT
#define OCPMPI_ENUM        MPI_UNSIGNED
#define OCPMPI_BOOL        MPI_UNSIGNED
#define OCPMPI_CHAR        MPI_CHAR
#define OCPMPI_BYTE        MPI_BYTE

#endif // __OCPDATATYPE_HEADER__

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/22/2023      Create file                          */
/*----------------------------------------------------------------------------*/
