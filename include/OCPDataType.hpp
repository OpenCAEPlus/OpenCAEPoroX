/*! \file    OCPDataType.hpp
 *  \brief   Data type used in OpenCAEPoro
 *  \author  Shizhe Li
 *  \date    Oct/22/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPDATATYPE_HEADER__
#define __OCPDATATYPE_HEADER__


#include <string>


#define OCPINTTYPEWIDTH 32


#if OCPINTTYPEWIDTH == 32
typedef int                OCP_SLL;  ///< Long long signed integer
typedef unsigned int       OCP_ULL;  ///< Long long unsigned integer

#define OCPMPI_SLL         MPI_INT
#define OCPMPI_ULL         MPI_UNSIGNED

#elif OCPINTTYPEWIDTH == 64
typedef long long          OCP_SLL;  ///< Long long signed integer
//typedef PetscInt64         OCP_SLL;  ///< Long long signed integer
typedef unsigned long long OCP_ULL;  ///< Long long unsigned integer

#define OCPMPI_SLL         MPI_LONG_LONG_INT
#define OCPMPI_ULL         MPI_UNSIGNED_LONG_LONG

#endif


#define OCPFLOATTYPEWIDTH 64


#if OCPFLOATTYPEWIDTH == 64

typedef double             OCP_DBL;  ///< float 64 precision
#define OCPMPI_DBL         MPI_DOUBLE

#elif  OCPFLOATTYPEWIDTH == 128
#ifdef BOOST_MATH_USE_FLOAT128

#include<boost/cstdfloat.hpp>

typedef _Quad              OCP_DBL;  ///< float 128 precision
#define OCPMPI_DBL         MPI_LONG_DOUBLE

#endif

#endif


 // Build-in data type
typedef unsigned int       USI;      ///< Generic unsigned integer
typedef unsigned int       OCP_USI;  ///< unsigned integer
typedef int                INT;      ///< Generic signed integer
typedef int                OCP_INT;  ///< integer
typedef float              OCP_SIN;  ///< Single precision
typedef unsigned int       OCP_BOOL; ///< OCP_BOOL in OCP
typedef char               OCP_CHAR; ///< Char


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