/*! \file    OCPUnits.hpp 
 *  \brief   OpenCAEPoro中使用的单位
 *  \author  Shizhe Li 
 *  \date    Oct/22/2023 
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OPENCAEPORO_UNITS_HEADER__
#define __OPENCAEPORO_UNITS_HEADER__

// OpenCAEPoroX 头文件
#include "OCPDataType.hpp"
#include "UtilError.hpp"

using namespace std;

enum class UnitType : USI{
	/// 场地单位
	FIELD,
	/// 公制单位
	METRIC
};

void SetUnit(const string& ut);

/// OCP单位类型
extern UnitType unitType;

//////////////////////////////////////////////////
/// 场地单位
//////////////////////////////////////////////////

/// 1 [bbl] = 5.61458 [ft3]
const OCP_DBL FIELD_CONV1           = 5.61458;

/// 1 [Mscf] = 1000 [ft3]
const OCP_DBL FIELD_CONV2           = 1000;

/// 0 [F] = 459.67 [R]
const OCP_DBL FIELD_CONV3           = 459.67;

/// 1 [btu] = 5.40395 [psia]��[ft3]
const OCP_DBL FIELD_CONV4           = 5.40395;

/// 1 [lbm]/[ft3] = 0.0160185 [gm-M]/[cc]
const OCP_DBL FIELD_CONV5           = 0.0160185;

/// 气体常数 ([ft3]��[psi]/[R]/[lb-M])
const OCP_DBL FIELD_GAS_CONSTANT    = 10.73159;

/// 重力常数 ([ft2]��[psi]/[lb])
const OCP_DBL FIELD_GRAVITY_FACTOR  = 0.00694444;

/// 达西常数 (([ft3/day]��[cp])/([md]��[ft]��[psi]))
const OCP_DBL FIELD_DARCY           = 0.001127 * FIELD_CONV1;

/// 标准条件下的水密度 ([lb]/[ft3])
const OCP_DBL FIELD_RHOW_STD        = 62.3664;

/// 标准条件下的空气密度 ([lb]/[ft3])
const OCP_DBL FIELD_RHOAIR_STD      = 0.076362;

/// 大气压力 ([psia])
const OCP_DBL FIELD_PRESSURE_STD    = 14.7;

/// 标准温度 ([F])
const OCP_DBL FIELD_TEMPERATURE_STD = 60;

/// 时间单位
const string  FIELD_TIME            = " Day ";

//////////////////////////////////////////////////
/// 公制单位
//////////////////////////////////////////////////

/// 1 [m3] = 1 [m3]
const OCP_DBL METRIC_CONV1           = 1.0;

/// 0 [C] = 273.15 [K]
const OCP_DBL METRIC_CONV2           = 273.15;

/// 1 [kj] = 0.01 [bar]��[m3]
const OCP_DBL METRIC_CONV3           = 0.01;

/// 1 [kg-M]/[m3] = 0.001 [gm-M]/[cc]
const OCP_DBL METRIC_CONV4           = 0.001;

/// 气体常数 ([m3]��[bars]/[K]/[kg-M])
const OCP_DBL METRIC_GAS_CONSTANT    = 0.083143;

/// 重力常数 ([m2]��[bars]/[kg])
const OCP_DBL METRIC_GRAVITY_FACTOR  = 0.0000980665;

/// 达西常数  (([m3/day]��[cP])/([md]��[m]��[bars]))
const OCP_DBL METRIC_DARCY           = 1.0;

/// 标准条件下的水密度 ([kg]/[m3])
const OCP_DBL METRIC_RHOW_STD        = 999.014;

/// 标准条件下的空气密度 ([kg]/[m3])
const OCP_DBL METRIC_RHOAIR_STD      = 1.2232;

/// 大气压力 ([barsa])
const OCP_DBL METRIC_PRESSURE_STD    = 1.01325;

/// 标准温度 ([C])
const OCP_DBL METRIC_TEMPERATURE_STD = 20;

/// 时间单位
const string METRIC_TIME = " Sec ";

//////////////////////////////////////////////////
/// 工作单位
//////////////////////////////////////////////////

/// 达西常数
extern OCP_DBL CONV_DARCY;

/// 气体常数
extern OCP_DBL GAS_CONSTANT;

/// 重力常数
extern OCP_DBL GRAVITY_FACTOR;

/// 标准条件下的水密度
extern OCP_DBL RHOW_STD;

/// 标准条件下的空气密度
extern OCP_DBL RHOAIR_STD;

/// 大气压力
extern OCP_DBL PRESSURE_STD;

/// 标准温度
extern OCP_DBL TEMPERATURE_STD;

/// [bbl] -> [ft3], [m3] -> [m3]
extern OCP_DBL CONV1;

/// [Mscf] -> [ft3], [m3] -> [m3]
extern OCP_DBL CONV2;

/// [Mscf] -> [bbl], [m3] -> [m3]
extern OCP_DBL CONV3;

/// [F] -> [R], [C] -> [K]
extern OCP_DBL CONV4;

/// [btu] -> [psia]��[ft3], [kj] -> [bar]��[m3]
extern OCP_DBL CONV5;

/// [lbm]/[ft3] -> [gm-M]/[cc], [kg-M]/[m3] -> [gm-M]/[cc]
extern OCP_DBL CONV6;

/// 时间单位
extern string TIMEUNIT;

#endif // __OPENCAEPORO_UNITS_HEADER__

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/22/2023      Create file                          */
/*----------------------------------------------------------------------------*/

