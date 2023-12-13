/*!
 * \file    OCPBoundary.hpp
 * \brief   OCPBoundary类的声明
 * \author  Shizhe Li
 * \date    2023年9月22日
 *
 * 这个文件包含了OCPBoundary类的声明。OCPBoundary类主要用于处理石油模拟中的边界问题。
 */

#ifndef __OCPBOUNDARY_HEADER__
#define __OCPBOUNDARY_HEADER__

#include "ParamReservoir.hpp"
#include "BulkVarSet.hpp"
#include "HeatLoss.hpp"
#include "BoundaryFlow.hpp"
#include <vector>

using namespace std;

/*!
 * \class   OCPBoundary
 * \brief   用于处理石油模拟中的边界问题的类
 *
 * OCPBoundary类包含了处理石油模拟中的边界问题所需要的所有成员变量和成员函数。
 */
class OCPBoundary{
	friend class BulkAccumuTerm01;

public:
	/*!
	 * \brief   设置OCPBoundary类的参数
	 * \param   rs_param  ParamReservoir类的实例，包含了储层的参数
	 * \param   bvs       BulkVarSet类的实例，包含了bulk变量的集合
	 */
	void Setup(const ParamReservoir& rs_param, const BulkVarSet& bvs) {
		heatLoss.Setup(rs_param, bvs, boundIndex);
		boundaryFlow.Setup(rs_param, bvs, boundName, boundIndex);
	}

	/*!
	 * \brief   重置到上一个时间步
	 */
	void ResetToLastTimeStep() { heatLoss.ResetToLastTimeStep(); }

	/*!
	 * \brief   更新到上一个时间步
	 */
	void UpdateLastTimeStep() { heatLoss.UpdateLastTimeStep(); }

public:
	/// Heat loss term
	HeatLoss     heatLoss;    ///< 热损失项
	/// Boundary flow term
	BoundaryFlow boundaryFlow;///< 边界流动项

public:
	auto& GetBoundaryIndex() { return boundIndex; } ///< 获取边界索引
	auto& GetBoundaryArea() { return boundArea; }   ///< 获取边界面积
	auto& GetBoundName() { return boundName; }      ///< 获取边界名称

protected:
	vector<string>  boundName;   ///< 边界名称，类型为string的vector
	vector<USI>     boundIndex;  ///< 边界索引，类型为USI的vector
	vector<OCP_DBL> boundArea;   ///< 边界面积，类型为OCP_DBL的vector
};

#endif /* end if __OCPBoundary_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/22/2023      Create file                          */
/*----------------------------------------------------------------------------*/
