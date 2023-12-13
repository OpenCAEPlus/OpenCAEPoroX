/*! \file    Domain.hpp
 *  \brief   域类声明
 *  \author  Shizhe Li
 *  \date    Feb/28/2023
 *
 *  \note    在OpenCAEPoroX中使用的参数大多与SLB的Eclipse兼容，
 *           但它有一些自己的规则以便于使用。它是可扩展的和用户友好的。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team.
 *  All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __DOMAIN_HEADER__
#define __DOMAIN_HEADER__

#include "OCPConst.hpp"
#include "Partition.hpp"
#include "UtilTiming.hpp"
#include "OCPTimeRecord.hpp"
#include <parmetis.h>
#include <vector>
#include <map>

using namespace std;

/*!
 * \class Domain
 * \brief Domain类负责在分区后的新分布。
 *        所有元素（井）都是活跃的，所有连接（射孔）都是活跃的，
 *        死亡元素被排除，因此全局索引也是针对活跃元素的。
 */
class Domain
{
	friend class Reservoir;
	friend class Bulk;
	friend class Well;
	friend class PeacemanWell;
	friend class LinearSystem;
	friend class OCPMatrix;
	// Method(tmp)
	friend class IsoT_IMPEC;
	friend class IsoT_FIM;
	friend class IsoT_AIMc;
	friend class T_FIM;
	// Linear Solver (tmp)
	friend class SamgSolver;

public:
	void Setup(const Partition& part, const PreParamGridWell& gridwell);
	auto GetNumGridTotal() const { return numElementTotal - numWellTotal; }
	auto GetNumGridInterior() const { return numGridInterior; }
	auto GetWell() const { return well; }
	const auto& GetGrid() const { return grid; }

public:
	MPI_Comm                     myComm; ///< MPI通信器
	OCP_INT                      numproc; ///< 进程总数
	OCP_INT                      myrank; ///< 当前进程的排名

protected:
	/// 元素总数（网格 + 井）
	OCP_ULL numElementTotal;
	/// 井的总数
	OCP_USI numWellTotal;
	/// 本地元素数：numGridInterior + numWellLocal
	OCP_USI numElementLocal;
	/// 当前进程存储的内部网格数
	OCP_USI numGridInterior;
	/// 当前进程存储的幽灵网格数
	OCP_USI numGridGhost;
	/// numGridInterior + numGridGhost
	OCP_USI numGridLocal;
	/// 当前进程存储的井数
	USI     numWellLocal;
	/// 从其他进程接收的元素CSR表示
	vector<vector<idx_t>>   elementCSR;
	///< 内部网格和幽灵网格的全局索引
	vector<OCP_ULL>         grid;
	/// 井的全局索引（索引从零开始）
	vector<OCP_USI>         well;
	/// 井的全局索引（从零开始），射孔的索引（从零开始）以及射孔的位置（块的本地索引）
	vector<vector<OCP_USI>> wellWPB;
	/// 网格的所有邻居数：包括井，包括自身
	vector<USI>             neighborNum;
	/// 全局索引 -> 本地索引（内部网格和幽灵网格）
	map<OCP_ULL, OCP_USI>   init_global_to_local;

	////////////////////////////////////////
	// 通信
	////////////////////////////////////////
public:
	const vector<OCP_ULL>* CalGlobalIndex(const USI& nw) const;

	////////////////////////////////////////
	// 隐式通信（首选，本地索引）
	////////////////////////////////////////
	USI numSendProc; ///< 发送进程数
	USI numRecvProc; ///< 接收进程数
	vector<vector<OCP_USI>> send_element_loc; ///< 发送元素的本地索引
	// vector<vector<OCP_CHAR>> send_buffer;
	vector<vector<OCP_USI>> recv_element_loc; ///< 接收元素的本地索引
	// vector<vector<OCP_CHAR>> recv_buffer;
	mutable vector<MPI_Request>  send_request; ///< MPI发送请求
	mutable vector<MPI_Request>  recv_request; ///< MPI接收请求

	////////////////////////////////////////
	// 全局索引通信
	////////////////////////////////////////
	mutable vector<OCP_ULL>     global_index;  ///< 方程中的内部网格 + 活跃井 + 幽灵网格

	// 井射孔
public:
	OCP_INT GetPerfLocation(const OCP_USI& wId, const USI& p) const; ///< 获取射孔位置
	/// 返回特定井（给定井全局索引）的射孔数量
	USI     GetPerfNum(const OCP_USI& wId) const;
};

#endif /* __DOMAIN_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Feb/28/2023      Create file                          */
/*----------------------------------------------------------------------------*/
