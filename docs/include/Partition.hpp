/*! \file    Partition.hpp
 *  \brief   本文件包含Partition类的声明。
 *  \author  Shizhe Li
 *  \date    Feb/21/2023
 *
 *  \note    OpenCAEPoroX使用的参数与SLB的Eclipse大体兼容，
 *           但它有自己的一些规则以便于使用。它是可扩展的和用户友好的。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARTITION_HEADER__
#define __PARTITION_HEADER__

// OpenCAEPoroX头文件
#include "PreParamGridWell.hpp"
#include "UtilTiming.hpp"
#include "OCPTimeRecord.hpp"

// ParMetis
#define rabs fabsf     // 与fasp冲突
#include <parmetis.h>
#undef  rabs

// 标准库
#include <algorithm>
#include <set>
using namespace std;

/*!
 * \class Partition
 * \brief 使用parmetis进行图分割，分割的数量与处理器的数量相同，在分割之后将变得无用。
 */
class Partition
{
	friend class Domain;
	friend class Reservoir;

public:
	/*!
	 * \brief 初始化MPI环境。
	 * \param comm MPI通信器。
	 */
	void InitMPI(MPI_Comm comm);

	/*!
	 * \brief 设置分区参数。
	 * \param grid 包含网格和井的预处理参数。
	 */
	void SetPartition(const PreParamGridWell& grid);

	/*!
	 * \brief 设置分布。
	 */
	void SetDistribution();

protected:
	/*!
	 * \brief 初始化分区参数。
	 */
	void InitParam();

	MPI_Comm    myComm;           ///< MPI通信器。
	OCP_INT     numproc, myrank;  ///< 进程数和当前进程的排名。
	idx_t       numWellTotal;     ///< 全部井的数量。
	idx_t       numElementTotal;  ///< 全部网格的数量（包括井）。
	idx_t       numElementLocal;  ///< 本地网格的数量（包括井）。
	idx_t       numEdgesLocal;    ///< 本地边的数量。
	idx_t*      numEdges;         ///< 每个进程边的数量。
	idx_t		maxNumVtx;        ///< 所有进程中顶点的最大数量。
	idx_t		maxNumEdges;      ///< 所有进程中边的最大数量。

	// parmetis的图数据结构
	idx_t*     vtxdist;          ///< 顶点分布。
	idx_t*     xadj;             ///< 相邻顶点数组。
	idx_t*     adjncy;           ///< 邻接顶点数组。
	idx_t*     vwgt;             ///< 顶点权重数组。
	idx_t*     adjwgt;           ///< 边权重数组。
	idx_t      wgtflag;          ///< 权重标志。
	idx_t      numflag;          ///< 数字标志。
	idx_t      ncon;             ///< 平衡约束的数量。
	idx_t      nparts;           ///< 分区的数量。
	real_t*    tpwgts;           ///< 目标分区权重。
	real_t     ubvec;            ///< 不平衡容忍度。
	idx_t*     options;          ///< 选项数组。
	idx_t      edgecut;          ///< 边切割数量。
	idx_t*     part;             ///< 分区数组。
	mutable vector<vector<idx_t>> elementCSR;  ///< 元素的CSR表示。

	// 辅助变量
	struct NeighborP
	{
		idx_t index, location;  ///< 邻居索引和位置。
		/*!
		 * \brief 构造函数。
		 * \param id 邻居的索引。
		 * \param lc 邻居的位置。
		 */
		NeighborP(const idx_t& id, const idx_t& lc) :index(id), location(lc) {};

		/*!
		 * \brief 比较函数。
		 * \param n1 第一个邻居。
		 * \param n2 第二个邻居。
		 * \return 比较结果。
		 */
		static bool less(const NeighborP& n1, const NeighborP& n2) { return n1.index < n2.index; }
	};
};

#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Feb/21/2023      Create file                          */
/*----------------------------------------------------------------------------*/
