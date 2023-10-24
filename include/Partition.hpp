/*! \file    Partition.hpp
 *  \brief   Partition class declaration
 *  \author  Shizhe Li
 *  \date    Feb/21/2023
 *
 *  \note    The params used in OpenCAEPoroX is mostly compatible with Eclipse by SLB,
 *           but it has some own rules for easy to use. It is extensible and friendly.
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARTITION_HEADER__
#define __PARTITION_HEADER__


// OpenCAEPoroX header files
#include "PreParamGridWell.hpp"
#include "UtilTiming.hpp"
#include "OCPTimeRecord.hpp"

// ParMetis
#define rabs fabsf     // conflict with fasp
#include <parmetis.h>
#undef  rabs

// standard
#include <algorithm>
#include <set>

using namespace std;



/// Graph Partition with parmetis, num of partiotion are the same as proc's num,
//  it's useless after partition.
class Partition
{
	friend class Domain;
	friend class Reservoir;

public:

	void InitMPI(MPI_Comm comm);
	void SetPartition(const PreParamGridWell& grid);
	void SetDistribution();
	
protected:
	void InitParam();

	MPI_Comm    myComm;
	OCP_INT     numproc, myrank;


	idx_t       numWellTotal;     ///< num of total wells
	idx_t       numElementTotal;  ///< num of total grid(wells are included)	
	idx_t       numElementLocal;  ///< num of local grid(wells are included)
	idx_t       numEdgesLocal;    ///< num of local edge
	idx_t*      numEdges;         ///< num of each process' edge
	idx_t		maxNumVtx;        ///< max num of vertex for all processes
	idx_t		maxNumEdges;      ///< max num of edges  for all processes

	

	// graph for parmetis
	idx_t*     vtxdist;
	idx_t*     xadj;
	idx_t*     adjncy;
	idx_t*     vwgt;
	idx_t*     adjwgt;
	idx_t      wgtflag;
	idx_t      numflag;
	idx_t      ncon;
	idx_t      nparts;
	real_t*    tpwgts;
	real_t     ubvec;
	idx_t*     options;
	idx_t      edgecut;
	idx_t*     part;


	mutable vector<vector<idx_t>> elementCSR;

	// Auxiliary vars
	struct NeighborP
	{
		idx_t index, location;
		NeighborP(const idx_t& id, const idx_t& lc) :index(id), location(lc) {};
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