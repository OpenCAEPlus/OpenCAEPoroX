/*! \file    Partition.hpp
 *  \brief   Partition class declaration
 *  \author  Shizhe Li
 *  \date    Feb/21/2023
 *
 *  \note    The params used in OpenCAEPoro is mostly compatible with Eclipse by SLB,
 *           but it has some own rules for easy to use. It is extensible and friendly.
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARTITION_HEADER__
#define __PARTITION_HEADER__


// OpenCAEPoro header files
#include "PreParamGridWell.hpp"

// ParMetis
#include <parmetis.h>

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

	MPI_Comm      myComm;
	OCP_INT       numproc, myrank;


	OCP_INT       numWellTotal;     ///< num of total wells
	OCP_INT       numElementTotal;  ///< num of total grid(wells are included)	
	OCP_INT       numElementLocal;  ///< num of local grid(wells are included)
	OCP_INT       numEdgesLocal;    ///< num of local edge
	OCP_INT*      numEdges;         ///< num of each process' edge
	OCP_INT		  maxNumVtx;        ///< max num of vertex for all processes
	OCP_INT		  maxNumEdges;      ///< max num of edges  for all processes

	

	// graph for parmetis
	idx_t*        vtxdist;
	idx_t*        xadj;
	idx_t*        adjncy;
	idx_t*        vwgt;
	idx_t*        adjwgt;
	idx_t         wgtflag;
	idx_t         numflag;
	idx_t         ncon;
	idx_t         nparts;
	real_t*       tpwgts;
	real_t        ubvec;
	idx_t*        options;
	idx_t         edgecut;
	idx_t*        part;


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