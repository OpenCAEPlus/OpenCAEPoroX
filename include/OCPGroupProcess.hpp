/*! \file    OCPGroupProcess.hpp
 *  \brief   OCPGroupProcess class declaration
 *  \author  Shizhe Li
 *  \date    May/02/2024
 *
 *  \note    The params used in OpenCAEPoroX is mostly compatible with Eclipse by SLB,
 *           but it has some own rules for easy to use. It is extensible and friendly.
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


#ifndef __OCPGROUPPROCESS_HEADER__
#define __OCPGROUPPROCESS_HEADER__


#include "OCPDataType.hpp"

#include <vector>
#include <iostream>
#include <set>
#include <unordered_map>
#include <mpi.h>

#ifdef WITH_PBGL
#include <boost/filesystem.hpp>
#include <boost/graph/use_mpi.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/distributed/connected_components_parallel_search.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/graph/parallel/distribution.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#endif

#ifdef WITH_BGL
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#endif


#ifdef WITH_PARMETIS
#define rabs fabsf     // conflict with fasp
#include <parmetis.h>
#undef  rabs
#endif


#ifdef WITH_METIS
#define rabs fabsf     // conflict with fasp
#include <metis.h>
#undef  rabs
#endif


enum class GroupMethod
{
	PBGL,
	PBGL_gather,
	BGL,
	ParMetis,
	ParMetis_gather,
	Metis,
	Self
};


// each process must know all its neighbors
void GroupProcess(const GroupMethod& method, std::set<OCP_INT>& cs_proc_group,
	const std::unordered_map<OCP_INT, OCP_INT>& proc_weight,
	MPI_Comm& cs_comm, const MPI_Comm& global_comm);


#endif



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           May/02/2024      Create file                          */
/*----------------------------------------------------------------------------*/