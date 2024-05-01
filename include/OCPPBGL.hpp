/*! \file    OCPPBGL.hpp
 *  \brief   OCPPBGL class declaration
 *  \author  Shizhe Li
 *  \date    May/01/2024
 *
 *  \note    The params used in OpenCAEPoroX is mostly compatible with Eclipse by SLB,
 *           but it has some own rules for easy to use. It is extensible and friendly.
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


#ifndef __OCPPBGL_HEADER__
#define __OCPPBGL_HEADER__

#include <vector>
#include <iostream>
#include <set>
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

// use parallel boost graph library to group processes
void GroupProcess(std::set<int>& cs_proc_group, MPI_Comm& cs_comm);

extern double TIME_PBGL;

#endif



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           May/01/2024      Create file                          */
/*----------------------------------------------------------------------------*/