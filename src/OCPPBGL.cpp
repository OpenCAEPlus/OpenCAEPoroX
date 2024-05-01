/*! \file    OCPPBGL.cpp
 *  \brief   OCPPBGL for OpenCAEPoroX simulator
 *  \author  Shizhe Li
 *  \date    May/01/2024
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPPBGL.hpp"


void GroupProcess(std::set<int>& cs_proc_group, MPI_Comm& cs_comm)
{

#ifdef WITH_PBGL

    {
        using namespace boost;
        using boost::graph::distributed::mpi_process_group;

        int world_rank, world_numproc;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_numproc);

        typedef adjacency_list<listS, distributedS<mpi_process_group, vecS>, undirectedS> Graph;
        Graph g(world_numproc);
        for (const auto& n : cs_proc_group) {
            add_edge(vertex(world_rank, g), vertex(n, g), g);
        }
        add_edge(vertex(world_rank, g), vertex(world_rank, g), g);
        synchronize(g);

        std::vector<int> local_components_vec(num_vertices(g));
        typedef iterator_property_map<std::vector<int>::iterator, property_map<Graph, vertex_index_t>::type> ComponentMap;
        ComponentMap component(local_components_vec.begin(), get(vertex_index, g));

        int num_components;
        // double start = MPI_Wtime();
        if (false) {         
            num_components = connected_components_ps(g, component);
        }
        else {
            num_components = connected_components(g, component);
        }
        //if (world_rank == 0) {
        //    std::cout << num_components << " components have been found! -- "
        //        << MPI_Wtime() - start << "s" << std::endl;
        //}


        if (cs_comm != MPI_COMM_NULL)   MPI_Comm_free(&cs_comm);
        MPI_Comm_split(MPI_COMM_WORLD, local_components_vec[0], world_rank, &cs_comm);

        int cs_numporc;
        MPI_Comm_size(cs_comm, &cs_numporc);
        std::vector<int> tmpR(cs_numporc);
        MPI_Allgather(&world_rank, 1, MPI_INT, tmpR.data(), 1, MPI_INT, cs_comm);

        cs_proc_group.clear();
        cs_proc_group.insert(tmpR.begin(), tmpR.end());
    }

#else
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // every rank is a group
    cs_proc_group.clear();
    cs_proc_group.insert(world_rank);
    if (cs_comm != MPI_COMM_NULL)   MPI_Comm_free(&cs_comm);
    MPI_Comm_split(MPI_COMM_WORLD, world_rank, world_rank, &cs_comm);

#endif
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           May/01/2024      Create file                          */
/*----------------------------------------------------------------------------*/