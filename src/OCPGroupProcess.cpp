/*! \file    OCPGroupProcess.cpp
 *  \brief   OCPGroupProcess for OpenCAEPoroX simulator
 *  \author  Shizhe Li
 *  \date    May/02/2024
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPGroupProcess.hpp"


#ifdef WITH_PBGL

 // each process with one vertex
static void GroupProcessPBGL01(std::set<int>& cs_proc_group, MPI_Comm& cs_comm, const MPI_Comm& global_comm)
{
    {
        double start0 = MPI_Wtime();

        using namespace boost;
        using boost::graph::distributed::mpi_process_group;

        int world_rank, world_numproc;
        MPI_Comm_rank(global_comm, &world_rank);
        MPI_Comm_size(global_comm, &world_numproc);

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
        double start = MPI_Wtime();
        if (false) {
            num_components = connected_components_ps(g, component);
        }
        else {
            num_components = connected_components(g, component);
        }
        if (world_rank == 0) {
            std::cout << num_components << " components have been found! -- "
                << MPI_Wtime() - start << "s" << std::endl;
        }

        if (cs_comm != MPI_COMM_NULL)   MPI_Comm_free(&cs_comm);
        MPI_Comm_split(global_comm, local_components_vec[0], world_rank, &cs_comm);

        int cs_numporc;
        MPI_Comm_size(cs_comm, &cs_numporc);
        std::vector<int> tmpR(cs_numporc);
        MPI_Allgather(&world_rank, 1, MPI_INT, tmpR.data(), 1, MPI_INT, cs_comm);

        cs_proc_group.clear();
        cs_proc_group.insert(tmpR.begin(), tmpR.end());
    }
}


static void GroupProcessPBGL02(std::set<int>& cs_proc_group, MPI_Comm& cs_comm, const MPI_Comm& global_comm)
{
    /// number of vertex in each process, remember that numproc = numvertex
    const int N = 1000;

    {
        int world_rank, world_numproc;
        MPI_Comm_rank(global_comm, &world_rank);
        MPI_Comm_size(global_comm, &world_numproc);


        MPI_Comm tmp_comm = MPI_COMM_NULL;
        int      tmp_color = world_rank / N;
        MPI_Comm_split(global_comm, tmp_color, world_rank, &tmp_comm);

        std::vector<int> all_num_neighbors;
        std::vector<int> all_neighbors;
        std::vector<int> displs;
        int tmp_rank, tmp_numproc;
        MPI_Comm_rank(tmp_comm, &tmp_rank);
        MPI_Comm_size(tmp_comm, &tmp_numproc);

        if (tmp_rank == 0) {
            // the local master process
            all_num_neighbors.resize(tmp_numproc, 0);
            displs.resize(tmp_numproc + 1, 0);
        }
        int local_num_neighbors = cs_proc_group.size();
        MPI_Gather(&local_num_neighbors, 1, MPI_INT, all_num_neighbors.data(), 1, MPI_INT, 0, tmp_comm);

        if (tmp_rank == 0) {
            for (int i = 1; i <= tmp_numproc; i++) {
                displs[i] = displs[i - 1] + all_num_neighbors[i - 1];
            }
            all_neighbors.resize(displs[tmp_numproc]);
        }

        std::vector<int> local_neighbors(cs_proc_group.begin(), cs_proc_group.end());
        MPI_Gatherv(local_neighbors.data(), local_num_neighbors, MPI_INT,
            all_neighbors.data(), all_num_neighbors.data(), displs.data(), MPI_INT, 0, tmp_comm);


        using namespace boost;
        using boost::graph::distributed::mpi_process_group;

        typedef adjacency_list<listS, distributedS<mpi_process_group, vecS>, undirectedS> Graph;

        MPI_Comm cp_comm = MPI_COMM_NULL;
        int cp_color = 1;
        if (tmp_rank == 0)  cp_color = 0;
        MPI_Comm_split(global_comm, cp_color, world_rank, &cp_comm);

        boost::mpi::environment  env;
        boost::mpi::communicator cp_world(cp_comm, boost::mpi::comm_take_ownership);

        mpi_process_group mpg(cp_world);

        Graph g(world_numproc, mpg);

        if (tmp_rank != 0) {
            g.clear();
        }

        // construct graph
        if (tmp_rank == 0) {
            // add edges
            for (int i = 0; i < tmp_numproc; i++) {
                // add self
                int cur_id = tmp_color * N + i;
                add_edge(vertex(cur_id, g), vertex(cur_id, g), g);
                for (int j = displs[i]; j < displs[i + 1]; j++) {
                    add_edge(vertex(cur_id, g), vertex(all_neighbors[j], g), g);
                }
            }
        }

        synchronize(g);

        // find connected components
        std::vector<int> local_components_vec(num_vertices(g));
        typedef iterator_property_map<std::vector<int>::iterator, property_map<Graph, vertex_index_t>::type> ComponentMap;
        ComponentMap component(local_components_vec.begin(), get(vertex_index, g));

        int num_components;
        //double start = MPI_Wtime();
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

        int color = 0;
        MPI_Scatter(local_components_vec.data(), 1, MPI_INT, &color, 1, MPI_INT, 0, tmp_comm);
        MPI_Comm_free(&tmp_comm);


        if (cs_comm != MPI_COMM_NULL)   MPI_Comm_free(&cs_comm);
        MPI_Comm_split(global_comm, color, world_rank, &cs_comm);

        int cs_numporc;
        MPI_Comm_size(cs_comm, &cs_numporc);
        std::vector<int> tmpR(cs_numporc);
        MPI_Allgather(&world_rank, 1, MPI_INT, tmpR.data(), 1, MPI_INT, cs_comm);

        cs_proc_group.clear();
        cs_proc_group.insert(tmpR.begin(), tmpR.end());
    }
}
#endif



#ifdef WITH_BGL

static void GroupProcessSeq(std::set<int>& cs_proc_group, MPI_Comm& cs_comm, const MPI_Comm& global_comm)
{
    /// gather the graph to one process first
    {
        int world_rank, world_numproc;
        MPI_Comm_rank(global_comm, &world_rank);
        MPI_Comm_size(global_comm, &world_numproc);

        std::vector<int> all_num_neighbors;
        std::vector<int> all_neighbors;
        std::vector<int> displs;


        if (world_rank == 0) {
            // the local master process
            all_num_neighbors.resize(world_numproc, 0);
            displs.resize(world_numproc + 1, 0);
        }
        int local_num_neighbors = cs_proc_group.size();
        MPI_Gather(&local_num_neighbors, 1, MPI_INT, all_num_neighbors.data(), 1, MPI_INT, 0, global_comm);

        if (world_rank == 0) {
            for (int i = 1; i <= world_numproc; i++) {
                displs[i] = displs[i - 1] + all_num_neighbors[i - 1];
            }
            all_neighbors.resize(displs[world_numproc]);
        }

        std::vector<int> local_neighbors(cs_proc_group.begin(), cs_proc_group.end());
        MPI_Gatherv(local_neighbors.data(), local_num_neighbors, MPI_INT,
            all_neighbors.data(), all_num_neighbors.data(), displs.data(), MPI_INT, 0, global_comm);


        std::vector<int> local_components_vec;

        if (world_rank == 0) {

            using namespace boost;

            typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

            Graph g(world_numproc);

            // add edges
            for (int i = 0; i < world_numproc; i++) {
                // add self
                int cur_id = i;
                add_edge(vertex(cur_id, g), vertex(cur_id, g), g);
                for (int j = displs[i]; j < displs[i + 1]; j++) {
                    add_edge(vertex(cur_id, g), vertex(all_neighbors[j], g), g);
                }
            }

            local_components_vec.resize(num_vertices(g));
            int num = boost::connected_components(g, &local_components_vec[0]);
        }


        int color = 0;
        MPI_Scatter(local_components_vec.data(), 1, MPI_INT, &color, 1, MPI_INT, 0, global_comm);


        if (cs_comm != MPI_COMM_NULL && cs_comm != MPI_COMM_WORLD && cs_comm != MPI_COMM_SELF) {
            MPI_Comm_free(&cs_comm);
        }
        MPI_Comm_split(global_comm, color, world_rank, &cs_comm);

        int cs_numporc;
        MPI_Comm_size(cs_comm, &cs_numporc);
        std::vector<int> tmpR(cs_numporc);
        MPI_Allgather(&world_rank, 1, MPI_INT, tmpR.data(), 1, MPI_INT, cs_comm);

        cs_proc_group.clear();
        cs_proc_group.insert(tmpR.begin(), tmpR.end());
    }
}

#endif


static void GroupProcessSelf(std::set<int>& cs_proc_group, MPI_Comm& cs_comm, const MPI_Comm& global_comm)
{
    int global_rank;
    MPI_Comm_rank(global_comm, &global_rank);
    // every rank is a group
    cs_proc_group.clear();
    cs_proc_group.insert(global_rank);
    cs_comm = MPI_COMM_SELF;
}



void GroupProcess(const OCP_INT& flag, std::set<OCP_INT>& cs_proc_group, MPI_Comm& cs_comm, const MPI_Comm& global_comm)
{
    switch (flag)
    {

#ifdef WITH_PBGL
    case 0:
        GroupProcessPBGL01(cs_proc_group, cs_comm, global_comm);
        break;
    case 1:
        GroupProcessPBGL02(cs_proc_group, cs_comm, global_comm);
        break;
#endif

#ifdef WITH_BGL
    case 2:
        GroupProcessSeq(cs_proc_group, cs_comm, global_comm);
        break;
#endif

    default:
        GroupProcessSelf(cs_proc_group, cs_comm, global_comm);
        break;
    }
}





/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           May/02/2024      Create file                          */
/*----------------------------------------------------------------------------*/