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

/// Average number of processes per process group when partitioning with METIS
static const int METIS_AVNP = 50;
static real_t ubvec = 1.01;

#ifdef WITH_PBGL

 // each process with one vertex
static int GroupProcessPBGL01(const std::set<int>& cs_proc_group, const MPI_Comm& work_comm)
{

	int work_rank, work_numproc;
	MPI_Comm_rank(work_comm, &work_rank);
	MPI_Comm_size(work_comm, &work_numproc);

	int color;
	{
		using namespace boost;
		using boost::graph::distributed::mpi_process_group;

		typedef adjacency_list<listS, distributedS<mpi_process_group, vecS>, undirectedS> Graph;

		// boost::mpi::environment  env;
		// boost::mpi::communicator cp_world(work_comm, boost::mpi::comm_attach);
		// mpi_process_group mpg(cp_world);

		// std::cout << work_rank << "*** " << mpg.size << "   " << mpg.rank << std::endl;


		Graph g(work_numproc);

		for (const auto& n : cs_proc_group) {
			add_edge(vertex(work_rank, g), vertex(n, g), g);
		}
		add_edge(vertex(work_rank, g), vertex(work_rank, g), g);

		synchronize(g);

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
		//if (work_rank == 0) {
		//    std::cout << num_components << " components have been found! -- "
		//        << MPI_Wtime() - start << "s" << std::endl;
		//}

		color = local_components_vec[0];
	}

	return color;
}

// gather version
static int GroupProcessPBGL02(const std::set<int>& cs_proc_group, const MPI_Comm& work_comm)
{
	/// number of vertex in each process, remember that numproc = numvertex
	const int N = 1000;


	int work_rank, work_numproc;
	MPI_Comm_rank(work_comm, &work_rank);
	MPI_Comm_size(work_comm, &work_numproc);


	MPI_Comm tmp_comm = MPI_COMM_NULL;
	int      tmp_color = work_rank / N;
	MPI_Comm_split(work_comm, tmp_color, work_rank, &tmp_comm);

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


	std::vector<int> local_components_vec;
	{
		using namespace boost;
		using boost::graph::distributed::mpi_process_group;

		typedef adjacency_list<listS, distributedS<mpi_process_group, vecS>, undirectedS> Graph;

		MPI_Comm cp_comm = MPI_COMM_NULL;
		int cp_color = 1;
		if (tmp_rank == 0)  cp_color = 0;
		MPI_Comm_split(work_comm, cp_color, work_rank, &cp_comm);

		boost::mpi::environment  env;
		boost::mpi::communicator cp_world(cp_comm, boost::mpi::comm_take_ownership);

		mpi_process_group mpg(cp_world);

		//std::cout << tmp_rank << "*** " << mpg.size << "   " << mpg.rank << std::endl;

		Graph g(work_numproc, mpg);

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
			synchronize(g);
		}

		//std::cout << tmp_rank << "*** "
		//	<< g.process_group().size << "   " << g.process_group().rank << "  " << num_vertices(g) << std::endl;

		// find connected components
		local_components_vec.resize(num_vertices(g));
		typedef iterator_property_map<std::vector<int>::iterator, property_map<Graph, vertex_index_t>::type> ComponentMap;
		ComponentMap component(local_components_vec.begin(), get(vertex_index, g));

		//double start = MPI_Wtime();
		int num_components;
		if (tmp_rank == 0) {
			if (false) {
				num_components = connected_components_ps(g, component);
			}
			else {
				num_components = connected_components(g, component);
			}

		}

		/*if (work_rank == 0) {
			std::cout << num_components << " components have been found! -- "
				<< MPI_Wtime() - start << "s" << std::endl;
		}*/


	}

	int color;
	MPI_Scatter(local_components_vec.data(), 1, MPI_INT, &color, 1, MPI_INT, 0, tmp_comm);
	MPI_Comm_free(&tmp_comm);

	return color;
}

#endif


#ifdef WITH_BGL
/// gather the graph to one process
static int GroupProcessSeq(const std::set<int>& cs_proc_group, const MPI_Comm& work_comm)
{

	int work_rank, work_numproc;
	MPI_Comm_rank(work_comm, &work_rank);
	MPI_Comm_size(work_comm, &work_numproc);

	std::vector<int> all_num_neighbors;
	std::vector<int> all_neighbors;
	std::vector<int> displs;


	if (work_rank == 0) {
		// the local master process
		all_num_neighbors.resize(work_numproc, 0);
		displs.resize(work_numproc + 1, 0);
	}
	int local_num_neighbors = cs_proc_group.size();
	MPI_Gather(&local_num_neighbors, 1, MPI_INT, all_num_neighbors.data(), 1, MPI_INT, 0, work_comm);

	if (work_rank == 0) {
		for (int i = 1; i <= work_numproc; i++) {
			displs[i] = displs[i - 1] + all_num_neighbors[i - 1];
		}
		all_neighbors.resize(displs[work_numproc]);
	}

	std::vector<int> local_neighbors(cs_proc_group.begin(), cs_proc_group.end());
	MPI_Gatherv(local_neighbors.data(), local_num_neighbors, MPI_INT,
		all_neighbors.data(), all_num_neighbors.data(), displs.data(), MPI_INT, 0, work_comm);


	std::vector<int> local_components_vec;
	{
		if (work_rank == 0) {

			using namespace boost;
			typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

			Graph g(work_numproc);

			// add edges
			for (int i = 0; i < work_numproc; i++) {
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
	}

	int color = 0;
	MPI_Scatter(local_components_vec.data(), 1, MPI_INT, &color, 1, MPI_INT, 0, work_comm);

	return color;
}

#endif


#ifdef WITH_PARMETIS

// each process with one vertex
static int GroupProcessParMetis01(const std::unordered_map<OCP_INT, OCP_INT>& proc_weight,
	const MPI_Comm& work_comm)
{

	int work_rank, work_numproc;
	MPI_Comm_rank(work_comm, &work_rank);
	MPI_Comm_size(work_comm, &work_numproc);

	const int NPARTS = (work_numproc - 1) / METIS_AVNP + 1;

	int color = 0;
	{
		// use Metis
		std::vector<idx_t> xadj(2, 0);
		std::vector<idx_t> adjncy;
		std::vector<idx_t> adjwgt;

		xadj[1] = proc_weight.size();
		int i = 0;
		for (const auto& n : proc_weight) {
			adjncy.push_back(n.first);
			adjwgt.push_back(n.second);
			i++;
		}

		std::vector<idx_t> vtxdist(work_numproc + 1, work_numproc);
		for (int i = 0; i < work_numproc; i++) {
			vtxdist[i] = i;
		}

		// print for check
		//if (false)
		//{
		//	for (int i = 0; i < work_numproc; i++) {
		//		if (i == work_rank) {
		//			std::this_thread::sleep_for(std::chrono::milliseconds(i * 500));
		//		}
		//	}

		//	std::cout << "****Rank " << work_rank << std::endl;
		//	for (int i = 0; i < 2; i++) {
		//		std::cout << xadj[i] << "   ";
		//	}
		//	std::cout << std::endl;

		//	for (int j = xadj[0]; j < xadj[1]; j++) {
		//		std::cout << "( " << adjncy[j] << ", " << adjwgt[j] << " )  ";
		//	}
		//	std::cout << std::endl;
		//}

		idx_t nparts = NPARTS;
		std::vector<idx_t> parts(1);
		idx_t* vwgt = nullptr;
		idx_t wgtflag = 1;
		idx_t numflag = 0;
		idx_t ncon = 1;
		idx_t edgecut;

		std::vector<real_t> tpwgts(ncon * nparts, 0);
		for (idx_t i = 0; i < ncon * nparts; i++)
			tpwgts[i] = 1.0 / nparts;
		std::vector<idx_t> options(3, 0);

		MPI_Comm tmp_comm{ MPI_COMM_NULL };
		MPI_Comm_dup(work_comm, &tmp_comm);

		ParMETIS_V3_PartKway(vtxdist.data(), xadj.data(), adjncy.data(), vwgt, adjwgt.data(), &wgtflag, &numflag, &ncon,
			&nparts, tpwgts.data(), &ubvec, options.data(), &edgecut, parts.data(), &tmp_comm);

		MPI_Comm_free(&tmp_comm);

		color = parts[0];
	}

	return color;
}

// gather version
static int GroupProcessParMetis02(const std::unordered_map<OCP_INT, OCP_INT>& proc_weight,
	const MPI_Comm& work_comm)
{
	/// number of vertex in each process, remember that numproc = numvertex
	const int N = 1000;

	int work_rank, work_numproc;
	MPI_Comm_rank(work_comm, &work_rank);
	MPI_Comm_size(work_comm, &work_numproc);

	const int NPARTS = (work_numproc - 1) / METIS_AVNP + 1;


	MPI_Comm tmp_comm = MPI_COMM_NULL;
	int      tmp_color = work_rank / N;
	MPI_Comm_split(work_comm, tmp_color, work_rank, &tmp_comm);

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
	int local_num_neighbors = proc_weight.size() * 2;
	MPI_Gather(&local_num_neighbors, 1, MPI_INT, all_num_neighbors.data(), 1, MPI_INT, 0, tmp_comm);

	if (tmp_rank == 0) {
		for (int i = 1; i <= tmp_numproc; i++) {
			displs[i] = displs[i - 1] + all_num_neighbors[i - 1];
		}
		all_neighbors.resize(displs[tmp_numproc]);
	}

	std::vector<int> local_neighbors(local_num_neighbors);
	{
		int num = proc_weight.size();
		int i = 0;
		for (const auto& n : proc_weight) {
			local_neighbors[i] = n.first;
			local_neighbors[i + num] = n.second;
			i++;
		}
	}

	MPI_Gatherv(local_neighbors.data(), local_num_neighbors, MPI_INT,
		all_neighbors.data(), all_num_neighbors.data(), displs.data(), MPI_INT, 0, tmp_comm);


	std::vector<int> parts;


	MPI_Comm cp_comm = MPI_COMM_NULL;
	int cp_color = 1;
	if (tmp_rank == 0)  cp_color = 0;
	MPI_Comm_split(work_comm, cp_color, work_rank, &cp_comm);

	if (tmp_rank == 0) {

		// use Metis
		std::vector<idx_t> xadj(tmp_numproc + 1, 0);
		std::vector<idx_t> adjncy(displs[tmp_numproc] / 2);
		std::vector<idx_t> adjwgt(displs[tmp_numproc] / 2);

		for (int i = 0; i < tmp_numproc; i++) {
			xadj[i + 1] = xadj[i] + all_num_neighbors[i] / 2;

			int* myptr = all_neighbors.data() + displs[i];
			std::copy(myptr, myptr + all_num_neighbors[i] / 2, &adjncy[xadj[i]]);

			myptr += all_num_neighbors[i] / 2;
			std::copy(myptr, myptr + all_num_neighbors[i] / 2, &adjwgt[xadj[i]]);
		}

		const int num_use_process = (work_numproc - 1) / N + 1;
		std::vector<idx_t> vtxdist(num_use_process + 1, 0);

		for (int i = 1; i < num_use_process; i++) {
			vtxdist[i] = vtxdist[i - 1] + N;
		}
		vtxdist[num_use_process] = work_numproc;

		// print for check
		//if (false)
		//{
		//	int cp_rank, cp_numproc;;
		//	MPI_Comm_rank(cp_comm, &cp_rank);
		//	MPI_Comm_size(cp_comm, &cp_numproc);

		//	for (int i = 0; i < cp_numproc; i++) {
		//		if (i == cp_rank) {
		//			std::this_thread::sleep_for(std::chrono::milliseconds(i * 500));
		//		}
		//	}

		//	std::cout << "****Rank " << cp_rank << std::endl;

		//	for (int i = 0; i < tmp_numproc + 1; i++) {
		//		std::cout << xadj[i] << "   ";
		//	}
		//	std::cout << std::endl;

		//	for (int i = 0; i < tmp_numproc; i++) {
		//		for (int j = xadj[i]; j < xadj[i + 1]; j++) {
		//			std::cout << "( " << adjncy[j] << ", " << adjwgt[j] << " )  ";
		//		}
		//		std::cout << std::endl;
		//	}
		//}

		idx_t nparts = NPARTS;
		std::vector<idx_t> tmp_parts(tmp_numproc);
		idx_t* vwgt = nullptr;
		idx_t wgtflag = 1;
		idx_t numflag = 0;
		idx_t ncon = 1;
		idx_t edgecut;

		real_t* tpwgts = new real_t[ncon * nparts]();
		for (idx_t i = 0; i < ncon * nparts; i++)
			tpwgts[i] = 1.0 / nparts;
		idx_t* options = new idx_t[3]();

		ParMETIS_V3_PartKway(vtxdist.data(), xadj.data(), adjncy.data(), vwgt, adjwgt.data(), &wgtflag, &numflag, &ncon,
			&nparts, tpwgts, &ubvec, options, &edgecut, tmp_parts.data(), &cp_comm);


		MPI_Comm_free(&cp_comm);
		parts.assign(tmp_parts.begin(), tmp_parts.end());
	}


	int color = 0;
	MPI_Scatter(parts.data(), 1, MPI_INT, &color, 1, MPI_INT, 0, tmp_comm);
	MPI_Comm_free(&tmp_comm);

	return color;
}

#endif


#ifdef WITH_METIS

void MyMetis(std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy, std::vector<idx_t>& adjwgt,
	idx_t& nParts, std::vector<idx_t>& part, decltype(METIS_PartGraphKway)* METIS_PartGraphFunc)
{
	idx_t nVertices = xadj.size() - 1;
	// idx_t nEdges    = adjncy.size() / 2;
	idx_t nWeights = 1;
	idx_t objval;

	int ret = METIS_PartGraphFunc(&nVertices, &nWeights, xadj.data(), adjncy.data(),
		NULL, NULL, adjwgt.data(), &nParts, NULL, &ubvec, NULL, &objval, part.data());

	if (ret != rstatus_et::METIS_OK) { std::cout << "METIS_ERROR" << std::endl; }
	//std::cout << "METIS_OK" << std::endl;
	//std::cout << "objval: " << objval << std::endl;
}

/// gather the graph to one process
static int GroupProcessMetis(const std::unordered_map<OCP_INT, OCP_INT>& proc_weight,
	const MPI_Comm& work_comm)
{
	int work_rank, work_numproc;
	MPI_Comm_rank(work_comm, &work_rank);
	MPI_Comm_size(work_comm, &work_numproc);

	//if (work_rank == 0) std::cout << "work_numproc: " << work_numproc << std::endl;

	const int NPARTS = (work_numproc - 1) / METIS_AVNP + 1;

	std::vector<int> all_num_neighbors;
	std::vector<int> all_neighbors;
	std::vector<int> displs;

	if (work_rank == 0) {
		// the local master process
		all_num_neighbors.resize(work_numproc, 0);
		displs.resize(work_numproc + 1, 0);
	}
	int local_num_neighbors = proc_weight.size() * 2;
	MPI_Gather(&local_num_neighbors, 1, MPI_INT, all_num_neighbors.data(), 1, MPI_INT, 0, work_comm);

	if (work_rank == 0) {
		for (int i = 1; i <= work_numproc; i++) {
			displs[i] = displs[i - 1] + all_num_neighbors[i - 1];
		}
		all_neighbors.resize(displs[work_numproc]);
	}

	std::vector<int> local_neighbors(local_num_neighbors);
	{
		int num = proc_weight.size();
		int i = 0;
		for (const auto& n : proc_weight) {
			local_neighbors[i] = n.first;
			local_neighbors[i + num] = n.second;
			i++;
		}
	}

	MPI_Gatherv(local_neighbors.data(), local_num_neighbors, MPI_INT,
		all_neighbors.data(), all_num_neighbors.data(), displs.data(), MPI_INT, 0, work_comm);



	std::vector<int> parts;

	{
		if (work_rank == 0) {
			// use Metis
			std::vector<idx_t> xadj(work_numproc + 1, 0);
			std::vector<idx_t> adjncy(displs[work_numproc] / 2);
			std::vector<idx_t> adjwgt(displs[work_numproc] / 2);

			for (int i = 0; i < work_numproc; i++) {
				xadj[i + 1] = xadj[i] + all_num_neighbors[i] / 2;

				int* myptr = all_neighbors.data() + displs[i];
				std::copy(myptr, myptr + all_num_neighbors[i] / 2, &adjncy[xadj[i]]);

				myptr += all_num_neighbors[i] / 2;
				std::copy(myptr, myptr + all_num_neighbors[i] / 2, &adjwgt[xadj[i]]);
			}


			// print for check
	/*            for (int i = 0; i < work_numproc + 1; i++) {
					std::cout << xadj[i] << "   ";
				}
				std::cout << std::endl;

				for (int i = 0; i < work_numproc; i++) {
					for (int j = xadj[i]; j < xadj[i + 1]; j++) {
						std::cout << "( " << adjncy[j] << ", " << adjwgt[j] << " )  ";
					}
					std::cout << std::endl;
				}    */


			std::vector<idx_t> tmp_parts(work_numproc);
			idx_t nparts = NPARTS;

			if (nparts > 8 || true) {
				MyMetis(xadj, adjncy, adjwgt, nparts, tmp_parts, METIS_PartGraphKway);
			}
			else {
				MyMetis(xadj, adjncy, adjwgt, nparts, tmp_parts, METIS_PartGraphRecursive);
			}

			parts.assign(tmp_parts.begin(), tmp_parts.end());
		}
	}


	int color = 0;
	MPI_Scatter(parts.data(), 1, MPI_INT, &color, 1, MPI_INT, 0, work_comm);

	return color;
}

#endif


static void GroupProcessSelf(std::set<int>& cs_proc_group, MPI_Comm& cs_comm, const MPI_Comm& global_comm)
{
	int global_rank;
	MPI_Comm_rank(global_comm, &global_rank);
	// every rank is a group
	cs_proc_group.clear();
	cs_proc_group.insert(global_rank);

	if (cs_comm != MPI_COMM_NULL)  MPI_Comm_free(&cs_comm);
	MPI_Comm_dup(MPI_COMM_SELF, &cs_comm);
}


static void CalCommGroup(const int& color, const int& global_rank, const MPI_Comm& work_comm, std::set<int>& cs_proc_group, MPI_Comm& cs_comm)
{

	if (cs_comm != MPI_COMM_NULL)  MPI_Comm_free(&cs_comm);
	MPI_Comm_split(work_comm, color, global_rank, &cs_comm);

	int cs_numporc;
	MPI_Comm_size(cs_comm, &cs_numporc);
	std::vector<int> tmpR(cs_numporc);
	MPI_Allgather(&global_rank, 1, MPI_INT, tmpR.data(), 1, MPI_INT, cs_comm);

	cs_proc_group.clear();
	cs_proc_group.insert(tmpR.begin(), tmpR.end());
}


// Since isolated points are removed first, only a subset of processes needs to participate in the partitioning
// algorithm, which speeds up the partitioning computation and also avoids null adjncy (which could cause errors) 
// in the ParMetis method. On the other hand, for PBGL, computation fails if only some processes are involved, 
// as it seems to require synchronization across all processes.
void GroupProcess(const GroupMethod& method, std::set<OCP_INT>& cs_proc_group,
	const std::unordered_map<OCP_INT, OCP_INT>& proc_weight,
	MPI_Comm& cs_comm, const MPI_Comm& global_comm)
{

	// remove the isolated points 
	int global_rank;
	MPI_Comm_rank(global_comm, &global_rank);

	GroupMethod flag = method;
	OCP_INT     color = 1;
	if (cs_proc_group.size() == 0 && proc_weight.size() == 0) {
		flag = GroupMethod::Self;
		color = 0;
	}
	else if (cs_proc_group.size() == 0 && proc_weight.size() > 0) {
		// use cs_proc_group's info to fill cs_proc_group if cs_proc_group is not given
		cs_proc_group.clear();
		for (const auto& p : proc_weight)  cs_proc_group.insert(p.first);
	}

	MPI_Comm work_comm = MPI_COMM_NULL;
	MPI_Comm_split(global_comm, color, global_rank, &work_comm);
	int work_rank;
	MPI_Comm_rank(work_comm, &work_rank);
	// recalulate index of points on each process

	std::set<OCP_INT>                    work_cs_proc_group;
	std::unordered_map<OCP_INT, OCP_INT> work_proc_weight;


	if (false) {

		MPI_Request work_request;
		MPI_Status  work_status;

		for (const auto& s : cs_proc_group) {
			MPI_Isend(&work_rank, 1, MPI_INT, s, 0, global_comm, &work_request);
		}
		for (const auto& r : cs_proc_group) {
			OCP_INT buffer = 0;
			MPI_Recv(&buffer, 1, MPI_INT, r, 0, global_comm, &work_status);

			work_cs_proc_group.insert(buffer);
			work_proc_weight[buffer] = proc_weight.at(r);
		}
	}
	else {

		const USI len = cs_proc_group.size();
		std::vector<MPI_Request>  send_request(len);
		std::vector<MPI_Request>  recv_request(len);
		std::vector<OCP_INT>      buffer(len);

		USI iter = 0;
		for (const auto& r : cs_proc_group) {
			MPI_Irecv(&buffer[iter], 1, MPI_INT, r, 0, global_comm, &recv_request[iter]);
			iter++;
		}

		iter = 0;
		for (const auto& s : cs_proc_group) {
			MPI_Isend(&work_rank, 1, MPI_INT, s, 0, global_comm, &send_request[iter]);
			iter++;
		}

		MPI_Waitall(iter, recv_request.data(), MPI_STATUS_IGNORE);
		MPI_Waitall(iter, send_request.data(), MPI_STATUS_IGNORE);

		iter = 0;
		for (const auto& r : cs_proc_group) {
			work_cs_proc_group.insert(buffer[iter]);
			work_proc_weight[buffer[iter]] = proc_weight.at(r);
			iter++;
		}
	}

	color = 0;

	switch (flag)
	{
#ifdef WITH_PBGL
	case GroupMethod::PBGL:
		color = GroupProcessPBGL01(work_cs_proc_group, work_comm);
		break;
	case GroupMethod::PBGL_gather:
		color = GroupProcessPBGL02(work_cs_proc_group, work_comm);
		break;
#endif

#ifdef WITH_BGL
	case GroupMethod::BGL:
		color = GroupProcessSeq(work_cs_proc_group, work_comm);
		break;
#endif

#ifdef WITH_PARMETIS
	case GroupMethod::ParMetis:
		color = GroupProcessParMetis01(work_proc_weight, work_comm);
		break;
	case GroupMethod::ParMetis_gather:
		color = GroupProcessParMetis02(work_proc_weight, work_comm);
		break;
#endif

#ifdef WITH_METIS
	case GroupMethod::Metis:
		color = GroupProcessMetis(work_proc_weight, work_comm);
		break;
#endif

	case GroupMethod::Self:
	default:
		GroupProcessSelf(cs_proc_group, cs_comm, global_comm);
		break;
	}


	if (flag != GroupMethod::Self) {
		CalCommGroup(color, global_rank, work_comm, cs_proc_group, cs_comm);
	}

	MPI_Comm_free(&work_comm);
}





/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           May/02/2024      Create file                          */
/*----------------------------------------------------------------------------*/