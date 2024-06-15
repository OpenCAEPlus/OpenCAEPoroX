/*! \file    Partition.cpp
 *  \brief   Partition for OpenCAEPoroX simulator
 *  \author  Shizhe Li
 *  \date    Feb/21/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "Partition.hpp"


void Partition::InitMPI(MPI_Comm comm)
{
	myComm = comm;
	MPI_Comm_rank(myComm, &myrank);
	MPI_Comm_size(myComm, &numproc);
}


void Partition::SetPartition(const PreParamGridWell& grid)
{
	if (numproc == 1) {
		numWellTotal    = grid.numWell;
		numElementTotal = grid.activeGridNum + grid.numWell;
		return;
	}
	
	if (CURRENT_RANK == MASTER_PROCESS)
		OCP_INFO("Set Initial Partition -- begin");


	///////////////////////////////////////////////////////
	// Calculate vtxdist
	///////////////////////////////////////////////////////

	vtxdist = new idx_t[numproc + 2]();  // the last element is well num
	if (myrank == MASTER_PROCESS) {
		numWellTotal       = grid.numWell;
		numElementTotal    = grid.activeGridNum + grid.numWell;
		numElementLocal    = numElementTotal / numproc;
		const OCP_INT rest = numElementTotal % numproc;

		vtxdist[0]           = 0;
		vtxdist[numproc]     = numElementTotal;
		vtxdist[numproc + 1] = numWellTotal;
		for (OCP_INT p = 1; p < numproc; p++) {
			vtxdist[p] = vtxdist[p - 1] + numElementLocal;
			if (p <= rest) {
				vtxdist[p]++;
			}
		}
	}

	MPI_Bcast(vtxdist, numproc + 2, IDX_T, MASTER_PROCESS, myComm);
	
	numElementTotal = vtxdist[numproc];
	numWellTotal    = vtxdist[numproc + 1];
	numElementLocal = vtxdist[myrank + 1] - vtxdist[myrank];
	

	///////////////////////////////////////////////////////
	// Calculate xadj adjncy adjwgt
	///////////////////////////////////////////////////////

	if (myrank == MASTER_PROCESS) {

		numEdges    = new idx_t[numproc]();
		maxNumVtx   = 0;
		maxNumEdges = 0;
		for (OCP_INT p = 0; p < numproc; p++) {			
			for (idx_t n = vtxdist[p]; n < vtxdist[p + 1]; n++) {
				numEdges[p] += grid.gNeighbor[n].size();
			}
			maxNumVtx   = max(maxNumVtx, (vtxdist[p + 1] - vtxdist[p]));
			maxNumEdges = max(maxNumEdges, numEdges[p]);			
		}
	}

	MPI_Scatter(numEdges, 1, IDX_T, &numEdgesLocal, 1, IDX_T, MASTER_PROCESS, myComm);
	MPI_Status status;

	if (myrank == MASTER_PROCESS) {

		xadj   = new idx_t[maxNumVtx + 1 + 2 * maxNumEdges]();  // zero is first

		for (OCP_INT p = 1; p < numproc; p++) {
			adjncy = xadj + vtxdist[p + 1] - vtxdist[p] + 1;
			adjwgt = adjncy + numEdges[p];

			for (idx_t n = 0; n < vtxdist[p + 1] - vtxdist[p]; n++) {
				const vector<ConnPair>& gNeigh = grid.gNeighbor[vtxdist[p] + n];
				xadj[n + 1] = xadj[n] + gNeigh.size();
				for (OCP_INT i = 0; i < gNeigh.size(); i++) {
					adjncy[xadj[n] + i] = gNeigh[i].ID();
					adjwgt[xadj[n] + i] = gNeigh[i].WGT();
				}
			}
			
			MPI_Send(xadj, vtxdist[p + 1] - vtxdist[p] + 1 + 2 * numEdges[p], IDX_T, p, 0, myComm);		
		}
		delete[] numEdges;

		// set master process
		adjncy = xadj + numElementLocal + 1;
		adjwgt = adjncy + numEdgesLocal;
		for (idx_t n = 0; n < numElementLocal; n++) {
			const vector<ConnPair>& gNeigh = grid.gNeighbor[n];
			xadj[n + 1] = xadj[n] + gNeigh.size();
			for (OCP_INT i = 0; i < gNeigh.size(); i++) {
				adjncy[xadj[n] + i] = gNeigh[i].ID();
				adjwgt[xadj[n] + i] = gNeigh[i].WGT();
			}
		}
	}
	else {
		xadj   = new idx_t[numElementLocal + 1 + 2 * numEdgesLocal](); // zero is first
		adjncy = xadj + numElementLocal + 1;
		adjwgt = adjncy + numEdgesLocal;
		
		MPI_Recv(xadj, numElementLocal + 1 + 2 * numEdgesLocal, IDX_T, MASTER_PROCESS, 0, myComm, &status);
	}

	if (CURRENT_RANK == MASTER_PROCESS) {
		OCP_INFO("Set Initial Partition -- end");
	}	
	
	CalPartition(grid);
}


void Partition::CalPartition(const PreParamGridWell& grid)
{
	GetWallTime timer;
	timer.Start();

	if (OCP_TRUE) {
		CalPartitionParMetis();
	}
	else {
		CalPartition2D(grid);
	}

	OCPTIME_PARMETIS += timer.Stop();
}


void Partition::CalPartition2D(const PreParamGridWell& grid)
{

	vector<idx_t> tmpPart;
	if (CURRENT_RANK == MASTER_PROCESS) {
		OCP_INFO("2D Partition -- begin");


		if (grid.gridType >= GridType::unstructured) {
			OCP_ABORT("This type of grid can not be partitioned!");
		}

		const USI nx = grid.nx;
		const USI ny = grid.ny;
		const USI nz = grid.nz;

		if ((nx * ny) % numproc != 0) {
			OCP_ABORT("The number of processes must be divisible by nx*ny !");
		}

		USI nPx, nPy;
		{
			OCP_DBL minDiff = 1E20;
			for (USI i = 1; i <= std::sqrt(numproc); i++) {
				if (numproc % i == 0) {
					USI p1 = i;
					USI p2 = numproc / i;

					for (USI j = 0; j < 2; j++) {
						USI px = j == 0 ? p1 : p2;
						USI py = j == 0 ? p2 : p1;

						if (nx % px == 0 && ny % py == 0) {
							OCP_DBL diff = fabs(nx * 1.0 / px - ny * 1.0 / py);

							if (diff < minDiff) {
								minDiff = diff;
								nPx = px;
								nPy = py;
							}
						}
					}
				}
			}
		}

		tmpPart.resize(numElementTotal);
		// for grid
		for (idx_t n = 0; n < numElementTotal - numWellTotal; n++) {

			const OCP_USI tmpID = ((grid.actGC.map_Act2All[n] % (nx * ny * nz)) % (nx * ny));

			const USI pI = (tmpID % nx) / (nx / nPx);
			const USI pJ = (tmpID / nx) / (ny / nPy);
			tmpPart[n] = pJ * nPx + pI;
		}
		// for well
		for (idx_t w = 0; w < numWellTotal; w++){
			// for well
			const OCP_USI tmpID = ((grid.actGC.map_Act2All[grid.connWellGrid[w][0]] % (nx * ny * nz)) % (nx * ny));
		    
		    const USI pI = (tmpID % nx) / (nx / nPx);
		    const USI pJ = (tmpID / nx) / (ny / nPy);
			tmpPart[w + numElementTotal - numWellTotal] = pJ * nPx + pI;
		}	
	}

	vector<OCP_INT> send_counts;
	vector<OCP_INT> displs;
	if (CURRENT_RANK == MASTER_PROCESS) {

		displs.resize(numproc + 1);
		copy(vtxdist, vtxdist + numproc + 1, displs.begin());
		send_counts.resize(numproc);
		for (OCP_USI p = 0; p < numproc; p++) {
			send_counts[p] = displs[p + 1] - displs[p];
		}
	}

	part = new idx_t[numElementLocal]();
	MPI_Scatterv(tmpPart.data(), send_counts.data(), displs.data(), IDX_T,
		         part, numElementLocal, IDX_T, MASTER_PROCESS, myComm);

	if (CURRENT_RANK == MASTER_PROCESS) {
		OCP_INFO("2D Partition -- end");
	}
}


void Partition::CalPartitionParMetis()
{
	if (CURRENT_RANK == MASTER_PROCESS) {
		OCP_INFO("ParMetis Partition -- begin");
	}

	InitParam();
	ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, &ncon,
		&nparts, tpwgts, &ubvec, options, &edgecut, part, &myComm);

	delete[] tpwgts;
	delete[] options;

	if (CURRENT_RANK == MASTER_PROCESS) {
		OCP_INFO("ParMetis Partition -- end");
	}
}


void Partition::InitParam()
{
	vwgt    = nullptr;
	wgtflag = 1;
	numflag = 0;
	ncon    = 1;
	nparts  = numproc;
	tpwgts  = new real_t[ncon * numproc]();
	for (idx_t i = 0; i < ncon * numproc; i++)
		tpwgts[i] = 1.0 / numproc;
	ubvec   = 1.05;
	options = new idx_t[3]();
	part    = new idx_t[numElementLocal]();
}


void Partition::SetDistribution()
{
	if (numproc == 1)  return;

	if (CURRENT_RANK == MASTER_PROCESS) {
		OCP_INFO("Set Distribution -- begin");
	}

	// now from partition to distribution
	// Setup elementCSR
	
	// Exchange neighbor's process first
	
	// right_neighbor stores all elements distributed in other proc now(may have duplicate elements)
	// which are ranged in ascending order with its initial global index.
	// right_neighbor_loc stores corresponding location within adjncy.
	// right_neighbor_proc stores the corresponding proc id, which will recive from other proc.
	// left_neighbor_proc stores the proc id of left_neighbor, arranged in ascending order with its global index
	// right_neighbor_proc will recive elements from proc right_neighbor_proc[0]'s left_neighbor_proc
	// left_neighbor_proc will send elements to proc left_neighbor_proc[0]'s right_neighbor_proc

	idx_t* adjproc = new idx_t[numEdgesLocal]();

	/// Get neighbors' process from other process
	vector<vector<idx_t>>     left_neighbor_proc;
	vector<vector<NeighborP>> right_neighbor;
	vector<vector<idx_t>>     right_neighbor_proc;

	for (idx_t i = 0; i < numElementLocal; i++) {
		for (idx_t j = xadj[i]; j < xadj[i + 1]; j++) {
			if (adjncy[j] < vtxdist[myrank] || adjncy[j] >= vtxdist[myrank + 1]) {
				// adjncy[j] is not in current proc, find its process
				idx_t p = 0;
				if (adjncy[j] >= vtxdist[p]) {
					while (adjncy[j] >= vtxdist[p + 1]) { p++; }
				}
				// adjncy[j] is in proc p now (p != global_rank)
				idx_t k = 0;
				for (k = 0; k < right_neighbor.size(); k++) {
					// indentify storage location
					if (right_neighbor[k][0].index == p) {
						// existing neighboe proc
						right_neighbor[k].push_back(NeighborP(adjncy[j], j));
						left_neighbor_proc[k].push_back(part[i]);
						break;
					}
				}
				if (k == right_neighbor.size()) {
					// new neighbor proc
					right_neighbor.push_back(vector<NeighborP>{NeighborP(p, p), NeighborP(adjncy[j], j)});
					left_neighbor_proc.push_back(vector<idx_t> {p, part[i]});
				}
			}
			else {
				// adjncy[j] is in current proc
				adjproc[j] = part[adjncy[j] - vtxdist[myrank]];
			}
		}
	}

	// Now allocate for right_neighbor_proc and reordering right_neighbor
	right_neighbor_proc.resize(right_neighbor.size());
	for (idx_t i = 0; i < right_neighbor.size(); i++) {
		right_neighbor_proc[i].resize(right_neighbor[i].size());
		right_neighbor_proc[i][0] = right_neighbor[i][0].index;
		sort(right_neighbor[i].begin() + 1, right_neighbor[i].end(), NeighborP::less);
	}

	MPI_Request request;
	MPI_Status  status;
	// Communicate to get neighbor elements' process
	for (auto& s : left_neighbor_proc) {
		// cout << "First stage:  " << global_rank << " send " << s.size() - 1 << "s to " << s[0] << endl;
		MPI_Isend(s.data() + 1, s.size() - 1, IDX_T, s[0], 0, myComm, &request);
	}
	for (auto& r : right_neighbor_proc) {
		// cout << "First stage:  " << global_rank << " recv " << r.size() - 1 << "s from " << r[0] << endl;
		MPI_Recv(r.data() + 1, r.size() - 1, IDX_T, r[0], 0, myComm, &status);
	}

	// Calculate adjncyProc
	for (idx_t n = 0; n < right_neighbor.size(); n++) {
		for (idx_t i = 1; i < right_neighbor[n].size(); i++) {
			adjproc[right_neighbor[n][i].location] = right_neighbor_proc[n][i];
		}
	}

	// Free memory
	vector<vector<NeighborP>>().swap(right_neighbor);
	vector<vector<idx_t>>().swap(right_neighbor_proc);


	/// Get elements and their neighbors belonging to current process,
	/// Also, process of them are needed.
	/// Setup elements those will be sent to corresponding process
	//  target process; num of elements, num of edges, global index of elements, xadj, adjncy, adjproc
	map<OCP_INT, vector<idx_t>> send_buffer;
	auto& recv_buffer = elementCSR;
	recv_buffer[myrank] = vector<idx_t>{ 0,0 };

	// Calculate necessary memory first
	for (idx_t i = 0; i < numElementLocal; i++) {
		if (part[i] == myrank) {
			recv_buffer[myrank][0]++;
			recv_buffer[myrank][1] += (xadj[i + 1] - xadj[i]);
		}
		else {
			if (!send_buffer.count(part[i])) {
				send_buffer[part[i]] = vector<idx_t>{ 1, xadj[i + 1] - xadj[i] };
			}
			else {
				send_buffer[part[i]][0]++;
				send_buffer[part[i]][1] += xadj[i + 1] - xadj[i];
			}
		}
	}

	// Allocate send_buffer and its pointer
	// global index of elements, xadj, adjncy, adjproc
	map<OCP_INT, vector<idx_t>> send_buffer_ptr;
	for (auto& s : send_buffer) {
		auto& sb = s.second;
		sb.resize(2 + sb[0] + (sb[0] + 1) + sb[1] + sb[1]);
		send_buffer_ptr[s.first].push_back(2);
		send_buffer_ptr[s.first].push_back(2 + sb[0]);
		send_buffer_ptr[s.first].push_back(2 + sb[0] + (sb[0] + 1));
		send_buffer_ptr[s.first].push_back(2 + sb[0] + (sb[0] + 1) + sb[1]);
	}
	// Allocate recv_buffer for self and its pointer
	// global index of elements, xadj, adjncy, adjproc
	vector<idx_t> recv_buffer_ptr(4);
	auto& rb = recv_buffer[myrank];
	rb.resize(2 + rb[0] + (rb[0] + 1) + rb[1] + rb[1]);
	recv_buffer_ptr[0] = 2;
	recv_buffer_ptr[1] = 2 + rb[0];
	recv_buffer_ptr[2] = 2 + rb[0] + (rb[0] + 1);
	recv_buffer_ptr[3] = 2 + rb[0] + (rb[0] + 1) + rb[1];

	// Setup send_buffer and recv_buffer
	for (idx_t i = 0; i < numElementLocal; i++) {

		if (part[i] == myrank) {
			auto& rb = recv_buffer[part[i]];

			rb[recv_buffer_ptr[0]++] = i + vtxdist[myrank];
			rb[recv_buffer_ptr[1]]   = rb[recv_buffer_ptr[1]++] + xadj[i + 1] - xadj[i];

			copy(&adjncy[xadj[i]], &adjncy[xadj[i + 1]], &rb[recv_buffer_ptr[2]]);
			copy(&adjproc[xadj[i]], &adjproc[xadj[i + 1]], &rb[recv_buffer_ptr[3]]);
			recv_buffer_ptr[2] += xadj[i + 1] - xadj[i];
			recv_buffer_ptr[3] += xadj[i + 1] - xadj[i];
		}
		else {
			auto& sb  = send_buffer[part[i]];
			auto& sbp = send_buffer_ptr[part[i]];

			sb[sbp[0]++] = i + vtxdist[myrank];
			sb[sbp[1]]   = sb[sbp[1]++] + xadj[i + 1] - xadj[i];

			copy(&adjncy[xadj[i]], &adjncy[xadj[i + 1]], &sb[sbp[2]]);
			copy(&adjproc[xadj[i]], &adjproc[xadj[i + 1]], &sb[sbp[3]]);
			sbp[2] += xadj[i + 1] - xadj[i];
			sbp[3] += xadj[i + 1] - xadj[i];
		}
	}

	// Free Memory
	map<OCP_INT, vector<idx_t>>().swap(send_buffer_ptr);
	vector<idx_t>().swap(recv_buffer_ptr);
	delete[] xadj;
	delete[] vtxdist;
	delete[] adjproc;
	delete[] part;


	vector<idx_t> send_size(numproc, 0);
	vector<idx_t> recv_size(numproc, 0);
	for (const auto& s : send_buffer) {
		send_size[s.first] = s.second.size();
	}

	MPI_Alltoall(&send_size[0], 1, IDX_T, &recv_size[0], 1, IDX_T, myComm);

	// Allocate for recv_buffer
	for (OCP_INT i = 0; i < numproc; i++) {
		if (recv_size[i] > 0) {
			recv_buffer[i].resize(recv_size[i]);
		}
	}

	// Free Memory
	vector<idx_t>().swap(send_size);
	vector<idx_t>().swap(recv_size);


	// Communicate for buffer
	for (auto& s : send_buffer) {
		// cout << "Second stage:  " << global_rank << " send " << s.size() - 1 << "s to " << s[0] << endl;
		MPI_Isend(s.second.data(), s.second.size(), IDX_T, s.first, 0, myComm, &request);
	}
	for (auto& r : recv_buffer) {
		if (r.first == myrank)  continue;
		// cout << "Second stage:  " << global_rank << " recv " << r.size() - 1 << "s from " << r[0] << endl;
		MPI_Recv(r.second.data(), r.second.size(), IDX_T, r.first, 0, myComm, &status);
	}

	MPI_Barrier(myComm);


	if (CURRENT_RANK == MASTER_PROCESS) {
		OCP_INFO("Set Distribution -- end");
	}
}


 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Feb/21/2023      Create file                          */
 /*----------------------------------------------------------------------------*/