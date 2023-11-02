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
				numEdges[p] += grid.numNeighbor[n];				
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
		OCP_INFO("ParMetis Partition -- begin");
	}
		


	GetWallTime timer;
	timer.Start();
	
	InitParam();
	ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, &ncon,
		&nparts, tpwgts, &ubvec, options, &edgecut, part, &myComm);

	OCPTIME_PARMETIS += timer.Stop() / TIME_S2MS;

	//////////////////////////////////////////////////////
	// Free Memory
	delete[] tpwgts;
	delete[] options;
	// Free Memory
	//////////////////////////////////////////////////////


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
				// estimate first
				idx_t p = adjncy[j] / numElementLocal;  // p must be an signed int
				if (adjncy[j] >= vtxdist[p]) {
					while (adjncy[j] >= vtxdist[p + 1]) { p++; }
				}
				else {
					while (adjncy[j] < vtxdist[--p]) {}
				}
				// adjncy[j] is in proc p now (p != myrank)
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
		// cout << "First stage:  " << myrank << " send " << s.size() - 1 << "s to " << s[0] << endl;
		MPI_Isend(s.data() + 1, s.size() - 1, IDX_T, s[0], 0, myComm, &request);
	}
	for (auto& r : right_neighbor_proc) {
		// cout << "First stage:  " << myrank << " recv " << r.size() - 1 << "s from " << r[0] << endl;
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
	//  target process, num of elements, num of edges, global index of elements, xadj, adjncy, adjproc
	vector<vector<idx_t>>  send_buffer;
	vector<vector<idx_t>>& recv_buffer = elementCSR;
	recv_buffer.push_back(vector<idx_t> {myrank, 0, 0});

	// Calculate necessary memory first
	for (idx_t i = 0; i < numElementLocal; i++) {
		if (part[i] == myrank) {
			recv_buffer[0][1] ++;
			recv_buffer[0][2] += (xadj[i + 1] - xadj[i]);
		}
		else {
			idx_t k;
			for (k = 0; k < send_buffer.size(); k++) {
				if (part[i] == send_buffer[k][0]) {
					// existing process
					send_buffer[k][1] ++;
					send_buffer[k][2] += (xadj[i + 1] - xadj[i]);
					break;
				}
			}
			if (k == send_buffer.size()) {
				// new process
				send_buffer.push_back(vector<idx_t>{part[i], 1, xadj[i + 1] - xadj[i]});
			}
		}
	}

	// Allocate send_buffer and its pointer
	vector<idx_t> send_buffer_ptr(send_buffer.size() * 4);
	idx_t         pIter = 0;
	for (auto& s : send_buffer) {
		s.resize(4 + 2 * (s[1] + s[2]));
		send_buffer_ptr[pIter++] = 3;
		send_buffer_ptr[pIter++] = 3 + s[1];
		send_buffer_ptr[pIter++] = 4 + s[1] + s[1];
		send_buffer_ptr[pIter++] = 4 + s[1] + s[1] + s[2];
	}
	// Allocate recv_buffer for self and its pointer
	// size is at most one now
	vector<idx_t> recv_buffer_ptr(recv_buffer.size() * 4);
	pIter = 0;
	for (auto& r : recv_buffer) {
		r.resize(4 + 2 * (r[1] + r[2]));
		recv_buffer_ptr[pIter++] = 3;
		recv_buffer_ptr[pIter++] = 3 + r[1];
		recv_buffer_ptr[pIter++] = 4 + r[1] + r[1];
		recv_buffer_ptr[pIter++] = 4 + r[1] + r[1] + r[2];
	}

	// Setup send_buffer and recv_buffer
	for (idx_t i = 0; i < numElementLocal; i++) {

		if (part[i] == myrank) {
			recv_buffer[0][recv_buffer_ptr[0]++] = i + vtxdist[myrank];
			recv_buffer[0][recv_buffer_ptr[1]] = recv_buffer[0][recv_buffer_ptr[1]++]
				+ xadj[i + 1] - xadj[i];

			copy(&adjncy[xadj[i]], &adjncy[xadj[i + 1]], &recv_buffer[0][recv_buffer_ptr[2]]);
			copy(&adjproc[xadj[i]], &adjproc[xadj[i + 1]], &recv_buffer[0][recv_buffer_ptr[3]]);
			recv_buffer_ptr[2] += xadj[i + 1] - xadj[i];
			recv_buffer_ptr[3] += xadj[i + 1] - xadj[i];
		}
		else {
			for (idx_t s = 0; s < send_buffer.size(); s++) {
				if (part[i] == send_buffer[s][0]) {
					send_buffer[s][send_buffer_ptr[s * 4 + 0]++] = i + vtxdist[myrank];
					send_buffer[s][send_buffer_ptr[s * 4 + 1]] = send_buffer[s][send_buffer_ptr[s * 4 + 1]++]
						+ xadj[i + 1] - xadj[i];

					copy(&adjncy[xadj[i]], &adjncy[xadj[i + 1]], &send_buffer[s][send_buffer_ptr[s * 4 + 2]]);
					copy(&adjproc[xadj[i]], &adjproc[xadj[i + 1]], &send_buffer[s][send_buffer_ptr[s * 4 + 3]]);
					send_buffer_ptr[s * 4 + 2] += xadj[i + 1] - xadj[i];
					send_buffer_ptr[s * 4 + 3] += xadj[i + 1] - xadj[i];
					break;
				}
			}
		}
	}

	// Free Memory
	vector<idx_t>().swap(send_buffer_ptr);
	vector<idx_t>().swap(recv_buffer_ptr);
	delete[] xadj;
	delete[] vtxdist;
	delete[] adjproc;
	delete[] part;


	vector<idx_t> send_size(numproc, 0);
	vector<idx_t> recv_size(numproc, 0);
	for (const auto& s : send_buffer) {
		send_size[s[0]] = s.size();
	}

	MPI_Alltoall(&send_size[0], 1, IDX_T, &recv_size[0], 1, IDX_T, myComm);

	// Allocate for recv_buffer
	for (OCP_INT i = 0; i < numproc; i++) {
		if (recv_size[i] > 0) {
			recv_buffer.push_back(vector<idx_t>(recv_size[i], i));
		}
	}

	// Free Memory
	vector<idx_t>().swap(send_size);
	vector<idx_t>().swap(recv_size);


	// Communicate for buffer
	for (auto& s : send_buffer) {
		// cout << "Second stage:  " << myrank << " send " << s.size() - 1 << "s to " << s[0] << endl;
		MPI_Isend(s.data() + 1, s.size() - 1, IDX_T, s[0], 0, myComm, &request);
	}
	for (auto& r : recv_buffer) {
		if (r[0] == myrank)  continue;
		// cout << "Second stage:  " << myrank << " recv " << r.size() - 1 << "s from " << r[0] << endl;
		MPI_Recv(r.data() + 1, r.size() - 1, IDX_T, r[0], 0, myComm, &status);
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