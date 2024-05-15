/*! \file    Partition.cpp
 *  \brief   Partition for OpenCAEPoroX simulator
 *  \author  Shizhe Li
 *  \date    Feb/28/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "Domain.hpp"



void Domain::Setup(const Partition& part, const PreParamGridWell& gridwell)
{
	InitComm(part);

	if (global_numproc == 1) {
		numElementTotal = part.numElementTotal;
		numElementLocal = numElementTotal;
		numWellTotal    = part.numWellTotal;
		numGridInterior = numElementTotal - numWellTotal;
		numGridGhost    = 0;
		numGridLocal    = numGridInterior;

		grid.resize(numGridInterior);
		for (OCP_USI n = 0; n < numGridInterior; n++) {
			grid[n] = n;
			init_global_to_local.insert(make_pair(n, n));
		}
		well.resize(numWellTotal);
		for (OCP_USI w = 0; w < numWellTotal; w++) {
			well[w] = w;
		}
		return;
	}

	if (CURRENT_RANK == MASTER_PROCESS) {
		OCP_INFO("Set Domain -- begin");
	}


	// Setup Domain
	elementCSR.swap(part.elementCSR);	
	numElementTotal = part.numElementTotal;
	numWellTotal    = part.numWellTotal;
	numElementLocal = 0;

	for (const auto& e : elementCSR) {
		numElementLocal += e.second[0];
	}

	// Traverse elementCSR in ascending order of process number
	OCP_ULL  global_well_start = numElementTotal - numWellTotal;
	OCP_USI  localIndex        = 0;
	map<OCP_INT, set<OCP_ULL>> ghostElement;
	grid.reserve(numElementLocal * 1.5); // preserved space
	for (const auto& e : elementCSR) {
		const auto&  ev           = e.second;
		const idx_t* my_vtx       = &ev[2];
		const idx_t* my_xadj      = &ev[2 + ev[0]];
		const idx_t* my_edge      = &ev[2 + ev[0] + (ev[0] + 1)];
		const idx_t* my_edge_proc = &ev[2 + ev[0] + (ev[0] + 1) + ev[1]];
		for (USI i = 0; i < ev[0]; i++) {
			if (my_vtx[i] >= global_well_start) {
				// well
				well.push_back(my_vtx[i] - global_well_start);
				continue;
			}
			grid.push_back(my_vtx[i]);
			init_global_to_local.insert(make_pair(static_cast<OCP_ULL>(my_vtx[i]), localIndex));
			for (USI j = my_xadj[i]; j < my_xadj[i + 1]; j++) {
				const OCP_INT proc = my_edge_proc[j];
				if (proc != global_rank) {
					// current interior grid is also ghost grid of other process
					send_element_loc[proc].insert(localIndex);
					ghostElement[proc].insert(my_edge[j]);
				}
			}
			localIndex++;
		}
	}

	numGridInterior = init_global_to_local.size();
	numWellLocal    = well.size();
	OCP_ASSERT(numGridInterior == numElementLocal - numWellLocal, "");

	numGridGhost = 0;
	localIndex   = init_global_to_local.size();

	for (const auto& g : ghostElement) {
		numGridGhost += g.second.size();
		recv_element_loc[g.first].push_back(localIndex);
		for (const auto& g1 : g.second) {
			grid.push_back(g1);
			init_global_to_local.insert(make_pair(g1, localIndex));
			localIndex++;
		}
		recv_element_loc[g.first].push_back(localIndex);
	}

	numGridLocal = numGridInterior + numGridGhost;
	send_request.resize(send_element_loc.size());
	recv_request.resize(recv_element_loc.size());
	grid.shrink_to_fit();

	//////////////////////////////////////////////////////////////
	// Output partition information
	//////////////////////////////////////////////////////////////

	if (false) {
		ofstream myFile;
		myFile.open("test/process" + to_string(global_rank) + ".txt");
		ios::sync_with_stdio(false);
		myFile.tie(0);

		for (const auto& e : elementCSR) {
			const auto& ev = e.second;
			// process, num grid, num edges
			myFile << setw(8) << e.first << setw(8) << ev[0] << setw(8) << ev[1] << endl;
			const idx_t* my_vtx       = &ev[2];
			const idx_t* my_xadj      = &ev[2 + ev[0]];
			const idx_t* my_edge      = &ev[2 + ev[0] + ev[0] + 1];
			const idx_t* my_edge_proc = &ev[2 + ev[0] + ev[0] + 1 + ev[1]];
			// vertex
			for (int i = 0; i < ev[0]; i++) {
				myFile << setw(8) << my_vtx[i];
			}
			myFile << endl;
			// vertex
			for (int i = 0; i < ev[0]; i++) {
				for (int j = my_xadj[i]; j < my_xadj[i + 1]; j++) {
					myFile << setw(8) << my_vtx[i];
				}
			}
			myFile << endl;
			// edges
			for (int i = 0; i < ev[1]; i++) {
				myFile << setw(8) << my_edge[i];
			}
			myFile << endl;
			// process
			for (int i = 0; i < ev[1]; i++) {
				myFile << setw(8) << my_edge_proc[i];
			}
			myFile << endl << endl << endl;
		}

		myFile << "init_global_to_local" << endl;
		for (const auto& m : init_global_to_local) {
			myFile << setw(8) << m.first;
		}
		myFile << endl;
		for (const auto& m : init_global_to_local) {
			myFile << setw(8) << m.second;
		}
		myFile << endl << endl;
		myFile << "local_to_init_global" << endl;
		for (int i = 0; i < grid.size(); i++) {
			myFile << setw(8) << i;
		}
		myFile << endl;
		for (const auto& e : grid) {
			myFile << setw(8) << e;
		}
		myFile << endl << endl << endl << "send" << endl;
		for (const auto& s : send_element_loc) {
			myFile << setw(8) << s.first;
			for (const auto& s1 : s.second) {
				myFile << setw(8) << s1;
			}
			myFile << endl;
		}
		myFile << endl << endl << endl << "recv" << endl;
		for (const auto& r : recv_element_loc) {
			myFile << setw(8) << r.first;
			for (const auto& r1 : r.second) {
				myFile << setw(8) << r1;
			}
			myFile << endl;
		}
		myFile << endl << endl << endl << "well" << endl;
		for (USI w = 0; w < well.size(); w++) {
			myFile << setw(8) << well[w];
		}
		myFile << endl << endl << endl;
		myFile << "Grid Num" << endl;
		myFile << setw(8) << numElementLocal << setw(8) << numGridInterior
			<< setw(8) << numWellLocal << setw(8) << numGridGhost << endl;
		myFile.close();
	}
	//////////////////////////////////////////////////////////////


	if (CURRENT_RANK == MASTER_PROCESS) {
		OCP_INFO("Set Domain -- end");
	}
}


void Domain::InitComm(const Partition& part)
{
	global_comm    = part.myComm;
	global_numproc = part.numproc;
	global_rank    = part.myrank;
	global_group_rank.clear();
	for (OCP_USI n = 0; n < global_numproc; n++) {
		global_group_rank.insert(n);
	}

	InitCSComm();
}


void Domain::InitCSComm()
{
	if (cs_comm != MPI_COMM_NULL)  MPI_Comm_free(&cs_comm);
	MPI_Comm_dup(global_comm, &cs_comm);

	cs_numproc           = global_numproc;
	cs_rank              = global_rank;
	cs_group_global_rank = global_group_rank;
	cs_group_local_rank  = cs_group_global_rank;
}


void Domain::SetCSComm(const unordered_map<OCP_USI, OCP_DBL>& bk_info)
{

	GetWallTime timer;
	timer.Start();

	unordered_map<OCP_INT, OCP_INT> proc_weight;

	SetCS01(bk_info, proc_weight);
	// SetCS02(bk_info, proc_weight);

	GroupProcess(GroupMethod::ParMetis_gather, cs_group_global_rank, proc_weight, cs_comm, global_comm);

	MPI_Comm_size(cs_comm, &cs_numproc);
	MPI_Comm_rank(cs_comm, &cs_rank);

	cs_group_local_rank.clear();
	for (OCP_USI n = 0; n < cs_numproc; n++) {
		cs_group_local_rank.insert(n);
	}

	OCPTIME_GROUPPROCESS += timer.Stop();
}


void Domain::SetCS01(const unordered_map<OCP_USI, OCP_DBL>& bk_info, unordered_map<OCP_INT, OCP_INT>& proc_wght)
{
	unordered_map<OCP_INT, OCP_DBL> tmp_proc_wght;

	cs_group_global_rank.clear();
	for (const auto& b : bk_info) {
		if (b.first < numGridInterior) {
			for (const auto& s : send_element_loc) {
				const auto& sv = s.second;
				if (sv.count(b.first)) {
					cs_group_global_rank.insert(s.first);
					if (tmp_proc_wght.count(s.first))   tmp_proc_wght[s.first] += b.second;
					else                                tmp_proc_wght[s.first] = b.second;
				}
			}
		}
		else {
			for (const auto& r : recv_element_loc) {
				const auto& rv = r.second;
				if (b.first >= rv[0] && b.first < rv[1]) {
					cs_group_global_rank.insert(r.first);
					if (tmp_proc_wght.count(r.first))   tmp_proc_wght[r.first] += b.second;
					else                                tmp_proc_wght[r.first] = b.second;
					break;
				}
			}
		}
	}

	ProcWeight_f2i(tmp_proc_wght, proc_wght);
}


void Domain::SetCS02(const unordered_map<OCP_USI, OCP_DBL>& bk_info, unordered_map<OCP_INT, OCP_INT>& proc_wght)
{
	unordered_map<OCP_INT, OCP_DBL> tmp_proc_wght;

	OCP_DBL selfW = 0;
	cs_group_global_rank.clear();
	for (const auto& b : bk_info) {
		if (b.first < numGridInterior) {
			for (const auto& s : send_element_loc) {
				const auto& sv = s.second;
				if (sv.count(b.first)) {
					cs_group_global_rank.insert(s.first);
					if (tmp_proc_wght.count(s.first))   tmp_proc_wght[s.first] += b.second;
					else                                tmp_proc_wght[s.first] = b.second;
				}
			}
			selfW += b.second;
		}
		else {
			for (const auto& r : recv_element_loc) {
				const auto& rv = r.second;
				if (b.first >= rv[0] && b.first < rv[1]) {
					cs_group_global_rank.insert(r.first);
					if (tmp_proc_wght.count(r.first))   tmp_proc_wght[r.first] += b.second;
					else                                tmp_proc_wght[r.first] = b.second;
					break;
				}
			}
		}
	}

	// add selfW to tmp_proc_wght
	for (const auto& r : recv_element_loc) {
		cs_group_global_rank.insert(r.first);

		if (tmp_proc_wght.count(r.first))  tmp_proc_wght[r.first] += selfW;
		else                               tmp_proc_wght[r.first] = selfW;
	}
	// receive other process' selfW
	vector<OCP_DBL> otherW(recv_element_loc.size());

	GetWallTime timer;
	timer.Start();
	// Get Ghost grid's global index by communication
	USI iter = 0;
	for (const auto& r : recv_element_loc) {

		MPI_Irecv(&otherW[iter], 1, OCPMPI_DBL, r.first, 0, global_comm, &recv_request[iter]);
		iter++;
	}

	iter = 0;
	for (const auto& s : send_element_loc) {

		MPI_Isend(&selfW, 1, OCPMPI_DBL, s.first, 0, global_comm, &send_request[iter]);
		iter++;
	}


	MPI_Waitall(iter, recv_request.data(), MPI_STATUS_IGNORE);
	MPI_Waitall(iter, send_request.data(), MPI_STATUS_IGNORE);

	OCPTIME_COMM_P2P += timer.Stop();

	iter = 0;
	for (const auto& r : recv_element_loc) {
		tmp_proc_wght[r.first] += otherW[iter++];
	}

	ProcWeight_f2i(tmp_proc_wght, proc_wght);

	//if (OCP_FALSE)
	//{
	//	MPI_Barrier(global_comm);
	//	std::this_thread::sleep_for(std::chrono::milliseconds(CURRENT_RANK * 200));
	//	cout << "Rank: " << CURRENT_RANK << endl;
	//	for (const auto& p : tmp_proc_wght) {
	//		cout << CURRENT_RANK << " - " << p.first << "  : "
	//			 << fixed << setprecision(12) << p.second << endl;
	//	}
	//}
}


void Domain::ProcWeight_f2i(const unordered_map<OCP_INT, OCP_DBL>& tmp_proc_wght, unordered_map<OCP_INT, OCP_INT>& proc_wght)
{
	OCP_DBL minW = 1E20;
	OCP_DBL maxW = 0.0;

	for (const auto& w : tmp_proc_wght) {
		minW = min(w.second, minW);
		maxW = max(w.second, maxW);
	}

	for (const auto& w : tmp_proc_wght) {
		proc_wght[w.first] = w.second / (minW + 1) + 1;
	}
}


OCP_BOOL Domain::IfIRankInLSCommGroup(const OCP_INT& p) const
{
	if (cs_numproc == global_numproc) {
		return OCP_TRUE;
	}
	else if (cs_group_global_rank.count(p)) {
		return OCP_TRUE;
	}
	else {
		return OCP_FALSE;
	}
}


OCP_INT Domain::GetPerfLocation(const OCP_USI& wId, const USI& p) const
{
	for (const auto& w : wellWPB) {
		if (w[0] == wId) {
			for (USI i = 1; i < w.size(); i += 2) {
				if (w[i] == p) {
					return w[i + 1];
				}
			}
		}
	}
	return -1;
}


USI Domain::GetPerfNum(const OCP_USI& wId) const
{
	for (const auto& w : wellWPB) {
		if (w[0] == wId) {
			return (w.size() - 1) / 2;
		}
	}
	OCP_ABORT("WRONG WELL SETUP!");
}


const vector<OCP_ULL>* Domain::CalGlobalIndex() const
{
	global_index.resize(numGridLocal + numActWellLocal);

	const OCP_ULL numElementLoc = numGridInterior + numActWellLocal;
	OCP_ULL       global_begin;
	OCP_ULL       global_end;

	GetWallTime timer;
	timer.Start();

	MPI_Scan(&numElementLoc, &global_end, 1, OCPMPI_ULL, MPI_SUM, cs_comm);

	OCPTIME_COMM_COLLECTIVE += timer.Stop();

	global_begin = global_end - numElementLoc;
	global_end   = global_end - 1;

	// Get Interior grid's global index
	for (OCP_USI n = 0; n < numElementLoc; n++)
		global_index[n] = n + global_begin;

	timer.Start();

	// Get Ghost grid's global index by communication
	USI iter = 0;
	for (const auto& r : recv_element_loc) {

		if (!IfIRankInLSCommGroup(r.first))  continue;

		const auto& rv = r.second;
		const auto  bId = rv[0] + numActWellLocal;
		MPI_Irecv(&global_index[bId], rv[1] - rv[0], OCPMPI_ULL, r.first, 0, global_comm, &recv_request[iter]);
		iter++;
	}

	iter = 0;
	vector<vector<OCP_ULL>> send_buffer(send_element_loc.size());
	for (const auto& s : send_element_loc) {

		if (!IfIRankInLSCommGroup(s.first))  continue;

		const auto& sv = s.second;
		auto&       sb = send_buffer[iter];
		sb.reserve(sv.size());
		for (const auto& sv1 : sv) {
			sb.push_back(global_index[sv1]);
		}
		MPI_Isend(sb.data(), sb.size(), OCPMPI_ULL, s.first, 0, global_comm, &send_request[iter]);
		iter++;
	}


	MPI_Waitall(iter, recv_request.data(), MPI_STATUS_IGNORE);
	MPI_Waitall(iter, send_request.data(), MPI_STATUS_IGNORE);

	OCPTIME_COMM_P2P += timer.Stop();

	return &global_index;
}




 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Feb/28/2023      Create file                          */
 /*----------------------------------------------------------------------------*/