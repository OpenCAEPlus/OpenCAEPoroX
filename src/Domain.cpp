/*! \file    Partition.cpp
 *  \brief   Partition for OpenCAEPoro simulator
 *  \author  Shizhe Li
 *  \date    Feb/28/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "Domain.hpp"



void Domain::Setup(const Partition& part)
{
	myComm  = part.myComm;	
	numproc = part.numproc;
	myrank  = part.myrank;

	if (numproc == 1) {
		numElementTotal = part.numElementTotal;
		numWellTotal    = part.numWellTotal;
		numGridInterior = numElementTotal - numWellTotal;
		numGridGhost    = 0;
		grid.resize(numGridInterior);
		for (OCP_USI n = 0; n < numGridInterior; n++) grid[n] = n;
		well.resize(numWellTotal);
		for (OCP_USI w = 0; w < numWellTotal; w++)    well[w] = w;
		return;
	}

	// Setup Domain
	elementCSR.swap(part.elementCSR);	
	numElementTotal = part.numElementTotal;
	numWellTotal    = part.numWellTotal;
	numElementLocal = 0;
	set<OCP_USI> recv_proc;
	for (const auto& e : elementCSR) {
		recv_proc.insert(e[0]);
		numElementLocal += e[1];
	}

	// Traverse elementCSR in ascending order of process number
	idx_t    global_well_start = numElementTotal - numWellTotal;
	USI      localIndex        = 0;
	vector<set<OCP_USI>> ghostElement;
	grid.reserve(numElementLocal * 1.5); // preserved space
	for (const auto& s : recv_proc) {
		for (const auto& e : elementCSR) {
			if (e[0] == s) {
				const idx_t* my_vtx  = &e[3];
				const idx_t* my_xadj = &e[3 + e[1]];
				const idx_t* my_edge = &e[3 + e[1] + e[1] + 1];
				const idx_t* my_edge_proc = &e[3 + e[1] + e[1] + 1 + e[2]];
				for (USI i = 0; i < e[1]; i++) {
					if (my_vtx[i] >= global_well_start) {
						// well
						well.push_back(my_vtx[i] - (numElementTotal - numWellTotal));
						well2Bulk.push_back(vector<OCP_USI>(my_edge + my_xadj[i], my_edge + my_xadj[i + 1]));
						continue;
					}
					grid.push_back(my_vtx[i]);
					init_global_to_local.insert(make_pair(static_cast<OCP_USI>(my_vtx[i]), localIndex));
					for (USI j = my_xadj[i]; j < my_xadj[i + 1]; j++) {
						if (my_edge_proc[j] != myrank) {
							// current interior grid is also ghost grid of other process
							USI s = 0;
							for (s = 0; s < send_element_loc.size(); s++) {
								if (my_edge_proc[j] == send_element_loc[s][0]) {
									if (localIndex != send_element_loc[s].back()) {
										send_element_loc[s].push_back(localIndex);
									}
									ghostElement[s].insert(static_cast<OCP_USI>(my_edge[j]));
									break;
								}
							}
							if (s == send_element_loc.size()) {
								send_element_loc.push_back(vector<OCP_USI>{static_cast<OCP_USI>(my_edge_proc[j]), localIndex});
								ghostElement.push_back(set<OCP_USI>{static_cast<OCP_USI>(my_edge[j])});
							}
						}
					}
					localIndex++;
				}
				break;
			}
		}
	}

	numGridInterior = init_global_to_local.size();
	numWellLocal    = well.size();
	OCP_ASSERT(numGridInterior == numElementLocal - numWellLocal, "");

	recv_element_loc.resize(send_element_loc.size());
	for (USI i = 0; i < send_element_loc.size(); i++) {
		recv_element_loc[i].reserve(3);
		recv_element_loc[i].push_back(send_element_loc[i][0]);
	}

	numGridGhost = 0;
	localIndex = init_global_to_local.size();
	USI sIter = 0;
	for (const auto& g : ghostElement) {
		numGridGhost += g.size();
		recv_element_loc[sIter].push_back(localIndex);
		for (const auto& g1 : g) {
			grid.push_back(g1);
			init_global_to_local.insert(make_pair(g1, localIndex));
			localIndex++;
		}
		recv_element_loc[sIter].push_back(localIndex);
		sIter++;
	}
	numGridLocal = numGridInterior + numGridGhost;


	MPI_Barrier(myComm);
	cout << myrank << " begin output" << endl;

	//////////////////////////////////////////////////////////////
	if (true) {
		ofstream myFile;
		myFile.open("test/process" + to_string(myrank) + ".txt");
		ios::sync_with_stdio(false);
		myFile.tie(0);

		for (const auto& e : elementCSR) {
			// num grid, num edges
			myFile << setw(8) << e[0] << setw(8) << e[1] << setw(8) << e[2] << endl;
			const idx_t* my_vtx = &e[3];
			const idx_t* my_xadj = &e[3 + e[1]];
			const idx_t* my_edge = &e[3 + e[1] + e[1] + 1];
			const idx_t* my_edge_proc = &e[3 + e[1] + e[1] + 1 + e[2]];

			for (int i = 0; i < e[1]; i++) {
				myFile << setw(8) << my_vtx[i];
			}
			myFile << endl;
			// vertex
			for (int i = 0; i < e[1]; i++) {
				for (int j = my_xadj[i]; j < my_xadj[i + 1]; j++) {
					myFile << setw(8) << my_vtx[i];
				}
			}
			myFile << endl;
			// edges
			for (int i = 0; i < e[2]; i++) {
				myFile << setw(8) << my_edge[i];
			}
			myFile << endl;
			// process
			for (int i = 0; i < e[2]; i++) {
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
			for (const auto& s1 : s) {
				myFile << setw(8) << s1;
			}
			myFile << endl;
		}
		myFile << endl << endl << endl << "recv" << endl;
		for (const auto& r : recv_element_loc) {
			for (const auto& r1 : r) {
				myFile << setw(8) << r1;
			}
			myFile << endl;
		}
		myFile << endl << endl << endl << "well" << endl;
		for (USI w = 0; w < well.size(); w++) {
			myFile << setw(8) << well[w];
			for (const auto& p : well2Bulk[w]) {
				myFile << setw(8) << p;
			}
		}
		myFile << endl << endl << endl;
		myFile << "Grid Num" << endl;
		myFile << setw(8) << numElementLocal << setw(8) << numGridInterior
			<< setw(8) << numWellLocal << setw(8) << numGridGhost << endl;
		myFile.close();
	}
	//////////////////////////////////////////////////////////////
}


OCP_INT Domain::GetPerfLocation(const OCP_USI& wId, const OCP_USI& tmpI, const OCP_USI& tmpJ, const OCP_USI& tmpK) const
{
	auto myFind = [](const OCP_USI pId, const OCP_USI* grid, OCP_USI lp, OCP_USI rp) {
		OCP_USI cp = (lp + rp) / 2;
		OCP_USI lcp = 0;
		while (lcp != cp) {
			lcp = cp;
			if (pId < grid[cp])   rp = cp;
			else                  lp = cp;
			cp = (lp + rp) / 2;
		}
		return cp;
	};

	// find if pid is in gridAllIndex
	const OCP_USI*    grid_range;
	if (allActive)    grid_range = grid.data();
	else              grid_range = gridAllIndex.data();


	OCP_USI pId = tmpK * nx * ny + tmpJ * nx + tmpI;
	const OCP_USI loc = myFind(pId, grid_range, 0, numGridInterior);

	if (pId != grid_range[loc]) {
		// this perf is in inactive bulk
		return -1;
	}
	else {
		// determine if this perf is active
		pId = grid[loc]; // active idex
		OCP_USI tmploc = myFind(pId, well2Bulk[wId].data(), 0, well2Bulk[wId].size());
		if (pId != well2Bulk[wId][tmploc]) {
			return -1;
		}
	}
	return loc;
}


const vector<OCP_USI>* Domain::CalGlobalIndex(const USI& nw) const
{
	numActWell = nw;
	global_index.resize(numGridLocal + numActWell);

	const OCP_INT numElementLoc = numGridInterior + numActWell;
	OCP_INT       global_begin;
	OCP_INT       global_end;
	MPI_Scan(&numElementLoc, &global_end, 1, MPI_INT, MPI_SUM, myComm);

	global_begin = global_end - numElementLoc;
	global_end   = global_end - 1;

	// Get Interior grid's global index
	for (OCP_USI n = 0; n < numElementLoc; n++)
		global_index[n] = n + global_begin;

	// Get Ghost grid's global index by communication
	vector<vector<OCP_INT>> recv_buffer(recv_element_loc.size());
	vector<vector<OCP_INT>> send_buffer(send_element_loc.size());
	for (USI i = 0; i < recv_buffer.size(); i++) {
		const vector<OCP_USI>& r = recv_element_loc[i];
		recv_buffer[i].resize(r[2] - r[1] + 1);
		recv_buffer[i][0] = r[0];
	}
	for (USI i = 0; i < send_buffer.size(); i++) {
		const vector<OCP_USI>& s = send_element_loc[i];
		send_buffer[i].resize(s.size());
		send_buffer[i][0] = s[0];
		for (USI j = 1; j < s.size(); j++) {
			send_buffer[i][j] = global_index[s[j]];
		}
	}

	MPI_Request request;
	MPI_Status  status;
	for (auto& s : send_buffer) {
		MPI_Isend(s.data() + 1, s.size() - 1, MPI_INT, s[0], 0, myComm, &request);
	}
	for (auto& r : recv_buffer) {
		MPI_Recv(r.data() + 1, r.size() - 1, MPI_INT, r[0], 0, myComm, &status);
	}

	// Assemble to global_index
	for (USI i = 0; i < recv_buffer.size(); i++) {
		const OCP_USI bId = recv_element_loc[i][1] + numActWell;
		copy(&recv_buffer[i][1], &recv_buffer[i].back() + 1, &global_index[bId]);
	}

	MPI_Barrier(myComm);
	return &global_index;
}




 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Feb/28/2023      Create file                          */
 /*----------------------------------------------------------------------------*/