/*! \file    Domain.hpp
 *  \brief   Domain class declaration
 *  \author  Shizhe Li
 *  \date    Feb/28/2023
 *
 *  \note    The params used in OpenCAEPoroX is mostly compatible with Eclipse by SLB,
 *           but it has some own rules for easy to use. It is extensible and friendly.
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __DOMAIN_HEADER__
#define __DOMAIN_HEADER__


#include "OCPConst.hpp"
#include "Partition.hpp"
#include "UtilTiming.hpp"
#include "OCPTimeRecord.hpp"
#include "OCPGroupProcess.hpp"

#include <parmetis.h>
#include <vector>
#include <map>
#include <unordered_map>

using namespace std;

/// Domain is the new distribution after the partition
// all elements(well) are active, all connections(perforations) are active
// dead elements are excluded, so global index is also for active elements
class Domain
{
	friend class Reservoir;
	friend class Bulk;
	friend class Well;
	friend class PeacemanWell;
	friend class LinearSystem;
	friend class OCPMatrix;

	// Method(tmp)
	friend class IsoT_IMPEC;
	friend class IsoT_FIM;
	friend class IsoT_FIMddm;
	friend class IsoT_AIMc;
	friend class T_FIM;

	// Linear Solver (tmp)
	friend class SamgSolver;

public:
	void Setup(const Partition& part, const PreParamGridWell& gridwell);
	auto GetNumGridTotal() const { return numElementTotal - numWellTotal; }
	auto GetNumGridInterior() const { return numGridInterior; }
	auto GetWell() const { return well; }
	const auto& GetGrid() const { return grid; }

protected:
	void InitComm(const Partition& part);

	// for global communication
public:
	MPI_Comm      global_comm;
	OCP_INT       global_numproc;
	OCP_INT       global_rank;
	set<OCP_INT>  global_group_rank;

public:
	// reset linear solver communication
	void InitCSComm();
	// for linear solver communication
	void SetCSComm();
	void SetCSComm(const unordered_map<OCP_USI, OCP_DBL>& bk_info);
	OCP_BOOL IfIRankInLSCommGroup(const OCP_INT& p) const;

protected:
	void SetCS01(const unordered_map<OCP_USI, OCP_DBL>& bk_info, unordered_map<OCP_INT, OCP_INT>& proc_wght);
	void SetCS02(const unordered_map<OCP_USI, OCP_DBL>& bk_info, unordered_map<OCP_INT, OCP_INT>& proc_wght);
	void ProcWeight_f2i(const unordered_map<OCP_INT, OCP_DBL>& tmp_proc_wght, unordered_map<OCP_INT, OCP_INT>& proc_wght);

public:
	MPI_Comm          cs_comm;
	OCP_INT           cs_numproc;
	OCP_INT           cs_rank;
	set<OCP_INT>      cs_group_global_rank;
	set<OCP_INT>      cs_group_local_rank;

protected:
	/// Num of Total Elements(grids + wells)
	OCP_ULL numElementTotal; 
	/// Num of Total Wells
	OCP_USI numWellTotal; 
	/// Num of local Elements: numGridInterior + numWellLocal
	OCP_USI numElementLocal; 
	/// Num of interior grid stored in current process
	OCP_USI numGridInterior; 
	/// Num of ghost grid stored in current process
	OCP_USI numGridGhost;   
	/// numGridInterior + numGridGhost
	OCP_USI numGridLocal; 
	/// Num of well stored in current process
	USI     numWellLocal;    

	/// elements received from other process in CSR
	map<OCP_INT, vector<idx_t>> elementCSR;
	/// initial global index of interior grid & ghost grid
	vector<OCP_ULL>             grid;
	/// well's global index(index starts from zero)
	vector<OCP_USI>             well; 

	/// Well's global index(start from zero), perfortions' index(start from zero) and location(bulks' local index) of perforation
	vector<vector<OCP_USI>> wellWPB;
	/// num of all neighbor of grids: well-included, self-included
	vector<USI>             neighborNum; 
	/// initial global index -> local index (interior grid & ghost grid)
	unordered_map<OCP_ULL, OCP_USI>   init_global_to_local;


	////////////////////////////////////////
	// Communication
	////////////////////////////////////////

public:
	const vector<OCP_ULL>* CalGlobalIndex() const;

	////////////////////////////////////////
	// Tacit Communication (Prefered, Local Index)
	////////////////////////////////////////

	map<OCP_INT, set<OCP_USI>>    send_element_loc;
	map<OCP_INT, vector<OCP_USI>> recv_element_loc;

	mutable vector<MPI_Request>  send_request;
	mutable vector<MPI_Request>  recv_request;

	////////////////////////////////////////
	// Global Index Communication
	////////////////////////////////////////

	mutable vector<OCP_ULL>     global_index;  ///< Interior grid + active well + ghost grid in equations

public:
	OCP_USI GetNumActElementForSolver() const { return numGridInterior + numActWellLocal; }
	void SetNumActWellLocal(const OCP_DBL& nw) const { numActWellLocal = nw; }

protected:
	/// number of active well
	mutable USI numActWellLocal;

	// Well perforations
public:
	OCP_INT GetPerfLocation(const OCP_USI& wId, const USI& p) const;
	/// Return number of perforations of specific well(given well global index)
	USI     GetPerfNum(const OCP_USI& wId) const;
};


#endif



 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Feb/28/2023      Create file                          */
 /*----------------------------------------------------------------------------*/