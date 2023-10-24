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

#include <parmetis.h>
#include <vector>
#include <map>

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

	// Method(tmp)
	friend class IsoT_IMPEC;
	friend class IsoT_FIM;
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

public:
	MPI_Comm                     myComm;
	OCP_INT                      numproc, myrank;

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
	vector<vector<idx_t>>   elementCSR;  
	///< global index of interior grid & ghost grid
	vector<OCP_ULL>         grid;
	/// well's global index(index starts from zero)
	vector<OCP_USI>         well; 

	/// Well's global index(start from zero), perfortions' index(start from zero) and location(bulks' local index) of perforation
	vector<vector<OCP_USI>> wellWPB;
	/// num of all neighbor of grids: well-included, self-included
	vector<USI>             neighborNum; 
	/// global index -> local index (interior grid & ghost grid)
	map<OCP_ULL, OCP_USI>   init_global_to_local;

	////////////////////////////////////////
	// Communication
	////////////////////////////////////////

public:
	const vector<OCP_ULL>* CalGlobalIndex(const USI& nw) const;

	////////////////////////////////////////
	// Tacit Communication (Prefered, Local Index)
	////////////////////////////////////////

	USI numSendProc, numRecvProc;

	vector<vector<OCP_USI>> send_element_loc;
	// vector<vector<OCP_CHAR>> send_buffer;
	vector<vector<OCP_USI>> recv_element_loc;
	// vector<vector<OCP_CHAR>> recv_buffer;

	mutable vector<MPI_Request>  send_request;
	mutable vector<MPI_Request>  recv_request;

	////////////////////////////////////////
	// Global Index Communication
	////////////////////////////////////////

	mutable vector<OCP_ULL>     global_index;  ///< Interior grid + active well + ghost grid in equations

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