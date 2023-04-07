/*! \file    Domain.hpp
 *  \brief   Domain class declaration
 *  \author  Shizhe Li
 *  \date    Feb/28/2023
 *
 *  \note    The params used in OpenCAEPoro is mostly compatible with Eclipse by SLB,
 *           but it has some own rules for easy to use. It is extensible and friendly.
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __DOMAIN_HEADER__
#define __DOMAIN_HEADER__


#include "OCPConst.hpp"
#include "Partition.hpp"

#include <parmetis.h>
#include <vector>
#include <map>

using namespace std;

/// Domain is the new distribution after the partition
class Domain
{
	friend class Reservoir;
	friend class Bulk;
	friend class Well;
	friend class LinearSystem;

	// tmp 
	friend class IsoT_IMPEC;
	friend class IsoT_FIM;
	friend class IsoT_FIMn;
	friend class IsoT_AIMc;
	friend class T_FIM;

public:
	void Setup(const Partition& part);
	vector<OCP_USI> GetWell() const { return well; }
	OCP_USI GetNumGridInterior()const { return numGridInterior; }

	MPI_Comm      myComm;
	OCP_INT       numproc, myrank;

protected:

	OCP_USI numElementTotal; ///< Num of Total Elements
	OCP_USI numWellTotal;    ///< Num of Total Wells
	OCP_USI numElementLocal; ///< Num of local Elements: numGridInterior + numWellLocal

	OCP_USI numGridInterior; ///< Num of interior grid stored in current process
	OCP_USI numGridGhost;    ///< Num of ghost grid stored in current process
	OCP_USI numGridLocal;    ///< numGridInterior + numGridGhost
	USI     numWellLocal;    ///< Num of well stored in current process, static


	vector<vector<idx_t>>   elementCSR;  ///< elements received from other process in CSR
	vector<OCP_USI>         grid;        ///< interior grid & ghost grid
	vector<OCP_USI>         well;        ///< well's index(index starts from zero)
	vector<vector<OCP_USI>> well2Bulk;   ///< wells' bulk neighbor, init global index
	vector<USI>             neighborNum; ///< num of neighbor of grids: well-included, self-included

	map<OCP_USI, OCP_USI>   init_global_to_local; ///< interior grid & ghost grid


	////////////////////////////////////////
	// Communication
	////////////////////////////////////////

public:
	const vector<OCP_USI>* CalGlobalIndex(const USI& nw) const;

protected:

	mutable USI             numActWell;  ///< num of active well

	////////////////////////////////////////
	// Tacit Communication (Prefered)
	////////////////////////////////////////

	vector<vector<OCP_USI>> send_element_loc;
	// vector<vector<OCP_CHAR>> send_buffer;
	vector<vector<OCP_USI>> recv_element_loc;
	// vector<vector<OCP_CHAR>> recv_buffer;

	////////////////////////////////////////
	// Global Index Communication
	////////////////////////////////////////

	mutable vector<USI>     global_index;  ///< Interior grid + active well + ghost grid


	// If all grids are active, if not, record their global index(include inactive grid)
	OCP_BOOL          allActive{ OCP_TRUE };
	vector<OCP_USI>   gridAllIndex;


	// Grid Info
public:
	void SetGirdDimens(const OCP_USI& NX, const OCP_USI& NY, const OCP_USI& NZ) {
		nx = NX; 
		ny = NY;
		nz = NZ;
	}

protected:
	OCP_USI nx, ny, nz;      ///< Dimensions of Grid (Orthogonal only)

	// well perf
public:
	OCP_INT GetPerfLocation(const OCP_USI& wId, const OCP_USI& tmpI, const OCP_USI& tmpJ, const OCP_USI& tmpK) const;
};




#endif



 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Feb/28/2023      Create file                          */
 /*----------------------------------------------------------------------------*/