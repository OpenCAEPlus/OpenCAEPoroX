/*! \file    PreProcess.hpp
 *  \brief   PreProcess for OpenCAEPoroX simulator
 *  \author  Shizhe Li
 *  \date    Feb/15/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PREPROCESS_HEADER__
#define __PREPROCESS_HEADER__

#include "OCPConst.hpp"
#include "ParamRead.hpp"
#include "PreParamGridWell.hpp"
#include "Partition.hpp"
#include "Domain.hpp"
#include "UtilTiming.hpp"
#include "OCPTimeRecord.hpp"

using namespace std;

/// Input grid param -> grid partion -> domain
class PreProcess
{
	friend class OpenCAEPoroX;
	friend class Reservoir;
	

public:

	PreProcess(const string& myFile, const OCP_INT& myRank, MPI_Comm comm);

protected:

	void GetFile(const string& myFile);

protected:

	string inputFile;    ///< Input file with its path (absolute or relative).
	string workdir;      ///< Current work directory.
	string filename;     ///< File name of input file.

	PreParamGridWell preParamGridWell; ///< Param of grids and wells

	Partition        partition;        ///< Partition with Parmetis 

	Domain           domain;   

};



#endif /* end if __PREPROCESS_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Feb/15/2023      Create file                          */
/*----------------------------------------------------------------------------*/