/*! \file    MixtureUnit.cpp
 *  \brief   MixtureUnit class declaration
 *  \author  Shizhe Li
 *  \date    Oct/04/2023 
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


#include "MixtureUnit.hpp"

MixtureUnit::MixtureUnit(const ParamReservoir& rs_param, const USI& i, BulkOptionalModules& opts)
{
	if (rs_param.comps) {
		mix = new OCPMixtureComp(rs_param, i, opts);
	}
	else {
		mix = new OCPMixtureK(rs_param, i, opts);
	}

	vs = &mix->GetVarSet();
	/// Optional Features
	// Calculate surface tension
	surTen = &opts.surTen;
	stMethodIndex = surTen->Setup(rs_param, i, opts.nb);
	// Miscible Factor
	misFac = &opts.misFac;
	mfMethodIndex = misFac->Setup(rs_param, i, opts.nb, surTen);
}


void MixtureUnit::Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) const
{
	mix->Flash(Pin, Tin, Niin);
}


void MixtureUnit::InitFlashIMPEC(const OCP_USI& bId, const BulkVarSet& bvs) const
{
	mix->InitFlash(bId, bvs);
	surTen->CalSurfaceTension(bId, stMethodIndex, *vs);
	misFac->CalMiscibleFactor(bId, mfMethodIndex);
}


void MixtureUnit::InitFlashFIM(const OCP_USI& bId, const BulkVarSet& bvs) const
{
	mix->InitFlashDer(bId, bvs);
	surTen->CalSurfaceTension(bId, stMethodIndex, *vs);
	misFac->CalMiscibleFactor(bId, mfMethodIndex);
}


void MixtureUnit::FlashIMPEC(const OCP_USI& bId, const BulkVarSet& bvs) const
{
	mix->Flash(bId, bvs);
	surTen->CalSurfaceTension(bId, stMethodIndex, *vs);
	misFac->CalMiscibleFactor(bId, mfMethodIndex);
}


void MixtureUnit::FlashFIM(const OCP_USI& bId, const BulkVarSet& bvs) const
{
	mix->FlashDer(bId, bvs);
	surTen->CalSurfaceTension(bId, stMethodIndex, *vs);
	misFac->CalMiscibleFactor(bId, mfMethodIndex);
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/04/2023      Create file                          */
/*----------------------------------------------------------------------------*/