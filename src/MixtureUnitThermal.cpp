/*! \file    MixtureUnitThermal_K.cpp
 *  \brief   MixtureUnitThermal_K class declaration
 *  \author  Shizhe Li
 *  \date    Nov/10/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// OpenCAEPoroX header files
#include "MixtureUnitThermal.hpp"

MixtureUnitThermal_OW::MixtureUnitThermal_OW(const ParamReservoir& param, const USI& tarId, OptionalModules& opts)
{
    OWTM = new OCPMixtureK(param, tarId, opts);
    vs   = &OWTM->GetVarSet();
}

void MixtureUnitThermal_OW::Flash(const OCP_DBL& Pin,
                               const OCP_DBL& Tin,
                               const OCP_DBL* Niin)
{
    OWTM->Flash(Pin, Tin, Niin);
}

void MixtureUnitThermal_OW::InitFlashIMPEC(const OCP_USI& bId, const BulkVarSet& bvs)
{
    OWTM->InitFlash(bId, bvs);
}

void MixtureUnitThermal_OW::InitFlashFIM(const OCP_USI& bId, const BulkVarSet& bvs)
{
    OWTM->InitFlashDer(bId, bvs);
}

void MixtureUnitThermal_OW::FlashIMPEC(const OCP_USI& bId, const BulkVarSet& bvs)
{
    OWTM->Flash(bId, bvs);
}

void MixtureUnitThermal_OW::FlashFIM(const OCP_USI& bId, const BulkVarSet& bvs)
{
    OWTM->FlashDer(bId, bvs);
}


OCP_DBL MixtureUnitThermal_OW::CalInjWellEnthalpy(const OCP_DBL& Tin, const OCP_DBL* Ziin)
{
    return OWTM->CalEnthalpy(Tin, Ziin);
}

OCP_DBL MixtureUnitThermal_OW::XiPhase(const OCP_DBL& Pin,
                                    const OCP_DBL& Tin,
                                    const vector<OCP_DBL>& Ziin,
                                    const PhaseType& pt)
{
    return OWTM->CalXi(Pin, Pin, Tin, &Ziin[0], pt);
}

OCP_DBL
MixtureUnitThermal_OW::RhoPhase(const OCP_DBL& Pin,
                             const OCP_DBL& Pbb,
                             const OCP_DBL& Tin,
                             const vector<OCP_DBL>& Ziin,
                             const PhaseType& pt)
{
    return OWTM->CalRho(Pin, Pbb, Tin, &Ziin[0], pt);
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           NOV/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/
