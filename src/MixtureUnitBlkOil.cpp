/*! \file    MixtureBlkOilUnit.cpp
 *  \brief   MixtureBlkOilUnit class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "MixtureUnitBlkOil.hpp"

 ///////////////////////////////////////////////
 // MixtureUnitBlkOil_OW
 ///////////////////////////////////////////////

MixtureUnitBlkOil_OW::MixtureUnitBlkOil_OW(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts)
{
    OWM = new OCPMixtureBlkOilOW(rs_param, i);

    vs          = &OWM->GetVarSet();
}

void MixtureUnitBlkOil_OW::Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{
    OWM->Flash(Pin, Tin, Niin);
}

void MixtureUnitBlkOil_OW::InitFlashIMPEC(const OCP_USI& bId, const BulkVarSet& bvs)
{
    OWM->InitFlash(bId, bvs);
}

void MixtureUnitBlkOil_OW::InitFlashFIM(const OCP_USI& bId, const BulkVarSet& bvs)
{
    OWM->InitFlashDer(bId, bvs);
}

void MixtureUnitBlkOil_OW::FlashIMPEC(const OCP_USI& bId, const BulkVarSet& bvs)
{
    OWM->Flash(bId, bvs);
}

void MixtureUnitBlkOil_OW::FlashFIM(const OCP_USI& bId, const BulkVarSet& bvs)
{
    OWM->FlashDer(bId, bvs);
}

OCP_DBL MixtureUnitBlkOil_OW::XiPhase(const OCP_DBL& Pin,
    const OCP_DBL& Tin,
    const vector<OCP_DBL>& Ziin,
    const PhaseType& pt)
{
    return OWM->CalXi(Pin, Pin, Tin, &Ziin[0], pt);
}

OCP_DBL MixtureUnitBlkOil_OW::RhoPhase(const OCP_DBL& Pin,
    const OCP_DBL& Pbb,
    const OCP_DBL& Tin,
    const vector<OCP_DBL>& Ziin,
    const PhaseType& pt)
{
    return OWM->CalRho(Pin, Pbb, Tin, &Ziin[0], pt);
}



///////////////////////////////////////////////
// MixtureUnitBlkOil_OW
///////////////////////////////////////////////

MixtureUnitBlkOil_GW::MixtureUnitBlkOil_GW(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts)
{
    GWM = new OCPMixtureBlkOilGW(rs_param, i);

    vs          = &GWM->GetVarSet();
}

void MixtureUnitBlkOil_GW::Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{
    GWM->Flash(Pin, Tin, Niin);
}

void MixtureUnitBlkOil_GW::InitFlashIMPEC(const OCP_USI& bId, const BulkVarSet& bvs)
{
    GWM->InitFlash(bId, bvs);
}

void MixtureUnitBlkOil_GW::InitFlashFIM(const OCP_USI& bId, const BulkVarSet& bvs)
{
    GWM->InitFlashDer(bId, bvs);
}

void MixtureUnitBlkOil_GW::FlashIMPEC(const OCP_USI& bId, const BulkVarSet& bvs)
{
    GWM->Flash(bId, bvs);
}

void MixtureUnitBlkOil_GW::FlashFIM(const OCP_USI& bId, const BulkVarSet& bvs)
{
    GWM->FlashDer(bId, bvs);
}

OCP_DBL MixtureUnitBlkOil_GW::XiPhase(const OCP_DBL& Pin,
    const OCP_DBL& Tin,
    const vector<OCP_DBL>& Ziin,
    const PhaseType& pt)
{
    return GWM->CalXi(Pin, Pin, Tin, &Ziin[0], pt);
}

OCP_DBL MixtureUnitBlkOil_GW::RhoPhase(const OCP_DBL& Pin,
    const OCP_DBL& Pbb,
    const OCP_DBL& Tin,
    const vector<OCP_DBL>& Ziin,
    const PhaseType& pt)
{
    return GWM->CalRho(Pin, Pbb, Tin, &Ziin[0], pt);
}




///////////////////////////////////////////////
// MixtureUnitBlkOil_OGW
///////////////////////////////////////////////

MixtureUnitBlkOil_OGW::MixtureUnitBlkOil_OGW(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts)
{
    OGWM = new OCPMixtureBlkOilOGW(rs_param, i);
    vs          = &OGWM->GetVarSet();
}

void MixtureUnitBlkOil_OGW::Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{

    OGWM->Flash(Pin, Tin, Niin);
}


void MixtureUnitBlkOil_OGW::InitFlashIMPEC(const OCP_USI& bId, const BulkVarSet& bvs)
{
    OGWM->InitFlash(bId, bvs);
}

void MixtureUnitBlkOil_OGW::InitFlashFIM(const OCP_USI& bId, const BulkVarSet& bvs)
{
    OGWM->InitFlashDer(bId, bvs);
}

void MixtureUnitBlkOil_OGW::FlashIMPEC(const OCP_USI& bId, const BulkVarSet& bvs)
{
    OGWM->Flash(bId, bvs);
}

void MixtureUnitBlkOil_OGW::FlashFIM(const OCP_USI& bId, const BulkVarSet& bvs)
{
    OGWM->FlashDer(bId, bvs);
}

OCP_DBL
MixtureUnitBlkOil_OGW::XiPhase(const OCP_DBL& Pin,
    const OCP_DBL& Tin,
    const vector<OCP_DBL>& Ziin,
    const PhaseType& pt)
{
    return OGWM->CalXi(Pin, Pin, Tin, &Ziin[0], pt);
}

OCP_DBL
MixtureUnitBlkOil_OGW::RhoPhase(const OCP_DBL& Pin,
    const OCP_DBL& Pbbin,
    const OCP_DBL& Tin,
    const vector<OCP_DBL>& Ziin,
    const PhaseType& pt)
{
    return OGWM->CalRho(Pin, Pbbin, Tin, &Ziin[0], pt);
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/