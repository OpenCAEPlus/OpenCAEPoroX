/*! \file    OCPMixtureBlkoilOW.hpp
 *  \brief   OCPMixtureBlkoilOW class declaration
 *  \author  Shizhe Li
 *  \date    Jul/13/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPMIXTUREBLKOILOW_HEADER__
#define __OCPMIXTUREBLKOILOW_HEADER__

#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "OCPFuncPVT.hpp"
#include "OCPMixture.hpp"

#include <vector>

using namespace std;


/// Calculate oil, gas, water relative permeability and capillary pressure
class OCPMixtureBlkOilOWMethod
{
public:
    OCPMixtureBlkOilOWMethod() = default;

protected:
    enum phaseIndex { oil, wat };
    enum componentIndex { OIL, WAT };
    OCPMixtureVarSet* vs;
};


/////////////////////////////////////////////////////
// OCPOWFMethod01
/////////////////////////////////////////////////////


/// Use SWOF
class OCPMixtureBlkOilOWMethod01 : public OCPMixtureBlkOilOWMethod
{
public:
    OCPMixtureBlkOilOWMethod01(const vector<vector<OCP_DBL>>& PVDOin, OCPMixtureVarSet* vsin);

protected:
    
    OCP_PVDO        PVDO;
    OCP_PVTW        PVTW;
};


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOW 
/////////////////////////////////////////////////////

class OCPMixtureBlkOilOW : public OCPMixture
{
public:
    OCPMixtureBlkOilOW() { mixtureType = OCPMIXTURE_BO_OW; }


protected:
    OCPMixtureBlkOilOWMethod* pmMethod;
};


#endif /* end if __OCPMIXTUREBLKOILOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/13/2023      Create file                          */
/*----------------------------------------------------------------------------*/