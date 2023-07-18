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


/////////////////////////////////////////////////////
// OCPMixtureBlkOilOWMethod
/////////////////////////////////////////////////////


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
// OCPMixtureBlkOilOWMethod01
/////////////////////////////////////////////////////


/// Use PVDO and PVTW
class OCPMixtureBlkOilOWMethod01 : public OCPMixtureBlkOilOWMethod
{
public:
    OCPMixtureBlkOilOWMethod01(const vector<vector<OCP_DBL>>& PVDOin,
                               const vector<vector<OCP_DBL>>& PVTWin,
                               const OCP_DBL& stdRhoO,
                               const OCP_DBL& stdRhoW,
                               OCPMixtureVarSet* vsin) {
        PVDO.Setup(PVDOin, stdRhoO);
        PVTW.Setup(PVTWin, stdRhoW);
        vs   = vsin;
    }


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
    void Setup(const ParamReservoir& rs_param, const USI& i);

protected:
    void GetStdRhoOW(const ParamReservoir& rs_param);
    void SetPN(const OCP_DBL& P, const OCP_DBL* Ni) { 
        vs.P = P;
        vs.Ni[0] = Ni[0];
        vs.Ni[1] = Ni[1];
    }
    void SetPN(const OCP_DBL& P, const OCP_DBL& Pb, const OCP_DBL& Sw) {
        vs.P     = P;
        vs.Pb    = Pb;
        vs.S[1]  = Sw;
    }

protected:
    OCPMixtureBlkOilOWMethod* pmMethod;
    OCP_DBL                   stdRhoO;
    OCP_DBL                   stdRhoW;
};


#endif /* end if __OCPMIXTUREBLKOILOW_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/13/2023      Create file                          */
/*----------------------------------------------------------------------------*/