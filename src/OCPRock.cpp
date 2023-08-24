/*! \file    OCPRock.cpp
 *  \brief   OCPRock class declaration
 *  \author  Shizhe Li
 *  \date    Nov/15/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// OpenCAEPoroX header files
#include "OCPRock.hpp"

///////////////////////////////////////////////
// OCPRock_Linear
///////////////////////////////////////////////

void OCPRockIsoT_Linear::CalPoro(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& poroInit, const BulkContent& bcType)
{
    const OCP_DBL dP = (P - Pref);
    poro             = poroInit * (1 + (cp1 + cp2 / 2 * dP) * dP);
    dPorodP          = poroInit * (cp1 + cp2 * dP);
}

///////////////////////////////////////////////
// OCPRockT
///////////////////////////////////////////////

inline void OCPRockT::CalRockHT(const OCP_DBL& T)
{
    const OCP_DBL Ta  = T + CONV5;
    const OCP_DBL Tra = Tref + CONV5;
    Hr                = hcp1 * (Ta - Tra) + 0.5 * hcp2 * (Ta * Ta - Tra * Tra);
    dHrdT             = hcp1 + hcp2 * Ta;
}

///////////////////////////////////////////////
// OCPRockT_Linear
///////////////////////////////////////////////

void OCPRockT_Linear ::CalPoro(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& poroInit, const BulkContent& bcType)
{
    if (bcType == BulkContent::rf) {
        // with fluid
        // calculate porosity
        const OCP_DBL dP = P - Pref;
        const OCP_DBL dT = T - Tref;
        poro = poroInit * (1 + (cp * dP - ct * dT + cpt * dP * dT));
        dPorodP = poroInit * (cp + cpt * dT);
        dPorodT = poroInit * (-ct + cpt * dP);

        if (ConstRock) {
            _poro = 1 - poroInit;
            _dPorodP = 0.0;
            _dPorodT = 0.0;
        }
        else {
            _poro = 1 - poro;
            _dPorodP = -dPorodP;
            _dPorodT = -dPorodT;
        }
    }

    // Calculate enthalpy of rock
    CalRockHT(T);
}

///////////////////////////////////////////////
// OCPRockT_Exp
///////////////////////////////////////////////

void OCPRockT_Exp::CalPoro(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& poroInit, const BulkContent& bcType)
{
    if (bcType == BulkContent::rf) {
        // with fluid
        const OCP_DBL dP = P - Pref;
        const OCP_DBL dT = T - Tref;
        poro = poroInit * exp(cp * dP - ct * dT + cpt * dP * dT);

        dPorodP = poro * (cp + cpt * dT);
        dPorodT = poro * (-ct + cpt * dP);
        if (ConstRock) {
            _poro = 1 - poroInit;
            _dPorodP = 0.0;
            _dPorodT = 0.0;
        }
        else {
            _poro = 1 - poro;
            _dPorodP = -dPorodP;
            _dPorodT = -dPorodT;
        }
    }

    CalRockHT(T);
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/15/2022      Create file                          */
/*----------------------------------------------------------------------------*/
