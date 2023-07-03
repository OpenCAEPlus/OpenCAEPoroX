/*! \file    PhasePermeability.cpp
 *  \brief   PhasePermeability class declaration
 *  \author  Shizhe Li
 *  \date    Dec/26/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPMiscible.hpp"



OCP_DBL SurfaceTensionMethod01::CalSurfaceTension()
{
    if (*NP == 1)
        return 100;
    else {
        const OCP_DBL Bv = *bv * CONV7;
        const OCP_DBL Bl = *bl * CONV7;
        OCP_DBL surTen = 0;
        for (USI i = 0; i < NC; i++)
            surTen += parachor[i] * (Bv * xv[i] - Bl * xl[i]);
        return pow(surTen, 4.0);
    }
}


/////////////////////////////////////////////////////////////////////
// Miscible For Compositional Model
/////////////////////////////////////////////////////////////////////

void Miscible::InputParam(const Miscstr& misterm)
{
    const USI len = misterm.surTenRef.size();
    if (len > 0) {
        ifMiscible = OCP_TRUE;
        surTenRef     = misterm.surTenRef[0];
        surTenPc      = 1;
        if (len > 2) {
            surTenPc = misterm.surTenRef[2] / surTenRef;
        }
        Fkexp = 0.25;
    }
}

USI Miscible::Setup(const OCP_USI& numBulk, const SurfaceTensionMethodParams& params)
{
    if (ifMiscible) {
        surTen.resize(numBulk);
        Fk.resize(numBulk);
        Fp.resize(numBulk);

        // Last time step
        lsurTen.resize(numBulk);

        surfaceTensionMethod.push_back(new SurfaceTensionMethod01(params));

        return surfaceTensionMethod.size() - 1;
    }
    return 0;
}


void Miscible::CalSurfaceTension(const OCP_USI& bId, const USI& mIndex)
{
    if (ifMiscible) {
        surTen[bId] = surfaceTensionMethod[mIndex]->CalSurfaceTension();
    }
}


OCP_BOOL Miscible::CalFkFp(const OCP_USI& n, OCP_DBL& fk, OCP_DBL& fp)
{
    if (surTen[n] >= surTenRef || surTen[n] <= TINY) {
        Fk[n] = 1;
        Fp[n] = 1;
        return OCP_FALSE; // InMiscible
    } else {
        Fk[n] = fk = min(1.0, pow(surTen[n] / surTenRef, Fkexp));
        Fp[n] = fp = min(surTenPc, surTen[n] / surTenRef);
        return OCP_TRUE; // Miscible
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/26/2022      Create file                          */
/*----------------------------------------------------------------------------*/