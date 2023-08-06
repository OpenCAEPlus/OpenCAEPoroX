/*! \file    AcceleratePVT.cpp
 *  \brief   AcceleratePVT class declaration
 *  \author  Shizhe Li
 *  \date    Dec/25/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "AcceleratePVT.hpp"


void SkipPSAVarset::ResetToLastTimeStep()
{
    flag     = lflag;
    minEigen = lminEigen;
    P        = lP;
    T        = lT;
    zi       = lzi;
}

void SkipPSAVarset::UpdateLastTimeStep()
{
    lflag     = flag;
    lminEigen = minEigen;
    lP        = P;
    lT        = T;
    lzi       = zi;
}


SkipPSAMethod01::SkipPSAMethod01(SkipPSAVarset* vsin)
{
    vs = vsin;
    if (!vs->ifSetup) {

        vs->ifSetup = OCP_TRUE;

        const OCP_USI& nb = vs->nb;
        const USI&     nc = vs->nc;

        zi.resize(nc);
      
        vs->flag.resize(nb);
        vs->P.resize(nb);
        vs->T.resize(nb);
        vs->minEigen.resize(nb);
        vs->zi.resize(nb * nc);

        vs->lflag.resize(nb);
        vs->lP.resize(nb);
        vs->lT.resize(nb);
        vs->lminEigen.resize(nb);
        vs->lzi.resize(nb * nc);
    }
}


void SkipPSAMethod01::SetPTZ(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{
    P = Pin;
    T = Tin;
    Nt = 0;
    for (USI i = 0; i < vs->nc; i++) {
        Nt += Niin[i];
    }
    for (USI i = 0; i < vs->nc; i++) {
        zi[i] = Niin[i] / Nt;
    }
}


OCP_BOOL SkipPSAMethod01::IfSkip(const OCP_USI& n) const
{
    if (vs->flag[n]) {
        if (fabs(1 - vs->P[n] / P) >= vs->minEigen[n] / 10) {
            return OCP_FALSE;
        }
        if (fabs(vs->T[n] - T) >= vs->minEigen[n] * 10) {
            return OCP_FALSE;
        }
        for (USI i = 0; i < vs->nc; i++) {
            if (fabs(zi[i] - vs->zi[n * vs->nc + i]) >= vs->minEigen[n] / 10) {
                return OCP_FALSE;
            }
        }
        return OCP_TRUE;
    }
    else {
        return OCP_FALSE;
    }
}


USI SkipPSAMethod01::CalFtype(const OCP_DBL& Pin,
                              const OCP_DBL& Tin,
                              const OCP_DBL* Niin,
                              const OCP_USI& n)
{
    SetPTZ(Pin, Tin, Niin);
	if (IfSkip(n)) {
		return 1;
	}
	else {
		return 0;
	}
}


USI SkipPSAMethod01::CalFtype(const OCP_DBL& Pin,
                     const OCP_DBL&          Tin,
                     const OCP_DBL*          Niin,
                     const OCP_DBL*          S,
                     const USI&              np,
                     const OCP_USI&          n)
{
    SetPTZ(Pin, Tin, Niin);
	if (IfSkip(n)) {
		return 1;
	}
	else if (np >= 2) {
        USI tmp = 0;
		for (USI j = 0; j < vs->np; j++) {
			if (S[j] >= 1E-4) {
                tmp++;
			}
		}
        // num of phases remains the same, then and flash from np phases directly
        // otherwise, flash from single phase
        if (tmp == np) return np;
        else           return 0;
	}
	else {
		return 0;
	}
}


void SkipPSAMethod01::CalSkipForNextStep(const OCP_USI& bId, const OCPPhaseEquilibrium& PE)
{
    const USI& ftype    = PE.GetFtype();
    EoSCalculation* eos = PE.GetEoS();

    if (ftype == 0) {
        // the range should be calculated, which also means last skip is unsatisfied
        const USI& nc = vs->nc;


        eos->CalLnPhiN(P, T, &zi[0], Nt, &lnphiN[0]);;

        for (USI i = 0; i < nc; i++) {
            for (USI j = 0; j <= i; j++) {
                skipMatSTA[i * nc + j] = delta(i, j) + Nt * sqrt(zi[i] * zi[j]) * lnphiN[i * nc + j];
                skipMatSTA[j * nc + i] = skipMatSTA[i * nc + j];
            }
        }

#ifdef DEBUG
        if (!CheckNan(skipMatSTA.size(), &skipMatSTA[0])) {
            OCP_WARNING("Nan in skipMatSTA!");
        }
#endif // DEBUG

        CalEigenSY(nc, &skipMatSTA[0], &eigenSkip[0], &eigenWork[0], 2 * nc + 1);
        vs->flag[bId]     = OCP_TRUE;
        vs->minEigen[bId] = eigenSkip[0];
        vs->P[bId]        = P;
        vs->T[bId]        = T;
        copy(zi.begin(), zi.end(), &vs->zi[bId * nc]);
    }
    else if (ftype == 1) {
        vs->flag[bId] = OCP_TRUE;
    }
    else if (ftype == 2) {
        vs->flag[bId] = OCP_FALSE;
    }
    else {
        OCP_ABORT("WRONG Ftype!");
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/25/2022      Create file                          */
/*----------------------------------------------------------------------------*/