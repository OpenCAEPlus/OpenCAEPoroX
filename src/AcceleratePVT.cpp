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


SkipPSAMethod01::SkipPSAMethod01(SkipPSAVarset* vsin, const OCPMixtureComp* compsin)
{
    vs     = vsin;
    comps  = compsin;   
    vs->np = comps->GetNPmaxPE();
    vs->nc = comps->GetNCPE();

    if (!vs->ifSetup) {

        vs->ifSetup = OCP_TRUE;

        vs->flag.resize(vs->nb);
        vs->P.resize(vs->nb);
        vs->T.resize(vs->nb);
        vs->minEigen.resize(vs->nb);
        vs->zi.resize(vs->nb * vs->nc);

        vs->lflag.resize(vs->nb);
        vs->lP.resize(vs->nb);
        vs->lT.resize(vs->nb);
        vs->lminEigen.resize(vs->nb);
        vs->lzi.resize(vs->nb * vs->nc);
    }

    lnphiN.resize(vs->nc * vs->nc);
    skipMatSTA.resize(vs->nc * vs->nc);
    eigenSkip.resize(vs->nc);
    eigenWork.resize(2 * vs->nc + 1);
}


OCP_BOOL SkipPSAMethod01::IfSkip(const OCP_DBL& Pin,
                                 const OCP_DBL& Tin,
                                 const OCP_DBL* Niin,
                                 const OCP_USI& bId) const
{
    if (vs->flag[bId]) {

        OCP_DBL Nt = 0;
        for (USI i = 0; i < vs->nc; i++) {
            Nt += Niin[i];
        }

        if (fabs(1 - vs->P[bId] / Pin) >= vs->minEigen[bId] / 10) {
            return OCP_FALSE;
        }
        if (fabs(vs->T[bId] - (Tin + CONV5)) >= vs->minEigen[bId] * 10) {
            return OCP_FALSE;
        }
        for (USI i = 0; i < vs->nc; i++) {
            if (fabs(Niin[i] / Nt - vs->zi[bId * vs->nc + i]) >= vs->minEigen[bId] / 10) {
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
                              const OCP_USI& bId)
{
	if (IfSkip(Pin, Tin, Niin, bId)) {
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
                     const OCP_USI&          bId)
{
	if (IfSkip(Pin, Tin, Niin, bId)) {
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


void SkipPSAMethod01::CalSkipForNextStep(const OCP_USI& bId)
{
    const USI& ftype    = comps->GetFtypePE();
    

    if (ftype == 0) {
        // the range should be calculated, which also means last skip is unsatisfied
        EoSCalculation* eos       = comps->GetEoSPE();
        const OCP_DBL&  P         = comps->GetP();
        const OCP_DBL&  T         = comps->GetT();
        const OCP_DBL&  Nt        = comps->GetNtPE();
        const vector<OCP_DBL>& zi = comps->GetZiPE();
        const USI&      nc  = vs->nc;

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


USI SkipPSA::Setup(const OCP_USI& nb, const OCPMixtureComp* compsin)
{
    if (ifUse) {
        vs.SetNb(nb);
        sm.push_back(new SkipPSAMethod01(&vs, compsin));
        return sm.size() - 1;
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/25/2022      Create file                          */
/*----------------------------------------------------------------------------*/