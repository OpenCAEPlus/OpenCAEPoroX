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

#include "AcceleratePEC.hpp"


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


SkipPSAMethod01::SkipPSAMethod01(SkipPSAVarset& svs, const OCPMixtureCompMethod* compMin)
{
    compM  = compMin;
    svs.np = compM->GetNPmax();
    svs.nc = compM->GetNC();

    if (!svs.ifSetup) {

        svs.ifSetup = OCP_TRUE;

        svs.flag.resize(svs.nb);
        svs.P.resize(svs.nb);
        svs.T.resize(svs.nb);
        svs.minEigen.resize(svs.nb);
        svs.zi.resize(svs.nb * svs.nc);

        svs.lflag.resize(svs.nb);
        svs.lP.resize(svs.nb);
        svs.lT.resize(svs.nb);
        svs.lminEigen.resize(svs.nb);
        svs.lzi.resize(svs.nb * svs.nc);
    }

    lnphiN.resize(svs.nc * svs.nc);
    skipMatSTA.resize(svs.nc * svs.nc);
    eigenSkip.resize(svs.nc);
    eigenWork.resize(2 * svs.nc + 1);
}


OCP_BOOL SkipPSAMethod01::IfSkip(const OCP_USI& bId, const SkipPSAVarset& svs, const OCPMixtureVarSet& mvs) const
{
    if (svs.flag[bId]) {

        OCP_DBL Nt = 0;
        for (USI i = 0; i < svs.nc; i++) {
            Nt += mvs.Ni[i];
        }

        if (fabs(1 - svs.P[bId] / mvs.P) >= svs.minEigen[bId] / 10) {
            return OCP_FALSE;
        }
        if (fabs(svs.T[bId] - mvs.T) >= svs.minEigen[bId] * 10) {
            return OCP_FALSE;
        }
        for (USI i = 0; i < svs.nc; i++) {
            if (fabs(mvs.Ni[i] / Nt - svs.zi[bId * svs.nc + i]) >= svs.minEigen[bId] / 10) {
                return OCP_FALSE;
            }
        }
        return OCP_TRUE;
    }
    else {
        return OCP_FALSE;
    }
}


USI SkipPSAMethod01::CalFtype01(const OCP_USI& bId, const SkipPSAVarset& svs, const OCPMixtureVarSet& mvs)
{
	if (IfSkip(bId, svs, mvs)) {
		return 1;
	}
	else {
		return 0;
	}
}


USI SkipPSAMethod01::CalFtype02(const OCP_USI& bId, const SkipPSAVarset& svs, const OCPMixtureVarSet& mvs, const USI& np)
{
    const USI& npPE = compM->GetNumPhasePE(np);
	if (IfSkip(bId, svs, mvs)) {
		return 1;
	}
	else if (npPE >= 2) {
        USI tmp = 0;
		for (USI j = 0; j < svs.np; j++) {
			if (mvs.S[j] >= 1E-4) {
                tmp++;
			}
		}
        // num of phases remains the same, then and flash from np phases directly
        // otherwise, flash from single phase
        if (tmp == npPE) return npPE;
        else             return 0;
	}
	else {
		return 0;
	}
}


void SkipPSAMethod01::CalSkipForNextStep(const OCP_USI& bId, SkipPSAVarset& svs)
{
    const USI& ftype    = compM->GetFtype();
    

    if (ftype == 0) {
        // the range should be calculated, which also means last skip is unsatisfied
        const EoSCalculation*  eos = compM->GetEoS();
        const OCP_DBL&         P   = compM->GetP();
        const OCP_DBL&         T   = compM->GetT();
        const OCP_DBL&         Nt  = compM->GetNt();
        const vector<OCP_DBL>& zi  = compM->GetZi();
        const USI&             nc  = svs.nc;

        eos->CalLnPhiN(P, T, &zi[0], Nt, &lnphiN[0]);

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
        svs.flag[bId]     = OCP_TRUE;
        svs.minEigen[bId] = eigenSkip[0];
        svs.P[bId]        = P;
        svs.T[bId]        = T;
        copy(zi.begin(), zi.end(), &svs.zi[bId * nc]);
    }
    else if (ftype == 1) {
        svs.flag[bId] = OCP_TRUE;
    }
    else if (ftype == 2) {
        svs.flag[bId] = OCP_FALSE;
    }
    else {
        OCP_ABORT("WRONG Ftype!");
    }
}


USI SkipPSA::Setup(const OCP_USI& nb, const OCPMixtureCompMethod* compM)
{
    if (ifUse) {
        vs.SetNb(nb);
        sm.push_back(new SkipPSAMethod01(vs, compM));
        return sm.size() - 1;
    }
    return 0;
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/25/2022      Create file                          */
/*----------------------------------------------------------------------------*/