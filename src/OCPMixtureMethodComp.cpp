/*! \file    OCPMixtureComp.cpp
 *  \brief   OCPMixtureComp class declaration
 *  \author  Shizhe Li
 *  \date    Jul/31/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPMixtureMethodComp.hpp"


////////////////////////////////////////////////////////////////
// OCPMixtureMethodComp
////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////
// Baisc Components Property
////////////////////////////////////////////////////////////////


void OCPMixtureMethodComp::Setup(const ComponentParam& param, const USI& tarId)
{
    NPmax = param.numPhase;
    NC    = param.numCom; 

    if (param.Tc.activity)   Tc = param.Tc.data[tarId];
    else                     OCP_ABORT("TCRIT hasn't been input!");
    if (param.Pc.activity)   Pc = param.Pc.data[tarId];
    else                     OCP_ABORT("PCRIT hasn't been input!");

    if (param.Vc.activity)   Vc = param.Vc.data[tarId];
    else if (param.Zc.activity) {
        const vector<OCP_DBL>& Zc = param.Zc.data[tarId];
        Vc.resize(NC);
        for (USI i = 0; i < NC; i++) {
            Vc[i] = GAS_CONSTANT * Zc[i] * Tc[i] / Pc[i];
        }
    }
    else                     OCP_ABORT("VCRIT or ZCRIT hasn't been input!");

    if (param.MW.activity)   MWC = param.MW.data[tarId];
    else                     OCP_ABORT("MW hasn't been input!");

    /// Setup modules
    eos.Setup(param, tarId);
    visCal.Setup(param, tarId);
    PE.Setup(param, tarId, &eos);

    /// Allocate memory
    zi.resize(NC);
    MW.resize(NPmax);
    vm.resize(NPmax);
    vmP.resize(NPmax);
    vmx.resize(NPmax);
    for (auto& v : vmx)
        v.resize(NC);
    phaseLabel.resize(NPmax);
    epIndex.resize(NPmax);
    rowork.resize(NC);

    pivot.resize(NPmax * NC);
    lnfugP.resize(NPmax);
    lnfugN.resize(NPmax);
    for (USI j = 0; j < NPmax; j++) {
        lnfugP[j].resize(NC);
        lnfugN[j].resize(NC * NC);
    }
    JmatDer.resize(NPmax * NC * NPmax  * NC);
    rhsDer.resize(NPmax * NC * (NC + 1));
}


////////////////////////////////////////////////////////////////
// EoS and Phase Equilibrium Calculation
////////////////////////////////////////////////////////////////


void OCPMixtureMethodComp::CopyPhaseFromPE(OCPMixtureVarSet& vs)
{
    NP = PE.GetNP();
    for (USI j = 0; j < NP; j++) {
        vs.nj[j] = Nt * PE.GetNu(j);
        copy(PE.GetX(j).begin(), PE.GetX(j).end(), &vs.x[j * vs.nc]);
    }
}


////////////////////////////////////////////////////////////////
// Baisc Phases Property £¨Phase Equilibrium Calculations£©
////////////////////////////////////////////////////////////////


void OCPMixtureMethodComp::CalMW(OCPMixtureVarSet& vs)
{
    for (USI j = 0; j < NP; j++) {
        MW[j] = 0;
        for (USI i = 0; i < NC; i++) {
            MW[j] += vs.x[j * vs.nc + i] * MWC[i];
        }
    }
}

void OCPMixtureMethodComp::CalVmVj(OCPMixtureVarSet& vs)
{
    for (USI j = 0; j < NP; j++) {
        vm[j]    = eos.CalVmDer(vs.P, vs.T, &vs.x[j * vs.nc], vmP[j], &vmx[j][0]);
        vs.vj[j] = vs.nj[j] * vm[j];
    }
}

void OCPMixtureMethodComp::CalProperty(OCPMixtureVarSet& vs)
{
    CopyPhaseFromPE(vs);
    CalMW(vs);
    CalVmVj(vs);
    CalXiRhoMu(vs);
    IdentifyPhase(vs);
    ReOrderPhase(vs);
}

void OCPMixtureMethodComp::CalXiRhoMu(OCPMixtureVarSet& vs)
{
    for (USI j = 0; j < NP; j++) {
        vs.xi[j] = 1 / vm[j];
        vs.rho[j] = MW[j] * vs.xi[j];
        vs.mu[j] = visCal.CalViscosity(ViscosityParams(&vs.P, &vs.T, &vs.x[j * vs.nc], &vs.xi[j]));
    }
}

void OCPMixtureMethodComp::CalPropertyDer(OCPMixtureVarSet& vs)
{
    CopyPhaseFromPE(vs);
    CalMW(vs);
    CalVmVj(vs);
    CalXiRhoMuDer(vs);
    IdentifyPhase(vs);
    ReOrderPhaseDer(vs);
}

void OCPMixtureMethodComp::CalXiRhoMuDer(OCPMixtureVarSet& vs)
{
    OCP_DBL   dummy;
    for (USI j = 0; j < NP; j++) {
        // molar density
        vs.xi[j] = 1 / vm[j];
        vs.xiP[j] = -vs.xi[j] * vs.xi[j] * vmP[j];
        for (USI i = 0; i < NC; i++) {
            vs.xix[j * vs.nc + i] = -vs.xi[j] * vs.xi[j] * vmx[j][i];
        }

        // mass density
        vs.rho[j] = MW[j] * vs.xi[j];
        vs.rhoP[j] = MW[j] * vs.xiP[j];
        for (USI i = 0; i < NC; i++) {
            vs.rhox[j * vs.nc + i] = vs.xix[j * vs.nc + i] * MW[j] + vs.xi[j] * MWC[i];
        }

        // viscosity
        vs.mu[j] = visCal.CalViscosity(ViscosityParams(&vs.P, &vs.T, &vs.x[j * vs.nc], &vs.xi[j], &vs.xiP[j], nullptr, &vs.xix[j * vs.nc]),
            vs.muP[j], dummy, &vs.mux[j * vs.nc]);
    }
}


///////////////////////////////////////////////
// Phase Identification
///////////////////////////////////////////////


void OCPMixtureMethodComp::IdentifyPhase(OCPMixtureVarSet& vs)
{
    vs.phaseExist[0] = OCP_FALSE;
    vs.phaseExist[1] = OCP_FALSE;
    if (NP == 1) {
        // Critical Temperature Method
        OCP_DBL A = 0;
        OCP_DBL B = 0;
        for (USI i = 0; i < NC; i++) {
            A += zi[i] * Vc[i] * Tc[i];
            B += zi[i] * Vc[i];
        }
        OCP_DBL Tc = A / B;
        if (vs.T > Tc) {
            phaseLabel[0]       = PhaseType::gas;
            vs.phaseExist[vs.g] = OCP_TRUE;
        }
        else {
            phaseLabel[0]       = PhaseType::oil;
            vs.phaseExist[vs.o] = OCP_TRUE;
        }
    }
    else {
        // Compare MW
        if (MW[0] > MW[1]) {
            phaseLabel[0] = PhaseType::oil;
            phaseLabel[1] = PhaseType::gas;
        }
        else {
            phaseLabel[0] = PhaseType::gas;
            phaseLabel[1] = PhaseType::oil;
        }
        vs.phaseExist[vs.o] = OCP_TRUE;
        vs.phaseExist[vs.g] = OCP_TRUE;
    }

    epIndex.clear();
    for (USI j = 0; j < NPmax; j++) {
        if (vs.phaseExist[j]) epIndex.push_back(j);
    }
}


void OCPMixtureMethodComp::ReOrderPhase(OCPMixtureVarSet& vs)
{
    // for NP <= 2 Now
    if (phaseLabel[0] != PhaseType::oil) {
        OCPSwap(&vs.nj[0], &vs.nj[1], 1, &rowork[0]);
        OCPSwap(&vs.x[0 * vs.nc], &vs.x[1 * vs.nc], NC, &rowork[0]);
        OCPSwap(&MW[0], &MW[1], 1, &rowork[0]);
        OCPSwap(&vs.vj[0], &vs.vj[1], 1, &rowork[0]);
        OCPSwap(&vm[0], &vm[1], 1, &rowork[0]);
        OCPSwap(&vmP[0], &vmP[1], 1, &rowork[0]);
        OCPSwap(&vmx[0][0], &vmx[1][0], NC, &rowork[0]);
        OCPSwap(&vs.xi[0], &vs.xi[1], 1, &rowork[0]);
        OCPSwap(&vs.rho[0], &vs.rho[1], 1, &rowork[0]);
        OCPSwap(&vs.mu[0], &vs.mu[1], 1, &rowork[0]);
    }
}


void OCPMixtureMethodComp::ReOrderPhaseDer(OCPMixtureVarSet& vs)
{
    // for NP <= 2 Now
    if (phaseLabel[0] != PhaseType::oil) {
        OCPSwap(&vs.nj[0], &vs.nj[1], 1, &rowork[0]);
        OCPSwap(&vs.x[0 * vs.nc], &vs.x[1 * vs.nc], NC, &rowork[0]);
        OCPSwap(&MW[0], &MW[1], 1, &rowork[0]);
        OCPSwap(&vs.vj[0], &vs.vj[1], 1, &rowork[0]);
        OCPSwap(&vm[0], &vm[1], 1, &rowork[0]);
        OCPSwap(&vmP[0], &vmP[1], 1, &rowork[0]);
        OCPSwap(&vmx[0][0], &vmx[1][0], NC, &rowork[0]);
        OCPSwap(&vs.xi[0], &vs.xi[1], 1, &rowork[0]);
        OCPSwap(&vs.xiP[0], &vs.xiP[1], 1, &rowork[0]);
        OCPSwap(&vs.xix[0 * vs.nc], &vs.xix[1 * vs.nc], NC, &rowork[0]);
        OCPSwap(&vs.rho[0], &vs.rho[1], 1, &rowork[0]);
        OCPSwap(&vs.rhoP[0], &vs.rhoP[1], 1, &rowork[0]);
        OCPSwap(&vs.rhox[0 * vs.nc], &vs.rhox[1 * vs.nc], NC, &rowork[0]);
        OCPSwap(&vs.mu[0], &vs.mu[1], 1, &rowork[0]);
        OCPSwap(&vs.muP[0], &vs.muP[1], 1, &rowork[0]);
        OCPSwap(&vs.mux[0 * vs.nc], &vs.mux[1 * vs.nc], NC, &rowork[0]);
    }
}


////////////////////////////////////////////////////////////////
// Derivatives Calculations
////////////////////////////////////////////////////////////////


void OCPMixtureMethodComp::CalVfiVfp_full01(OCPMixtureVarSet& vs)
{
    if (NP == 1) {
        const USI j = epIndex[0];
        vs.vjP[j] = vs.nj[j] * vmP[j];
        vs.vfP    = vs.vjP[j];
        for (USI i = 0; i < NC; i++) {
            vs.vji[j][i] = vm[j] + vmx[j][i];
            for (USI k = 0; k < NC; k++) {
                vs.vji[j][i] -= vmx[j][k] * vs.x[j * vs.nc + k];
            }
            vs.vfi[i] = vs.vji[j][i];
        }
    }
    else {
        // NP > 1
        for (const auto& j : epIndex) {
            eos.CalLnFugN(vs.P, vs.T, &vs.x[j * vs.nc], vs.nj[j], &lnfugN[j][0]);
            eos.CalLnFugP(vs.P, vs.T, &vs.x[j * vs.nc], &lnfugP[j][0]);
        }
        AssembleMatVfiVfp_full01();
        AssembleRhsVfiVfp_full01();
        LUSolve(NC + 1, NC * NP, &JmatDer[0], &rhsDer[0], &pivot[0]);

        // now d nm0 / dP(dNk) has been available
        const OCP_DBL* dnkjdP = &rhsDer[0];
        // Calculate Vfp
        vs.vfP = 0;
        vector<OCP_DBL> tmp(NP, 0);
        for (const auto& j : epIndex) {
            for (USI i = 0; i < NC; i++) {
                tmp[j] -= vmx[j][i] * vs.x[j * vs.nc + i];
            }
            vs.vjP[j] = vs.nj[j] * vmP[j];
            for (USI k = 0; k < NC; k++) {
                vs.vjP[j] += (vm[j] + tmp[j] + vmx[j][k]) * dnkjdP[k];
            }
            dnkjdP += NC;
            vs.vfP += vs.vjP[j];
        }

        // Calculate Vfi
        const OCP_DBL* dnkjdN = dnkjdP;
        for (USI i = 0; i < NC; i++) {
            vs.vfi[i] = 0;
            for (const auto& j : epIndex) {
                vs.vji[j][i] = 0;
                for (USI k = 0; k < NC; k++) {
                    vs.vji[j][i] += (vm[j] + tmp[j] + vmx[j][k]) * dnkjdN[k];
                }
                vs.vfi[i] += vs.vji[j][i];
                dnkjdN += NC;
            }
        }
    }

#ifdef OCP_NANCHECK
    if (!CheckNan(vs.vfi.size(), &vs.vfi[0])) {
        OCP_ABORT("INF or NAN in vfi !");
    }
    if (!CheckNan(1, &vs.vfP)) {
        OCP_ABORT("INF or NAN in vfP !");
    }
#endif // NANCHECK
}


void OCPMixtureMethodComp::AssembleMatVfiVfp_full01()
{
    fill(JmatDer.begin(), JmatDer.end(), 0.0);
    // Attention 1: JmatDer should be sorted by column
    // Attention 2: d ln fij / d nkj is symetric for each j;
    OCP_DBL* bId = &JmatDer[0];
    for (USI j = 0; j < NP - 1; j++) {
        // for jth phase
        OCP_DBL* fugNj = &lnfugN[epIndex[j]][0];
        for (USI i = 0; i < NC; i++) {
            // for ith components
            Dcopy(NC, bId, fugNj);
            bId    += (NP - 1 - j) * NC;
            bId[i] = 1.0;
            bId    += (1 + j) * NC;
            fugNj  += NC;
        }
        bId += NC;
    }
    // NP - 1 phase
    bId            = &JmatDer[(NP - 1) * (NP * NC * NC)];
    OCP_DBL* fugNj = &lnfugN[epIndex[NP - 1]][0];
    for (USI i = 0; i < NC; i++) {
        for (USI j = 0; j < NP - 1; j++) {
            Dcopy(NC, bId, fugNj);
            bId += NC;
        }
        Dscalar((NP - 1) * NC, -1.0, bId - (NP - 1) * NC);
        fugNj  += NC;
        bId[i] = 1.0;
        bId    += NC;
    }
}


void OCPMixtureMethodComp::AssembleRhsVfiVfp_full01()
{
    fill(rhsDer.begin(), rhsDer.end(), 0.0);
    OCP_DBL* rhstmp = &rhsDer[0];
    const USI refP  = epIndex[NP - 1];
    // d P
    for (USI j = 0; j < NP - 1; j++) {
        for (USI i = 0; i < NC; i++) {
            rhstmp[i] = lnfugP[refP][i] - lnfugP[epIndex[j]][i];
        }
        rhstmp += NC;
    }
    // d Nk
    for (USI k = 0; k < NC; k++) {       
        rhstmp[NC * (NP - 1) + k] = 1;
        rhstmp += NP * NC;
    }

#ifdef OCP_NANCHECK
    if (!CheckNan(rhsDer.size(), &rhsDer[0])) {
        OCP_ABORT("INF or NAN in rhsDer !");
    }
#endif // NANCHECK
}


void OCPMixtureMethodComp::CaldXsdXp01(OCPMixtureVarSet& vs)
{
    // Calculate Sj and Vf before
    // Note that vf is the total fluid volume
    // And all vji, vjP are available(for all phase, all components)
    
    // dS / dP
    // S = Sj, xij
    // P = P, Ni
    // water is included
    fill(vs.dXsdXp.begin(), vs.dXsdXp.end(), 0);
    const USI ncol = vs.nc + 1;

    // dS / dP, dS / dN
    for (const auto& j : epIndex) {
		OCP_DBL* bId = &vs.dXsdXp[j * ncol];
		// dS / dP
		bId[0] = (vs.vjP[j] - vs.vfP * vs.S[j]) / vs.Vf;
		bId++;
		// dS / dN
		for (USI m = 0; m < NC; m++) {
			bId[m] = (vs.vji[j][m] - vs.vfi[m] * vs.S[j]) / vs.Vf;
		}
    }

    if (NP == 1) {
        // dxij / dNm
        OCP_DBL* bId = &vs.dXsdXp[(vs.np + epIndex[0] * vs.nc) * ncol + 1];
        for (USI i = 0; i < NC; i++) {
            for (USI m = 0; m < NC; m++) {
                bId[m] = (delta(i, m) * Nt - vs.Ni[i]) / (Nt * Nt);
            }
            bId += ncol;
        }
    } 
    else {
        // NP >= 2
        OCP_DBL* const bId = &vs.dXsdXp[vs.np * ncol];
        // water component is only in water phase in current case
        // dxij / dP
        const OCP_DBL* dnkjdP = &rhsDer[0];
        for (const auto& j : epIndex) {
            OCP_DBL*  bId0     = bId + j * vs.nc * ncol;
            OCP_DBL   njDerSum = 0;
            for (USI i = 0; i < NC; i++) njDerSum += dnkjdP[i];
            for (USI i = 0; i < NC; i++) {
                bId0[0] = (dnkjdP[i]  - njDerSum * vs.x[j * vs.nc + i]) / vs.nj[j];
                bId0    += ncol;
            }
            dnkjdP += NC;
        }
        // dxij / dNm
        const OCP_DBL* dnkjdN = dnkjdP;
        for (USI m = 0; m < NC; m++) {
            for (const auto& j : epIndex) {
                OCP_DBL   njDerSum = 0;
                for (USI i = 0; i < NC; i++) njDerSum += dnkjdN[i];
                OCP_DBL* bId0 = bId + j * vs.nc * ncol + m + 1;
                for (USI i = 0; i < NC; i++) {
                    bId0[0] = (dnkjdN[i] - njDerSum * vs.x[j * vs.nc + i]) / vs.nj[j];
                    bId0    += ncol;
                }
                dnkjdN += NC;
            }
        }
    }
}


void OCPMixtureMethodComp::CalVfiVfp_full02(OCPMixtureVarSet& vs)
{
    // Attention!
    // NP = 1 or NP = 2

    if (NP == 1) {
        const USI j = epIndex[0];
        vs.vjP[j]   = vs.nj[j] * vmP[j];
        vs.vfP      = vs.vjP[j];
		for (USI i = 0; i < NC; i++) {
            vs.vji[j][i] = vm[j] + vmx[j][i];
			for (USI k = 0; k < NC; k++) {
                vs.vji[j][i] -= vmx[j][k] * vs.x[j * vs.nc + k];
			}
            vs.vfi[i] = vs.vji[j][i];
		}
    } 
    else if (NP == 2) {
        // JUST FOR NP = 2
        for (const auto& j : epIndex) {
            eos.CalLnFugN(vs.P, vs.T, &vs.x[j * vs.nc], vs.nj[j], &lnfugN[j][0]);
            eos.CalLnFugP(vs.P, vs.T, &vs.x[j * vs.nc], &lnfugP[j][0]);
        }

        AssembleMatVfiVfp_full02();
        AssembleRhsVfiVfp_full02();
        LUSolve(NC + 1, NC, &JmatDer[0], &rhsDer[0], &pivot[0]);

        // now d nm0 / dP(dNk) has been available
        OCP_DBL tmp0 = 0;
        OCP_DBL tmp1 = 0;
        for (USI i = 0; i < NC; i++) {
            tmp0 -= vmx[0][i] * vs.x[0 * vs.nc + i];
            tmp1 -= vmx[1][i] * vs.x[1 * vs.nc + i];
        }

        // vfP
        const OCP_DBL* dnkjdP = &rhsDer[0];
        vs.vjP[0] = vs.nj[0] * vmP[0];
        vs.vjP[1] = vs.nj[1] * vmP[1];
        for (USI k = 0; k < NC; k++) {
            vs.vjP[0] += (vm[0] + tmp0 + vmx[0][k]) * dnkjdP[k];
            vs.vjP[1] -= (vm[1] + tmp1 + vmx[1][k]) * dnkjdP[k];
        }
        vs.vfP = vs.vjP[0] + vs.vjP[1];

        // Calculate Vfi
        const OCP_DBL* dnkjdN = dnkjdP + NC;
        for (USI i = 0; i < NC; i++) {
            vs.vji[0][i] = 0;
            vs.vji[1][i] = 0;
            for (USI k = 0; k < NC; k++) {
                vs.vji[0][i] += (vm[0] + tmp0 + vmx[0][k]) * dnkjdN[k];
                vs.vji[1][i] += (vm[1] + tmp1 + vmx[1][k]) * (delta(i, k) - dnkjdN[k]);
            }
            vs.vfi[i] = vs.vji[0][i] + vs.vji[1][i];
            dnkjdN += NC;
        }
    }
    else {
        OCP_ABORT("USE CalVfiVfp_full01 !");
    }
}


void OCPMixtureMethodComp::AssembleMatVfiVfp_full02()
{
    // NP = 2
    fill(JmatDer.begin(), JmatDer.end(), 0.0);
    // Attention 1: JmatDer should be sorted by column
    // Attention 2: d ln fij / d nkj is symetric for each j;
    OCP_DBL* tmpMat = &JmatDer[0];
    for (USI i = 0; i < NC; i++) {
        for (USI m = 0; m < NC; m++) {
            tmpMat[m] = lnfugN[0][i * NC + m] + lnfugN[1][i * NC + m];
        }
        tmpMat += NC;
    }
}


void OCPMixtureMethodComp::AssembleRhsVfiVfp_full02()
{
    // NP = 2
    fill(rhsDer.begin(), rhsDer.end(), 0.0);
    OCP_DBL* rhstmp = &rhsDer[0];

    // dP
    for (USI i = 0; i < NC; i++) {
        rhstmp[i] = lnfugP[1][i] - lnfugP[0][i];
    }
    rhstmp += NC;

    // dNk
    for (USI k = 0; k < NC; k++) {
        for (USI i = 0; i < NC; i++) {
            rhstmp[i] = lnfugN[1][k * NC + i]; // d lnfij / d nkj = d lnfkj / d nij
        }
        rhstmp += NC;
    }
}


void OCPMixtureMethodComp::CaldXsdXp02(OCPMixtureVarSet& vs)
{
    // Calculate Sj and Vf before
    // Note that vf is the total fluid volume
    // And all vji, vjP are available
    // for phases and components that are in PE


    // Attention!
    // NP = 1 or NP = 2

    // dS / dP
    // S = Sj, xij
    // P = P, Ni
    // water is included
    fill(vs.dXsdXp.begin(), vs.dXsdXp.end(), 0);
    const USI     ncol = vs.nc + 1;

    for (const auto& j : epIndex) {
        OCP_DBL* bId = &vs.dXsdXp[j * ncol];
        // dS / dP
        bId[0] = (vs.vjP[j] - vs.vfP * vs.S[j]) / vs.Vf;
        bId++;
        // dS / dN
        for (USI m = 0; m < NC; m++) {
            bId[m] = (vs.vji[j][m] - vs.vfi[m] * vs.S[j]) / vs.Vf;
        }
    }

    if (NP == 1) {
        // dxij / dNm
		OCP_DBL* bId = &vs.dXsdXp[(vs.np + epIndex[0] * vs.nc) * ncol + 1];
		for (USI i = 0; i < NC; i++) {
			for (USI m = 0; m < NC; m++) {
				bId[m] = (delta(i, m) * Nt - vs.Ni[i]) / (Nt * Nt);
			}
			bId += ncol;
		}
    } 
    else if (NP == 2)  {
        // NP = 2
        OCP_DBL* bId            = &vs.dXsdXp[vs.np * ncol];
        // dxij / dP, dxij / dNm
        const OCP_DBL* dnkjdNP  = &rhsDer[0];
        const OCP_DBL* x0       = &vs.x[0 * vs.nc];
        const OCP_DBL* x1       = &vs.x[1 * vs.nc];
        OCP_DBL*       bId0     = bId + 0 * vs.nc * ncol;
        OCP_DBL*       bId1     = bId + 1 * vs.nc * ncol;
        OCP_DBL        njDerSum = 0;
        // dxij / dP
        for (USI i = 0; i < NC; i++) njDerSum += dnkjdNP[i];
        for (USI i = 0; i < NC; i++) {
            bId0[0] = (dnkjdNP[i] - njDerSum * x0[i]) / vs.nj[0];
            bId1[0] = (-dnkjdNP[i] + njDerSum * x1[i]) / vs.nj[1];
            bId0    += ncol;
            bId1    += ncol;
        }
        dnkjdNP += NC;
       
        // dxij / dNm
        for (USI m = 0; m < NC; m++) {
            njDerSum = 0;
            for (USI i = 0; i < NC; i++) njDerSum += dnkjdNP[i];
            bId0 = bId + 0 * vs.nc * ncol + m + 1;
            bId1 = bId + 1 * vs.nc * ncol + m + 1;
            for (USI i = 0; i < NC; i++) {
                bId0[0] = (dnkjdNP[i] - njDerSum * x0[i]) / vs.nj[0];
                bId1[0] = ((delta(i, m) - dnkjdNP[i]) - (1 - njDerSum) * x1[i]) / vs.nj[1];
                bId0    += ncol;
                bId1    += ncol;
            }
            dnkjdNP += NC;
        }
    }
    else {
        OCP_ABORT("USE CaldXsdXp01 !");
    }
}


////////////////////////////////////////////////////////////////
// OCPMixtureMethodComp01
////////////////////////////////////////////////////////////////


OCPMixtureMethodComp01::OCPMixtureMethodComp01(const ParamReservoir& rs_param, const USI& i, OCPMixtureVarSet& vs)
{
    vs.Init(OCPMixtureType::COMP, rs_param.comsParam.numPhase + 1, rs_param.comsParam.numCom + 1);

    Setup(rs_param.comsParam, i);
    // water property
    OCP_DBL std_RhoW;
    if (rs_param.PVTW_T.data.size() != 0) {
        if (rs_param.gravity.activity)
            std_RhoW = RHOW_STD * rs_param.gravity.data[1];
        if (rs_param.density.activity) std_RhoW = RHOW_STD;
        PVTW.Setup(rs_param.PVTW_T.data[i], std_RhoW, stdVw);
    }
    else {
        OCP_ABORT("PVTW is Missing!");
    }
    wIdP = vs.np - 1;
    wIdC = vs.nc - 1;
}


void OCPMixtureMethodComp01::SetVarSet(const OCP_USI& bId, const BulkVarSet& bvs, OCPMixtureVarSet& mvs) const
{
    mvs.P = bvs.P[bId];
    mvs.T = bvs.T[bId] + CONV5;
    copy(&bvs.Ni[bId * bvs.nc], &bvs.Ni[bId * bvs.nc] + bvs.nc, mvs.Ni.begin());
    copy(&bvs.S[bId * bvs.np], &bvs.S[bId * bvs.np] + bvs.np, mvs.S.begin());
    copy(&bvs.phaseExist[bId * bvs.np], &bvs.phaseExist[bId * bvs.np] + bvs.np, mvs.phaseExist.begin());
    copy(&bvs.xij[bId * bvs.np * bvs.nc], &bvs.xij[bId * bvs.np * bvs.nc] + bvs.np * bvs.nc, mvs.x.begin());
       
    mvs.phaseNum = 0;
    for (OCP_USI j = 0; j < mvs.np; j++) {
        if (mvs.phaseExist[j]) {
            mvs.phaseNum++;
        }
    }
}


void OCPMixtureMethodComp01::SetVarSet(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* Ni, OCPMixtureVarSet& mvs) const
{
    mvs.P = P;
    mvs.T = T + CONV5;
    copy(Ni, Ni + mvs.nc, mvs.Ni.begin());
}


void OCPMixtureMethodComp01::Flash(OCPMixtureVarSet& vs)
{
    InitNtZ(vs);
    PE.PhaseEquilibrium(vs.P, vs.T, &zi[0]);
    CalProperty(vs);
    CalPropertyW(-1.0, vs);
    vs.CalVfS();
}


void OCPMixtureMethodComp01::InitFlash(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
    InitNtZ(vs);
    Nt = 1;
    PE.PhaseEquilibrium(vs.P, vs.T, &zi[0]);
    CalPropertyDer(vs);
    CalPropertyW(Vp * vs.S[wIdP], vs);
    CorrectNt(Vp * (1 - vs.S[wIdP]), vs);
    vs.CalVfS();
    CalVfiVfp_full02(vs);
}


void OCPMixtureMethodComp01::Flash(OCPMixtureVarSet& vs, const USI& ftype)
{
    InitNtZ(vs);
    PE.PhaseEquilibrium(vs.P, vs.T, &zi[0], ftype, vs.phaseNum - 1, &vs.x[0], vs.nc);
    CalProperty(vs);
    CalPropertyW(-1.0, vs);
    vs.CalVfS();
    CalVfiVfp_full02(vs);  
}


void OCPMixtureMethodComp01::InitFlashDer(const OCP_DBL& Vp, OCPMixtureVarSet& vs)
{
    InitNtZ(vs);
    Nt = 1;
    PE.PhaseEquilibrium(vs.P, vs.T, &zi[0]);
    CalPropertyDer(vs);
    CalPropertyW(Vp * vs.S[wIdP], vs);
    CorrectNt(Vp * (1 - vs.S[wIdP]), vs);
    vs.CalVfS();
    CalVfiVfp_full02(vs);
    CaldXsdXp02(vs);
}


void OCPMixtureMethodComp01::FlashDer(OCPMixtureVarSet& vs, const USI& ftype)
{
    InitNtZ(vs);
    PE.PhaseEquilibrium(vs.P, vs.T, &zi[0], ftype, vs.phaseNum - 1, &vs.x[0], vs.nc);
    CalPropertyDer(vs);
    CalPropertyW(-1.0, vs);
    vs.CalVfS();
    CalVfiVfp_full02(vs);    
    CaldXsdXp02(vs);
}


OCP_DBL OCPMixtureMethodComp01::CalRho(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt)
{
    // assume that only single phase exists here
    if (pt == PhaseType::wat) {
        // water phase
        return PVTW.CalRhoW(P);
    }
    else {
        // hydrocarbon phase
        OCP_DBL xitmp = CalXi(P, T, &z[0], pt);
        OCP_DBL MWtmp = 0;
        for (USI i = 0; i < NC; i++) MWtmp += z[i] * MWC[i];
        return MWtmp * xitmp;
    }
}


OCP_DBL OCPMixtureMethodComp01::CalXi(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt)
{
    // assume that only single phase exists here
    if (pt == PhaseType::wat) {
        // water phase
        return PVTW.CalXiW(P);
    }
    else {
        // oil phase
        return 1 / eos.CalVm(P, T + CONV5, &z[0]);
    }
}


void OCPMixtureMethodComp01::CalVStd(OCPMixtureVarSet& vs)
{
    InitNtZ(vs);
    PE.PhaseEquilibrium(vs.P, vs.T, &zi[0]);
    CalProperty(vs);
    // Calculate water volume
    vs.vj[wIdP] = vs.Ni[wIdC] * stdVw * CONV1;
}


OCP_DBL OCPMixtureMethodComp01::CalVmStd(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL* z, const PhaseType& pt)
{
    if (pt == PhaseType::wat) {
        // water phase
        return (stdVw * CONV1);
    }
    else {
        // oil phase
        return eos.CalVm(P, T + CONV5, &z[0]);
    }
}


void OCPMixtureMethodComp01::InitNtZ(OCPMixtureVarSet& vs)
{
    Nt = 0;
    for (USI i = 0; i < NC; i++)    Nt += vs.Ni[i];
    for (USI i = 0; i < NC; i++)    zi[i] = vs.Ni[i] / Nt;
    vs.Nt = Nt + vs.Ni[wIdC];
}


void OCPMixtureMethodComp01::CorrectNt(const OCP_DBL& vh, OCPMixtureVarSet& vs)
{
    // correct Nt
    OCP_DBL vtmp = 0;
    for (const auto& j : epIndex) {
        vtmp += vs.vj[j];
    }
    Nt = vh / vtmp;
    for (const auto& j : epIndex) {
        vs.nj[j] *= Nt;
        vs.vj[j] *= Nt;
    }
    for (USI i = 0; i < NC; i++) {
        vs.Ni[i] = zi[i] * Nt;
    }
    vs.Nt = Nt + vs.Ni[wIdC];
}


///////////////////////////////////////////////
// Water Property
///////////////////////////////////////////////



void OCPMixtureMethodComp01::CalPropertyW(const OCP_DBL& vw, OCPMixtureVarSet& vs)
{
    // water property
    vs.phaseExist[wIdP]       = OCP_TRUE;
    vs.x[wIdP * vs.nc + wIdC] = 1.0;
    // if vw >= 0(water volume is given), then correct Nw
    PVTW.CalRhoXiMuDer(vs.P, vs.rho[wIdP], vs.xi[wIdP], vs.mu[wIdP], vs.rhoP[wIdP], vs.xiP[wIdP], vs.muP[wIdP]);
    if (vw >= 0) {
        vs.Ni[wIdC] = vw * vs.xi[wIdP];
    }
    vs.nj[wIdP]        = vs.Ni[wIdC];
    vs.vj[wIdP]        = vs.nj[wIdP] / vs.xi[wIdP];
    vs.vji[wIdP][wIdC] = 1 / vs.xi[wIdP];
    vs.vjP[wIdP]       = -vs.nj[wIdP] * vs.xiP[wIdP] / (vs.xi[wIdP] * vs.xi[wIdP]);
}


////////////////////////////////////////////////////////////
// Derivatives Calculations
////////////////////////////////////////////////////////////


void OCPMixtureMethodComp01::CalVfiVfp_full01(OCPMixtureVarSet& vs)
{
    OCPMixtureMethodComp::CalVfiVfp_full01(vs);
    AddVfpVfiW(vs);
}


void OCPMixtureMethodComp01::CaldXsdXp01(OCPMixtureVarSet& vs)
{
    OCPMixtureMethodComp::CaldXsdXp01(vs);
    CaldXsdXpW(vs);
}


void OCPMixtureMethodComp01::CalVfiVfp_full02(OCPMixtureVarSet& vs)
{
    OCPMixtureMethodComp::CalVfiVfp_full02(vs);
    AddVfpVfiW(vs);
}


void OCPMixtureMethodComp01::CaldXsdXp02(OCPMixtureVarSet& vs)
{
    OCPMixtureMethodComp::CaldXsdXp02(vs);
    CaldXsdXpW(vs);
}


void OCPMixtureMethodComp01::AddVfpVfiW(OCPMixtureVarSet& vs)
{
    // for water
    vs.vfP       += vs.vjP[wIdP];
    vs.vfi[wIdC] =  vs.vji[wIdP][wIdC];
}


void OCPMixtureMethodComp01::CaldXsdXpW(OCPMixtureVarSet& vs)
{
    // for water
    // d Sj / dNw
    const USI ncol = vs.nc + 1;   
    for (const auto& j : epIndex) {
        OCP_DBL* bId = &vs.dXsdXp[j* ncol + wIdC + 1];
        bId[0]       = (vs.vji[j][wIdC] - vs.vfi[wIdC] * vs.S[j]) / vs.Vf;
    }
    // d Sw / d (P, Ni, Nw)
    OCP_DBL* bId = &vs.dXsdXp[wIdP * ncol];
    bId[0]       = (vs.vjP[wIdP] - vs.vfP * vs.S[wIdP]) / vs.Vf;
    bId++;
    for (USI i = 0; i < vs.nc; i++) {
        bId[i] = (vs.vji[wIdP][i] - vs.vfi[i] * vs.S[wIdP]) / vs.Vf;
    }
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/31/2023      Create file                          */
/*----------------------------------------------------------------------------*/