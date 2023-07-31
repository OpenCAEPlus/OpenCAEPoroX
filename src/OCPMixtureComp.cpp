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

#include "OCPMixtureComp.hpp"


////////////////////////////////////////////////////////////////
// OCPMixtureCompMethod
////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////
// Baisc Components Property £¨Phase Equilibrium Calculations£©
////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////
// EoS and Phase Equilibrium Calculation
////////////////////////////////////////////////////////////////


void OCPMixtureCompMethod::CopyPhaseFromPE(OCPMixtureVarSet& vs)
{
    NP = PE.GetNP();
    ftype = PE.GetFtype();
    flagSkip = PE.GetFlagSkip();
    for (USI j = 0; j < NP; j++) {
        vs.nj[j] = Nh * PE.GetNu(j);
        copy(PE.GetX(j).begin(), PE.GetX(j).end(), &vs.x[j * vs.nc]);
    }
}


////////////////////////////////////////////////////////////////
// Baisc Phases Property £¨Phase Equilibrium Calculations£©
////////////////////////////////////////////////////////////////


void OCPMixtureCompMethod::CalMW(OCPMixtureVarSet& vs)
{
    for (USI j = 0; j < NP; j++) {
        MW[j] = 0;
        for (USI i = 0; i < NC; i++) {
            MW[j] += vs.x[j * vs.nc + i] * MWC[i];
        }
    }
}

void OCPMixtureCompMethod::CalVmVj(OCPMixtureVarSet& vs)
{
    for (USI j = 0; j < NP; j++) {
        vm[j] = eos.CalVmDer(vs.P, vs.T, &vs.x[j * vs.nc], vmP[j], &vmx[j][0]);
        vs.vj[j] = vs.nj[j] * vm[j];
    }
}

void OCPMixtureCompMethod::CalProperty(OCPMixtureVarSet& vs)
{
    CopyPhaseFromPE(vs);
    CalMW(vs);
    CalVmVj(vs);
    CalXiRhoMu(vs);
    IdentifyPhase(vs);
    ReOrderPhase(vs);
}

void OCPMixtureCompMethod::CalXiRhoMu(OCPMixtureVarSet& vs)
{
    for (USI j = 0; j < NP; j++) {
        vs.xi[j] = 1 / vm[j];
        vs.rho[j] = MW[j] * vs.xi[j];
        vs.mu[j] = visCal.CalViscosity(ViscosityParams(&vs.P, &vs.T, &vs.x[j * vs.nc], &vs.xi[j]));
    }
}

void OCPMixtureCompMethod::CalPropertyDer(OCPMixtureVarSet& vs)
{
    CopyPhaseFromPE(vs);
    CalMW(vs);
    CalVmVj(vs);
    CalXiRhoMuDer(vs);
    IdentifyPhase(vs);
    ReOrderPhaseDer(vs);
}

void OCPMixtureCompMethod::CalXiRhoMuDer(OCPMixtureVarSet& vs)
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


////////////////////////////////////////////////////////////////
// Derivatives Calculations
////////////////////////////////////////////////////////////////


void OCPMixtureCompMethod::CalVfiVfp_full01()
{

}


void OCPMixtureCompMethod::AssembleMatVfiVfp_full01()
{

}


void OCPMixtureCompMethod::AssembleRhsVfiVfp_full01()
{

}


void OCPMixtureCompMethod::CaldXsdXp01()
{

}


void OCPMixtureCompMethod::CalVfiVfp_full02()
{

}


void OCPMixtureCompMethod::AssembleMatVfiVfp_full02()
{

}


void OCPMixtureCompMethod::AssembleRhsVfiVfp_full02()
{


}


void OCPMixtureCompMethod::CaldXsdXp02()
{

}




///////////////////////////////////////////////
// Phase Identification
///////////////////////////////////////////////


void OCPMixtureCompMethod::IdentifyPhase(OCPMixtureVarSet& vs)
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
            phaseLabel[0]      = GAS;
            vs.phaseExist[GAS] = OCP_TRUE;
        }
        else {
            phaseLabel[0]      = OIL;
            vs.phaseExist[OIL] = OCP_TRUE;
        }
    }
    else {
        // Compare MW
        if (MW[0] > MW[1]) {
            phaseLabel[0] = OIL;
            phaseLabel[1] = GAS;
        }
        else {
            phaseLabel[0] = GAS;
            phaseLabel[1] = OIL;
        }
        vs.phaseExist[OIL] = OCP_TRUE;
        vs.phaseExist[GAS] = OCP_TRUE;
    }

    epIndexH.clear();
    for (USI j = 0; j < NP; j++) {
        epIndexH.push_back(phaseLabel[j]);
    }
}


void OCPMixtureCompMethod::ReOrderPhase(OCPMixtureVarSet& vs)
{

}


void OCPMixtureCompMethod::ReOrderPhaseDer(OCPMixtureVarSet& vs)
{

}


///////////////////////////////////////////////
// Water Property
///////////////////////////////////////////////

void OCPMixtureCompMethod01::CalPropertyW(const OCP_DBL& vw, OCPMixtureVarSet& vs)
{
    // if vw >= 0(water volume is given), then correct Nw
    vs.phaseExist[wIdP] = OCP_TRUE;
    vs.x[vs.nc + wIdC] = 1.0;

    PVTW.CalRhoXiMuDer(vs.P, vs.rho[wIdP], vs.xi[wIdP], vs.mu[wIdP], vs.rhoP[wIdP], vs.xiP[wIdP], vs.muP[wIdP]);
    if (vw >= 0) {
        vs.Ni[wIdC] = vw * vs.xi[wIdP];
    }
    vs.nj[wIdP]        = vs.Ni[wIdC];
    vs.vj[wIdP]        = vs.nj[wIdP] / vs.xi[wIdP];
    vs.vji[wIdP][wIdC] = 1 / vs.xi[wIdP];
    vs.vjP[wIdP]       = -vs.nj[wIdP] * vs.xiP[wIdP] / (vs.xi[wIdP] * vs.xi[wIdP]);
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/31/2023      Create file                          */
/*----------------------------------------------------------------------------*/