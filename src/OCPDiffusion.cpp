/*! \file    OCPDiffusion.cpp
 *  \brief   OCPDiffusion class declaration
 *  \author  Shizhe Li
 *  \date    Feb/28/2024
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPDiffusion.hpp"


OCPDiffusionMethod01::OCPDiffusionMethod01(const ParamReservoir& rs_param, OCPDiffusionVarSet& dvs)
{
    bcd.Setup();
    // temp
    diffuCP.resize(2);
    diffuCP[dvs.g].resize(dvs.nc, 1.6E-1);
    diffuCP[dvs.w].resize(dvs.nc, 1.0E-5);
}


void OCPDiffusionMethod01::CalFlux(const BulkConnPair& bp, const OCPDiffusionVarSet& dvs, const BulkVarSet& bvs, FluxVarSet& fvs) const
{  
    const USI& np = fvs.np;
    const USI& nc = fvs.nc;
    auto& flux_ni = fvs.flux_ni;

    const OCP_USI& bId = bp.BId();
    const OCP_USI& eId = bp.EId();
    const OCP_DBL& Akd = bp.Diffu();
    
    OCP_USI       bId_np_j, eId_np_j, uId_np_j;
    OCP_BOOL      exbegin, exend;
    OCP_DBL       weight;

    for (USI j = 0; j < np; j++) {
        bId_np_j = bId * np + j;
        eId_np_j = eId * np + j;

        exbegin = bvs.phaseExist[bId_np_j];
        exend   = bvs.phaseExist[eId_np_j];

        if ((!exbegin) && (!exend))  continue;

        if ((exbegin) && (exend)) {
            for (USI i = 0; i < nc; i++) {
                if (bvs.xij[bId_np_j * nc + i] > bvs.xij[eId_np_j * nc + i])  uId_np_j = bId_np_j;
                else                                                          uId_np_j = eId_np_j;
                flux_ni[i] += Akd * diffuCP[j][i] * bvs.S[uId_np_j] * bvs.xi[uId_np_j] * (bvs.xij[bId_np_j * nc + i] - bvs.xij[eId_np_j * nc + i]);
            }
        }
        else {
			if (exbegin) {
				uId_np_j = bId_np_j;
				weight   = 1.0;
			}
			else {
				uId_np_j = eId_np_j;
				weight   = -1.0;
			}
			for (USI i = 0; i < nc; i++) {
				flux_ni[i] += Akd * diffuCP[j][i] * bvs.S[uId_np_j] * bvs.xi[uId_np_j] * bvs.xij[uId_np_j * nc + i] * weight;
			}
        }
    }
}


void OCPDiffusionMethod01::AssembleMatFIM(const BulkConnPair& bp, const OCPDiffusionVarSet& dvs, const BulkVarSet& bvs, FluxVarSet& fvs) const
{
    const USI& np     = fvs.np;
    const USI& nc     = fvs.nc;
    const USI& ncol1  = fvs.ncol1;
    const USI& ncol2  = fvs.ncol2;
    auto&      dFdXpB = fvs.dFdXpB;
    auto&      dFdXpE = fvs.dFdXpE;
    auto&      dFdXsB = fvs.dFdXsB;
    auto&      dFdXsE = fvs.dFdXsE;
    OCP_DBL*   dFdXpU;  // up    bulk: dF / dXp
    OCP_DBL*   dFdXsU;  // up    bulk: dF / dXs

    const OCP_USI& bId = bp.BId();
    const OCP_USI& eId = bp.EId();
    const OCP_DBL& Akd = bp.Diffu();

    OCP_USI  bId_np_j, eId_np_j, uId_np_j;
    OCP_BOOL exbegin, exend;
    OCP_DBL  weight;

    for (USI j = 0; j < np; j++) {
        bId_np_j = bId * np + j;
        eId_np_j = eId * np + j;

        exbegin = bvs.phaseExist[bId_np_j];
        exend   = bvs.phaseExist[eId_np_j];

        if ((!exbegin) && (!exend))  continue;

        if ((exbegin) && (exend)) {
            // Akd * diffuCP[j][i] * bvs.S[uId_np_j] * bvs.xi[uId_np_j] * (bvs.xij[bId_np_j * nc + i] - bvs.xij[eId_np_j * nc + i]);
            for (USI i = 0; i < nc; i++) {
                const OCP_DBL tmpC = Akd * diffuCP[j][i];
                const OCP_DBL dxij = bvs.xij[bId_np_j * nc + i] - bvs.xij[eId_np_j * nc + i];

                if (bvs.xij[bId_np_j * nc + i] > bvs.xij[eId_np_j * nc + i]) {
                    uId_np_j = bId_np_j;
                    dFdXpU   = &dFdXpB[0];
                    dFdXsU   = &dFdXsB[0];
                }
                else {
                    uId_np_j = eId_np_j;
                    dFdXpU   = &dFdXpE[0];
                    dFdXsU   = &dFdXsE[0];
                }

                // dP
                dFdXpU[(i + 1) * ncol1] += tmpC * bvs.S[uId_np_j] * bvs.xiP[uId_np_j] * dxij;

                // dS
                dFdXsU[(i + 1) * ncol2 + j] += tmpC * bvs.xi[uId_np_j] * dxij;

                // dx
                for (USI k = 0; k < nc; k++) {
                    dFdXsU[(i + 1) * ncol2 + np + j * nc + k] += tmpC * bvs.S[uId_np_j] * bvs.xix[uId_np_j * nc + k] * dxij;
                }
                // k = i
                dFdXsB[(i + 1) * ncol2 + np + j * nc + i] += tmpC * bvs.S[uId_np_j] * bvs.xi[uId_np_j];
                dFdXsE[(i + 1) * ncol2 + np + j * nc + i] -= tmpC * bvs.S[uId_np_j] * bvs.xi[uId_np_j];
            }
        }
        else {
            if (exbegin) {
                uId_np_j = bId_np_j;
                weight   = 1.0;
                dFdXpU   = &dFdXpB[0];
                dFdXsU   = &dFdXsB[0];
            }
            else {
                uId_np_j = eId_np_j;
                weight   = -1.0;
                dFdXpU   = &dFdXpE[0];
                dFdXsU   = &dFdXsE[0];
            }
            // Akd * diffuCP[j][i] * bvs.S[uId_np_j] * bvs.xi[uId_np_j] * (bvs.xij[bId_np_j * nc + i] - bvs.xij[eId_np_j * nc + i]);
            // Akd * diffuCP[j][i] * bvs.S[uId_np_j] * bvs.xi[uId_np_j] * bvs.xij[uId_np_j * nc + i] * weight;     
            for (USI i = 0; i < nc; i++) {

                const OCP_DBL tmpC = Akd * diffuCP[j][i];
                const OCP_DBL dxij = bvs.xij[uId_np_j * nc + i] * weight;

                // dP
                dFdXpU[(i + 1) * ncol1] += tmpC * bvs.S[uId_np_j] * bvs.xiP[uId_np_j] * dxij;

                // dS
                dFdXsU[(i + 1) * ncol2 + j] += tmpC * bvs.xi[uId_np_j] * dxij;

                // dx
                for (USI k = 0; k < nc; k++) {
                    dFdXsU[(i + 1) * ncol2 + np + j * nc + k] += tmpC * bvs.S[uId_np_j] * bvs.xix[uId_np_j * nc + k] * dxij;
                }
                // k = i
                dFdXsB[(i + 1) * ncol2 + np + j * nc + i] += tmpC * bvs.S[uId_np_j] * bvs.xi[uId_np_j];
                dFdXsE[(i + 1) * ncol2 + np + j * nc + i] -= tmpC * bvs.S[uId_np_j] * bvs.xi[uId_np_j];
            }
        }
    }
}


void OCPDiffusion::Setup(const ParamReservoir& rs_param, const BulkVarSet& bvs)
{
    if (ifUse) {
        vs.SetUp(bvs.np, bvs.nc, bvs.o, bvs.g, bvs.w);
        dM.push_back(new OCPDiffusionMethod01(rs_param, vs));
    }
}


void OCPDiffusion::CalDiffu(BulkConnPair& bp, const Bulk& bk)
{ 
    if (ifUse) {
        dM[0]->CalDiffu(bp, bk);
    }
}


void OCPDiffusion::CalFlux(const BulkConnPair& bp, const BulkVarSet& bvs, FluxVarSet& fvs)
{
    if (ifUse) {
        dM[0]->CalFlux(bp, vs, bvs, fvs);
    }
}


void OCPDiffusion::AssembleMatFIM(const BulkConnPair& bp, const BulkVarSet& bvs, FluxVarSet& fvs) const
{
    if (ifUse) {
        dM[0]->AssembleMatFIM(bp, vs, bvs, fvs);
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Feb/28/2024      Create file                          */
/*----------------------------------------------------------------------------*/