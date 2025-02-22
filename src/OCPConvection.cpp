/*! \file    OCPConvection.cpp
 *  \brief   OCPConvection class declaration
 *  \author  Shizhe Li
 *  \date    May/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


// OpenCAEPoroX header files
#include "OCPConvection.hpp"


////////////////////////////////////////////
// OCPConvection01
////////////////////////////////////////////


void OCPConvection01::CalFlux(const BulkConnPair& bp, const Bulk& bk, FluxVarSet& fvs)
{
	// Calculte upblock, rho, vj, flux_ni

    const BulkVarSet& bvs = bk.vs;

    auto& flux_ni = fvs.flux_ni;

    const OCP_USI bId = bp.BId();
    const OCP_USI eId = bp.EId();
    const OCP_DBL Akd = bp.Trans(); 
    OCP_USI       bId_np_j, eId_np_j, uId_np_j;
    OCP_BOOL      exbegin, exend;
    OCP_DBL       rho;

    for (USI j = 0; j < np; j++) {
        bId_np_j = bId * np + j;
        eId_np_j = eId * np + j;

        exbegin  = bvs.phaseExist[bId_np_j];
        exend    = bvs.phaseExist[eId_np_j];

        if ((exbegin) && (exend)) {
            rho = (bvs.rho[bId_np_j] + bvs.rho[eId_np_j]) / 2;
        }
        else if (exbegin && (!exend)) {
            rho = bvs.rho[bId_np_j];
        }
        else if ((!exbegin) && (exend)) {
            rho = bvs.rho[eId_np_j];
        }
        else {
            upblock[j] = bId;
            dP[j]      = 0;
            vj[j] = 0;
            continue;
        }

        upblock[j] = bId;
        dP[j] = (bvs.Pj[bId_np_j] - GRAVITY_FACTOR * rho * bvs.depth[bId]) -
            (bvs.Pj[eId_np_j] - GRAVITY_FACTOR * rho * bvs.depth[eId]);

        if (dP[j] < 0)  upblock[j] = eId;

        uId_np_j = upblock[j] * np + j;

        if (!bvs.phaseExist[uId_np_j]) {
            vj[j] = 0;
            continue;
        }

        vj[j] = Akd * bvs.kr[uId_np_j] / bvs.mu[uId_np_j] * dP[j];
       
        for (USI i = 0; i < nc; i++) {
            flux_ni[i] += vj[j] * bvs.xi[uId_np_j] * bvs.xij[uId_np_j * nc + i];
        }
    }
}


void OCPConvection01::AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk, FluxVarSet& fvs)
{
    const USI& ncol1  = fvs.ncol1;
    const USI& ncol2  = fvs.ncol2;
    auto&      dFdXpB = fvs.dFdXpB;
    auto&      dFdXpE = fvs.dFdXpE;
    auto&      dFdXsB = fvs.dFdXsB;
    auto&      dFdXsE = fvs.dFdXsE;

    const BulkVarSet& bvs = bk.vs;

    const OCP_USI bId    = bp.BId();
    const OCP_USI eId    = bp.EId();
    const OCP_DBL Akd    = bp.Trans();
    const OCP_DBL dGamma = GRAVITY_FACTOR * (bvs.depth[bId] - bvs.depth[eId]);

    OCP_USI  bId_np_j, eId_np_j, uId_np_j, dId_np_j;
    OCP_BOOL phaseExistDj;
    OCP_DBL  rhoWghtU, rhoWghtD;
    OCP_DBL  dP, transJ, transIJ, kr, mu, xi, xij, xiP, muP, rhox, xix, mux, tmp;

    OCP_DBL* dFdXpU;     // up    bulk: dF / dXp
    OCP_DBL* dFdXpD;     // down  bulk: dF / dXp
    OCP_DBL* dFdXsU;     // up    bulk: dF / dXs
    OCP_DBL* dFdXsD;     // down  bulk: dF / dXs

    for (USI j = 0; j < np; j++) {
        uId_np_j = bcvs.upblock[c * np + j] * np + j;
        if (!bvs.phaseExist[uId_np_j]) continue;
        bId_np_j     = bId * np + j;
        eId_np_j     = eId * np + j;

        if (bId_np_j == uId_np_j) {
            dId_np_j     = eId_np_j;
            phaseExistDj = bvs.phaseExist[dId_np_j];                      
            dFdXpU       = &dFdXpB[0];
            dFdXpD       = &dFdXpE[0];
            dFdXsU       = &dFdXsB[0];
            dFdXsD       = &dFdXsE[0];            
        }
        else {
            dId_np_j     = bId_np_j;
            phaseExistDj = bvs.phaseExist[dId_np_j];                      
            dFdXpU       = &dFdXpE[0];
            dFdXpD       = &dFdXpB[0];
            dFdXsU       = &dFdXsE[0];
            dFdXsD       = &dFdXsB[0];            
        }
        if (phaseExistDj) {
            rhoWghtU = 0.5;
            rhoWghtD = 0.5;
        }
        else {
            rhoWghtU = 1;
            rhoWghtD = 0;
        }

        dP      = bcvs.dP[c * np + j];
        xi      = bvs.xi[uId_np_j];
        kr      = bvs.kr[uId_np_j];
        mu      = bvs.mu[uId_np_j];
        muP     = bvs.muP[uId_np_j];
        xiP     = bvs.xiP[uId_np_j];
        transJ  = Akd * kr / mu;

        for (USI i = 0; i < nc; i++) {
            xij    = bvs.xij[uId_np_j * nc + i];
            transIJ = xij * xi * transJ;

            // dP
            dFdXpB[(i + 1) * ncol1] += transIJ;
            dFdXpE[(i + 1) * ncol1] -= transIJ;

            tmp = transJ * xiP * xij * dP;
            tmp += -transIJ * muP / mu * dP;
            dFdXpU[(i + 1) * ncol1] +=
                (tmp - transIJ * rhoWghtU * bvs.rhoP[uId_np_j] * dGamma);
            dFdXpD[(i + 1) * ncol1] +=
                -transIJ * rhoWghtD * bvs.rhoP[dId_np_j] * dGamma;

            // dS
            for (USI k = 0; k < np; k++) {
                dFdXsB[(i + 1) * ncol2 + k] +=
                    transIJ * bvs.dPcdS[bId_np_j * np + k];
                dFdXsE[(i + 1) * ncol2 + k] -=
                    transIJ * bvs.dPcdS[eId_np_j * np + k];
                dFdXsU[(i + 1) * ncol2 + k] +=
                    Akd * bvs.dKrdS[uId_np_j * np + k] / mu * xi * xij * dP;
            }
            // dxij
            for (USI k = 0; k < nc; k++) {
                rhox = bvs.rhox[uId_np_j * nc + k];
                xix  = bvs.xix[uId_np_j * nc + k];
                mux  = bvs.mux[uId_np_j * nc + k];
                tmp  = -transIJ * rhoWghtU * rhox * dGamma;
                tmp += transJ * xix * xij * dP;
                tmp += -transIJ * mux / mu * dP;
                dFdXsU[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                dFdXsD[(i + 1) * ncol2 + np + j * nc + k] +=
                    -transIJ * rhoWghtD * bvs.rhox[dId_np_j * nc + k] * dGamma;
            }
            // k = i
            dFdXsU[(i + 1) * ncol2 + np + j * nc + i] += transJ * xi * dP;
        }
    }
}


void OCPConvection01::AssembleMatAIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk, FluxVarSet& fvs)
{
    const USI& ncol1  = fvs.ncol1;
    const USI& ncol2  = fvs.ncol2;
    auto&      dFdXpB = fvs.dFdXpB;
    auto&      dFdXpE = fvs.dFdXpE;
    auto&      dFdXsB = fvs.dFdXsB;
    auto&      dFdXsE = fvs.dFdXsE;

    const BulkVarSet& bvs = bk.vs;

    const OCP_USI bId    = bp.BId();
    const OCP_USI eId    = bp.EId();
    const OCP_DBL Akd    = bp.Trans();
    const OCP_DBL dGamma = GRAVITY_FACTOR * (bvs.depth[bId] - bvs.depth[eId]);

    OCP_USI  bId_np_j, eId_np_j, uId_np_j, dId_np_j;
    OCP_BOOL phaseExistDj;
    OCP_DBL  rhoWghtU, rhoWghtD;
    OCP_DBL  dP, transJ, transIJ, kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux, tmp;

    OCP_DBL* dFdXpU;     // up    bulk: dF / dXp
    OCP_DBL* dFdXpD;     // down  bulk: dF / dXp
    OCP_DBL* dFdXsU;     // up    bulk: dF / dXs
    OCP_DBL* dFdXsD;     // down  bulk: dF / dXs

    OCP_BOOL bIdFIM{ OCP_FALSE }, eIdFIM{ OCP_FALSE };
    if (bk.bulkTypeAIM.IfFIMbulk(bId)) bIdFIM = OCP_TRUE;
    if (bk.bulkTypeAIM.IfFIMbulk(eId)) eIdFIM = OCP_TRUE;

    if (!bIdFIM && !eIdFIM) {
        // both are explicit
        for (USI j = 0; j < np; j++) {
            uId_np_j = bcvs.upblock[c * np + j] * np + j;
            if (!bvs.phaseExist[uId_np_j]) continue;
            xi = bvs.xi[uId_np_j];
            kr = bvs.kr[uId_np_j];
            mu = bvs.mu[uId_np_j];
            transJ = Akd * xi * kr / mu;
            for (USI i = 0; i < nc; i++) {
                xij = bvs.xij[uId_np_j * nc + i];
                transIJ = xij * transJ;
                // Pressure -- Primary var
                dFdXpB[(i + 1) * ncol1] += transIJ;
                dFdXpE[(i + 1) * ncol1] -= transIJ;

                // maybe more derivatives should be considered  -- xiP, rhoP, muP
                // test
            }
        }
    }
    else if ((bIdFIM && !eIdFIM) || (!bIdFIM && eIdFIM)) {
        // one is explicit, one is implicit

        OCP_BOOL uIdFIM;
        OCP_DBL  flag_be;

        for (USI j = 0; j < np; j++) {
            uId_np_j = bcvs.upblock[c * np + j] * np + j;
            if (!bvs.phaseExist[uId_np_j]) continue;
            bId_np_j     = bId * np + j;
            eId_np_j     = eId * np + j;

            if (bId_np_j == uId_np_j) {
                dId_np_j     = eId_np_j;
                phaseExistDj = bvs.phaseExist[dId_np_j];                
                uIdFIM       = bIdFIM;               
                flag_be      = 1;
                dFdXpU       = &dFdXpB[0];
                dFdXpD       = &dFdXpE[0];
                dFdXsU       = &dFdXsB[0];
                dFdXsD       = &dFdXsE[0];  
            }
            else {
                dId_np_j     = bId_np_j;
                phaseExistDj = bvs.phaseExist[dId_np_j];               
                uIdFIM       = eIdFIM;               
                flag_be      = -1;
                dFdXpU       = &dFdXpE[0];
                dFdXpD       = &dFdXpB[0];
                dFdXsU       = &dFdXsE[0];
                dFdXsD       = &dFdXsB[0];
            }
            if (phaseExistDj) {
                rhoWghtU = 0.5;
                rhoWghtD = 0.5;
            }
            else {
                rhoWghtU = 1;
                rhoWghtD = 0;
            }

            dP     = bcvs.dP[c * np + j];
            xi     = bvs.xi[uId_np_j];
            kr     = bvs.kr[uId_np_j];
            mu     = bvs.mu[uId_np_j];
            transJ = Akd * kr / mu;            

            if (uIdFIM) {
                // upblock is implicit
                muP  = bvs.muP[uId_np_j];
                xiP  = bvs.xiP[uId_np_j];
                rhoP = bvs.rhoP[uId_np_j];
                for (USI i = 0; i < nc; i++) {
                    xij     = bvs.xij[uId_np_j * nc + i];
                    transIJ = xij * xi * transJ;

                    // Pressure -- Primary var
                    dFdXpB[(i + 1) * ncol1] += transIJ;
                    dFdXpE[(i + 1) * ncol1] -= transIJ;

                    tmp = transIJ * (-rhoP * dGamma) * rhoWghtU;
                    tmp += xij * transJ * xiP * dP;
                    tmp += -transIJ * muP / mu * dP;
                    dFdXpU[(i + 1) * ncol1] += tmp;

                    // maybe more derivatives should be considered  -- xiP, rhoP, muP

                    // Second var
                    // dS
                    for (USI k = 0; k < np; k++) {
                        tmp = transIJ * bvs.dPcdS[uId_np_j * np + k] * flag_be;
                        tmp += Akd * xij * xi / mu * bvs.dKrdS[uId_np_j * np + k] * dP;
						dFdXsU[(i + 1) * ncol2 + k] += tmp;
                    }
                    // d xij
                    for (USI k = 0; k < nc; k++) {
                        rhox = bvs.rhox[uId_np_j * nc + k] * rhoWghtU;
                        xix  = bvs.xix[uId_np_j * nc + k];
                        mux  = bvs.mux[uId_np_j * nc + k];
                        tmp  = -transIJ * rhox * dGamma;
                        tmp  += xij * transJ * xix * dP;
                        tmp  += -transIJ * mux / mu * dP;
                        dFdXsU[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                    }
                    dFdXsU[(i + 1) * ncol2 + np + j * nc + i] += transJ * xi * dP;
                }
            }
            else {
                // downblock is implicit
                rhoP = bvs.rhoP[dId_np_j];

                for (USI i = 0; i < nc; i++) {
                    xij     = bvs.xij[uId_np_j * nc + i];
                    transIJ = xij * xi * transJ;

                    // Pressure -- Primary var
                    dFdXpB[(i + 1) * ncol1] += transIJ;
                    dFdXpE[(i + 1) * ncol1] -= transIJ;

                    dFdXpD[(i + 1) * ncol1] -= transIJ * rhoP * dGamma * rhoWghtD;

                    // maybe more derivatives should be considered  -- xiP, rhoP, muP

                    // Second var
                    // dS
                    for (USI k = 0; k < np; k++) {
                        dFdXsD[(i + 1) * ncol2 + k] +=
                            transIJ * bvs.dPcdS[dId_np_j * np + k] * (-flag_be);
                    }
                    // d xij
                    for (USI k = 0; k < nc; k++) {
                        rhox = bvs.rhox[dId_np_j * nc + k] * rhoWghtD;
                        tmp  = -transIJ * rhox * dGamma;
                        dFdXsD[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                    }
                }               
            }
        }
    }
    else {

        // both are implicit
        for (USI j = 0; j < np; j++) {
            uId_np_j     = bcvs.upblock[c * np + j] * np + j;
            if (!bvs.phaseExist[uId_np_j]) continue;
            bId_np_j     = bId * np + j;
            eId_np_j     = eId * np + j;

            if (bId_np_j == uId_np_j) {
                dId_np_j     = eId_np_j;
                phaseExistDj = bvs.phaseExist[dId_np_j];
                dFdXpU       = &dFdXpB[0];
                dFdXpD       = &dFdXpE[0];
                dFdXsU       = &dFdXsB[0];
                dFdXsD       = &dFdXsE[0];
            }
            else {
                dId_np_j     = bId_np_j;
                phaseExistDj = bvs.phaseExist[dId_np_j];
                dFdXpU       = &dFdXpE[0];
                dFdXpD       = &dFdXpB[0];
                dFdXsU       = &dFdXsE[0];
                dFdXsD       = &dFdXsB[0];
            }
            if (phaseExistDj) {
                rhoWghtU = 0.5;
                rhoWghtD = 0.5;
            }
            else {
                rhoWghtU = 1;
                rhoWghtD = 0;
            }

            dP     = bcvs.dP[c * np + j];
            xi     = bvs.xi[uId_np_j];
            kr     = bvs.kr[uId_np_j];
            mu     = bvs.mu[uId_np_j];
            muP    = bvs.muP[uId_np_j];
            xiP    = bvs.xiP[uId_np_j];
            transJ = Akd * kr / mu;

            for (USI i = 0; i < nc; i++) {
                xij     = bvs.xij[uId_np_j * nc + i];
                transIJ = xij * xi * transJ;

                // dP
                dFdXpB[(i + 1) * ncol1] += transIJ;
                dFdXpE[(i + 1) * ncol1] -= transIJ;

                tmp = transJ * xiP * xij * dP;
                tmp += -transIJ * muP / mu * dP;
                dFdXpU[(i + 1) * ncol1] +=
                    (tmp - transIJ * rhoWghtU * bvs.rhoP[uId_np_j] * dGamma);
                dFdXpD[(i + 1) * ncol1] +=
                    -transIJ * rhoWghtD * bvs.rhoP[dId_np_j] * dGamma;

                // dS
                for (USI k = 0; k < np; k++) {
                    dFdXsB[(i + 1) * ncol2 + k] +=
                        transIJ * bvs.dPcdS[bId_np_j * np + k];
                    dFdXsE[(i + 1) * ncol2 + k] -=
                        transIJ * bvs.dPcdS[eId_np_j * np + k];
                    dFdXsU[(i + 1) * ncol2 + k] +=
                        Akd * bvs.dKrdS[uId_np_j * np + k] / mu * xi * xij * dP;
                }
                // dxij
                for (USI k = 0; k < nc; k++) {
                    rhox = bvs.rhox[uId_np_j * nc + k];
                    xix  = bvs.xix[uId_np_j * nc + k];
                    mux  = bvs.mux[uId_np_j * nc + k];
                    tmp  = -transIJ * rhoWghtU * rhox * dGamma;
                    tmp  += transJ * xix * xij * dP;
                    tmp  += -transIJ * mux / mu * dP;
                    dFdXsU[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                    dFdXsD[(i + 1) * ncol2 + np + j * nc + k] +=
                        -transIJ * rhoWghtD * bvs.rhox[dId_np_j * nc + k] * dGamma;
                }
                dFdXsU[(i + 1) * ncol2 + np + j * nc + i] += transJ * xi * dP;
            }
        }
    }
}


void OCPConvection01::AssembleMatIMPEC(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk, FluxVarSet& fvs)
{

    auto& valbb = fvs.valbb;   
    auto& valee = fvs.valee;
    auto& rhsb  = fvs.rhsb;   
    auto& rhse  = fvs.rhse;

    const BulkVarSet& bvs = bk.vs;
    
    const OCP_USI bId = bp.BId();
    const OCP_USI eId = bp.EId();
    const OCP_DBL Akd = bp.Trans();   
    const OCP_DBL dP  = bvs.P[bId] - bvs.P[eId];
     

    for (USI j = 0; j < np; j++) {
        const OCP_USI uId_np_j = bcvs.upblock[c * np + j] * np + j;        
        if (!bvs.phaseExist[uId_np_j]) continue;

        OCP_DBL valbi = 0;
        OCP_DBL valei = 0;

        for (USI i = 0; i < nc; i++) {
            valbi += bvs.vfi[bId * nc + i] * bvs.xij[uId_np_j * nc + i];
            valei += bvs.vfi[eId * nc + i] * bvs.xij[uId_np_j * nc + i];
        }

        OCP_DBL tmp = bvs.xi[uId_np_j] * Akd * bvs.kr[uId_np_j] / bvs.mu[uId_np_j];
        valbb += tmp * valbi;
        valee += tmp * valei;
        tmp   *= -(bcvs.dP[c * np + j] - dP);
        rhsb  += tmp * valbi;
        rhse  -= tmp * valei;
    }
}


////////////////////////////////////////////
// OCPConvection02
////////////////////////////////////////////


void OCPConvection02::CalFlux(const BulkConnPair& bp, const Bulk& bk, FluxVarSet& fvs)
{
    // Calculte upblock, rho, vj, flux_ni

    const BulkVarSet& bvs = bk.vs;

    auto& flux_ni = fvs.flux_ni;

    const OCP_USI bId = bp.BId();
    const OCP_USI eId = bp.EId();
    const OCP_DBL Akd = bp.Trans();
    OCP_USI       bId_np_j, eId_np_j, uId_np_j;
    OCP_BOOL      exbegin, exend;

    for (USI j = 0; j < np; j++) {
        bId_np_j = bId * np + j;
        eId_np_j = eId * np + j;

        exbegin = bvs.phaseExist[bId_np_j];
        exend = bvs.phaseExist[eId_np_j];

        if ((!exbegin) && (!exend)) {
            upblock[j] = bId;
            dP[j] = 0;
            vj[j] = 0;
            continue;
        }

        upblock[j] = bId;

        if (j == 0) {
            dP[j] = (bvs.Pj[bId_np_j] - bvs.Pj[eId_np_j])
                - bvs.dzMtrx[bId] * (bvs.S[bId_np_j + 1] - bvs.S[eId_np_j + 1])
                * (bvs.rho[bId_np_j + 1] - bvs.rho[bId_np_j]) * GRAVITY_FACTOR * 0.5;
        }
        else if (j == 1) {
            dP[j] = (bvs.Pj[bId_np_j] - bvs.Pj[eId_np_j])
                + bvs.dzMtrx[bId] * (bvs.S[bId_np_j] - bvs.S[eId_np_j])
                * (bvs.rho[bId_np_j] - bvs.rho[bId_np_j - 1]) * GRAVITY_FACTOR * 0.5;
        }
        else {
            OCP_ABORT("Inavailable!");
        }

        if (dP[j] < 0)  upblock[j] = eId;

        uId_np_j = upblock[j] * np + j;

        if (!bvs.phaseExist[uId_np_j]) {
            vj[j] = 0;
            continue;
        }

        vj[j] = Akd * bvs.kr[uId_np_j] / bvs.mu[uId_np_j] * dP[j];

        for (USI i = 0; i < nc; i++) {
            flux_ni[i] += vj[j] * bvs.xi[uId_np_j] * bvs.xij[uId_np_j * nc + i];
        }
    }
}

void OCPConvection02::AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk, FluxVarSet& fvs)
{
    const USI& ncol1  = fvs.ncol1;
    const USI& ncol2  = fvs.ncol2;
    auto&      dFdXpB = fvs.dFdXpB;
    auto&      dFdXpE = fvs.dFdXpE;
    auto&      dFdXsB = fvs.dFdXsB;
    auto&      dFdXsE = fvs.dFdXsE;

    const BulkVarSet& bvs = bk.vs;

    const OCP_USI bId = bp.BId();
    const OCP_USI eId = bp.EId();
    const OCP_DBL Akd = bp.Trans();
    const OCP_DBL dGamma = GRAVITY_FACTOR * (bvs.depth[bId] - bvs.depth[eId]);

    OCP_USI  bId_np_j, eId_np_j, uId_np_j, dId_np_j;
    OCP_BOOL phaseExistDj;
    OCP_DBL  rhoWghtU, rhoWghtD;
    OCP_DBL  dP, transJ, transIJ, kr, mu, xi, xij, xiP, muP, rhox, xix, mux, tmp;

    OCP_DBL* dFdXpU;     // up    bulk: dF / dXp
    OCP_DBL* dFdXpD;     // down  bulk: dF / dXp
    OCP_DBL* dFdXsU;     // up    bulk: dF / dXs
    OCP_DBL* dFdXsD;     // down  bulk: dF / dXs

    for (USI j = 0; j < np; j++) {
        uId_np_j = bcvs.upblock[c * np + j] * np + j;
        if (!bvs.phaseExist[uId_np_j]) continue;
        bId_np_j = bId * np + j;
        eId_np_j = eId * np + j;

        if (bId_np_j == uId_np_j) {
            dId_np_j = eId_np_j;
            phaseExistDj = bvs.phaseExist[dId_np_j];
            dFdXpU = &dFdXpB[0];
            dFdXpD = &dFdXpE[0];
            dFdXsU = &dFdXsB[0];
            dFdXsD = &dFdXsE[0];
        }
        else {
            dId_np_j = bId_np_j;
            phaseExistDj = bvs.phaseExist[dId_np_j];

            dFdXpU = &dFdXpE[0];
            dFdXpD = &dFdXpB[0];
            dFdXsU = &dFdXsE[0];
            dFdXsD = &dFdXsB[0];
        }
        if (phaseExistDj) {
            rhoWghtU = 0.5;
            rhoWghtD = 0.5;
        }
        else {
            rhoWghtU = 1;
            rhoWghtD = 0;
        }

        dP = bcvs.dP[c * np + j];

        xi = bvs.xi[uId_np_j];
        kr = bvs.kr[uId_np_j];
        mu = bvs.mu[uId_np_j];
        muP = bvs.muP[uId_np_j];
        xiP = bvs.xiP[uId_np_j];
        transJ = Akd * kr / mu;

        for (USI i = 0; i < nc; i++) {
            xij = bvs.xij[uId_np_j * nc + i];
            transIJ = xij * xi * transJ;

            // dP
            dFdXpB[(i + 1) * ncol1] += transIJ;
            dFdXpE[(i + 1) * ncol1] -= transIJ;

            tmp = transJ * xiP * xij * dP;
            tmp += -transIJ * muP / mu * dP;
            dFdXpU[(i + 1) * ncol1] +=
                (tmp - transIJ * rhoWghtU * bvs.rhoP[uId_np_j] * dGamma);
            dFdXpD[(i + 1) * ncol1] +=
                -transIJ * rhoWghtD * bvs.rhoP[dId_np_j] * dGamma;

            // dS
            for (USI k = 0; k < np; k++) {
                dFdXsB[(i + 1) * ncol2 + k] +=
                    transIJ * bvs.dPcdS[bId_np_j * np + k];
                dFdXsE[(i + 1) * ncol2 + k] -=
                    transIJ * bvs.dPcdS[eId_np_j * np + k];
                dFdXsU[(i + 1) * ncol2 + k] +=
                    Akd * bvs.dKrdS[uId_np_j * np + k] / mu * xi * xij * dP;
            }
            // dxij
            for (USI k = 0; k < nc; k++) {
                rhox = bvs.rhox[uId_np_j * nc + k];
                xix = bvs.xix[uId_np_j * nc + k];
                mux = bvs.mux[uId_np_j * nc + k];
                tmp = -transIJ * rhoWghtU * rhox * dGamma;
                tmp += transJ * xix * xij * dP;
                tmp += -transIJ * mux / mu * dP;
                dFdXsU[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                dFdXsD[(i + 1) * ncol2 + np + j * nc + k] +=
                    -transIJ * rhoWghtD * bvs.rhox[dId_np_j * nc + k] * dGamma;
            }
            dFdXsU[(i + 1) * ncol2 + np + j * nc + i] += transJ * xi * dP;
        }
    }
}


////////////////////////////////////////////
// OCPConvectionT01
////////////////////////////////////////////


void OCPConvectionT01::CalFlux(const BulkConnPair& bp, const Bulk& bk, FluxVarSet& fvs)
{
    // Calculte upblock, rho, vj, flux_ni, conduct_H
    auto& flux_ni = fvs.flux_ni;

    const OCP_USI bId = bp.BId();
    const OCP_USI eId = bp.EId();

    const BulkVarSet& bvs = bk.vs; 

    if (bvs.cType[bId] == BulkContent::rf && bvs.cType[eId] == BulkContent::rf) {
        OCP_USI       bId_np_j, eId_np_j, uId_np_j;
        OCP_BOOL      exbegin, exend;
        OCP_DBL       rho;

        const OCP_DBL Akd = bp.Trans();

        for (USI j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            exbegin = bvs.phaseExist[bId_np_j];
            exend = bvs.phaseExist[eId_np_j];

            if ((exbegin) && (exend)) {
                rho = (bvs.rho[bId_np_j] + bvs.rho[eId_np_j]) / 2;
            }
            else if (exbegin && (!exend)) {
                rho = bvs.rho[bId_np_j];
            }
            else if ((!exbegin) && (exend)) {
                rho = bvs.rho[eId_np_j];
            }
            else {
                upblock[j] = bId;
                dP[j]      = 0;
                vj[j] = 0;
                Hj[j] = 0;
                continue;
            }

            upblock[j] = bId;
            dP[j]      = (bvs.Pj[bId_np_j] - GRAVITY_FACTOR * rho * bvs.depth[bId]) -
                         (bvs.Pj[eId_np_j] - GRAVITY_FACTOR * rho * bvs.depth[eId]);

            if (dP[j] < 0)  upblock[j] = eId;

            uId_np_j = upblock[j] * np + j;

            if (!bvs.phaseExist[uId_np_j]) {
                vj[j] = 0;
                Hj[j] = 0;
                continue;
            }

            vj[j] = Akd * bvs.kr[uId_np_j] / bvs.mu[uId_np_j] * dP[j];
            Hj[j] = vj[j] * bvs.xi[uId_np_j] * bvs.H[uId_np_j];

            for (USI i = 0; i < nc; i++) {
                flux_ni[i] += vj[j] * bvs.xi[uId_np_j] * bvs.xij[uId_np_j * nc + i];
            }
        }
    }
}


void OCPConvectionT01::AssembleMatFIM(const BulkConnPair& bp, const OCP_USI& c, const BulkConnVarSet& bcvs, const Bulk& bk, FluxVarSet& fvs)
{
    const USI& ncol1  = fvs.ncol1;
    const USI& ncol2  = fvs.ncol2;
    auto&      dFdXpB = fvs.dFdXpB;
    auto&      dFdXpE = fvs.dFdXpE;
    auto&      dFdXsB = fvs.dFdXsB;
    auto&      dFdXsE = fvs.dFdXsE;

    const OCP_USI bId = bp.BId();
    const OCP_USI eId = bp.EId();

    const BulkVarSet& bvs = bk.vs;

    if (bvs.cType[bId] == BulkContent::rf && bvs.cType[eId] == BulkContent::rf) {
        // Fluid Bulk Connection
        const OCP_DBL Akd = bp.Trans();
        const OCP_DBL dGamma = GRAVITY_FACTOR * (bvs.depth[bId] - bvs.depth[eId]);

        OCP_USI  bId_np_j, eId_np_j, uId_np_j, dId_np_j;
        OCP_BOOL phaseExistDj;
        OCP_DBL  xi, xij, kr, mu, rhox, xiP, xiT, xix, muP, muT, mux, H, HT, Hx;
        OCP_DBL  dP;
        OCP_DBL  transJ, transIJ, transH;
        OCP_DBL  rhoWghtU, rhoWghtD;
        OCP_DBL  tmp;

        OCP_DBL* dFdXpU;     // up    bulk: dF / dXp
        OCP_DBL* dFdXpD;     // down  bulk: dF / dXp
        OCP_DBL* dFdXsU;     // up    bulk: dF / dXs
        OCP_DBL* dFdXsD;     // down  bulk: dF / dXs

        for (USI j = 0; j < np; j++) {
            uId_np_j = bcvs.upblock[c * np + j] * np + j;
            if (!bvs.phaseExist[uId_np_j]) continue;
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            if (bId_np_j == uId_np_j) {
                dId_np_j     = eId_np_j;
                phaseExistDj = bvs.phaseExist[dId_np_j];
                dFdXpU       = &dFdXpB[0];
                dFdXpD       = &dFdXpE[0];
                dFdXsU       = &dFdXsB[0];
                dFdXsD       = &dFdXsE[0];
            }
            else {
                dId_np_j     = bId_np_j;
                phaseExistDj = bvs.phaseExist[dId_np_j];
                dFdXpU       = &dFdXpE[0];
                dFdXpD       = &dFdXpB[0];
                dFdXsU       = &dFdXsE[0];
                dFdXsD       = &dFdXsB[0];
            }
            if (phaseExistDj) {
                rhoWghtU = 0.5;
                rhoWghtD = 0.5;
            }
            else {
                rhoWghtU = 1;
                rhoWghtD = 0;
            }

            dP     = bcvs.dP[c * np + j];
            xi     = bvs.xi[uId_np_j];
            kr     = bvs.kr[uId_np_j];
            mu     = bvs.mu[uId_np_j];
            xiP    = bvs.xiP[uId_np_j];
            xiT    = bvs.xiT[uId_np_j];
            muP    = bvs.muP[uId_np_j];
            muT    = bvs.muT[uId_np_j];
            H      = bvs.H[uId_np_j];
            HT     = bvs.HT[uId_np_j];
            transJ = Akd * kr / mu;

            // Mass Conservation
            for (USI i = 0; i < nc; i++) {
                xij     = bvs.xij[uId_np_j * nc + i];
                transIJ = transJ * xi * xij;

                // dP
                dFdXpB[(i + 1) * ncol1] += transIJ;
                dFdXpE[(i + 1) * ncol1] -= transIJ;

                tmp = transJ * xiP * xij * dP;
                tmp += -transIJ * muP / mu * dP;
                dFdXpU[(i + 1) * ncol1] +=
                    (tmp - transIJ * rhoWghtU * bvs.rhoP[uId_np_j] * dGamma);
                dFdXpD[(i + 1) * ncol1] +=
                    -transIJ * rhoWghtD * bvs.rhoP[dId_np_j] * dGamma;

                // dT
                tmp = transJ * xiT * xij * dP;
                tmp += -transIJ * muT / mu * dP;
                dFdXpU[(i + 2) * ncol1 - 1] +=
                    (tmp - transIJ * rhoWghtU * bvs.rhoT[uId_np_j] * dGamma);
                dFdXpD[(i + 2) * ncol1 - 1] +=
                    -transIJ * rhoWghtD * bvs.rhoT[dId_np_j] * dGamma;

                // dS
                for (USI k = 0; k < np; k++) {
                    dFdXsB[(i + 1) * ncol2 + k] +=
                        transIJ * bvs.dPcdS[bId_np_j * np + k];
                    dFdXsE[(i + 1) * ncol2 + k] -=
                        transIJ * bvs.dPcdS[eId_np_j * np + k];
                    dFdXsU[(i + 1) * ncol2 + k] +=
                        Akd * bvs.dKrdS[uId_np_j * np + k] / mu * xi * xij * dP;
                }
                // dxij
                for (USI k = 0; k < nc; k++) {
                    rhox = bvs.rhox[uId_np_j * nc + k];
                    xix = bvs.xix[uId_np_j * nc + k];
                    mux = bvs.mux[uId_np_j * nc + k];
                    tmp = -transIJ * rhoWghtU * rhox * dGamma;
                    tmp += transJ * xix * xij * dP;
                    tmp += -transIJ * mux / mu * dP;
                    dFdXsU[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                    dFdXsD[(i + 1) * ncol2 + np + j * nc + k] +=
                        -transIJ * rhoWghtD * bvs.rhox[dId_np_j * nc + k] * dGamma;
                }
                dFdXsU[(i + 1) * ncol2 + np + j * nc + i] += transJ * xi * dP;
            }

            // Energy Conservation
            transH = transJ * xi * H;
            // dP
            dFdXpB[(ncol1 - 1) * ncol1] += transH;
            dFdXpE[(ncol1 - 1) * ncol1] -= transH;

            tmp = transJ * xiP * H * dP;
            tmp += -transJ * xi * muP / mu * dP * H;
            dFdXpU[(ncol1 - 1) * ncol1] +=
                (tmp - transH * rhoWghtU * bvs.rhoP[uId_np_j] * dGamma);
            dFdXpD[(ncol1 - 1) * ncol1] +=
                -transH * rhoWghtD * bvs.rhoP[dId_np_j] * dGamma;

            // dT
            tmp = transJ * xiT * H * dP;
            tmp += transJ * xi * HT * dP;
            tmp += -transH * muT / mu * dP;
            dFdXpU[ncol1 * ncol1 - 1] +=
                (tmp - transH * rhoWghtU * bvs.rhoT[uId_np_j] * dGamma);
            dFdXpD[ncol1 * ncol1 - 1] +=
                -transH * rhoWghtD * bvs.rhoT[dId_np_j] * dGamma;

            // dS
            for (USI k = 0; k < np; k++) {
                dFdXsB[(nc + 1) * ncol2 + k] +=
                    transH * bvs.dPcdS[bId_np_j * np + k];
                dFdXsE[(nc + 1) * ncol2 + k] -=
                    transH * bvs.dPcdS[eId_np_j * np + k];
                dFdXsU[(nc + 1) * ncol2 + k] +=
                    Akd * bvs.dKrdS[uId_np_j * np + k] / mu * xi * H * dP;
            }
            // dxij
            for (USI k = 0; k < nc; k++) {
                rhox = bvs.rhox[uId_np_j * nc + k];
                xix = bvs.xix[uId_np_j * nc + k];
                mux = bvs.mux[uId_np_j * nc + k];
                Hx = bvs.Hx[uId_np_j * nc + k];
                tmp = -transH * rhoWghtU * rhox * dGamma;
                tmp += transJ * xix * H * dP;
                tmp += transJ * xi * Hx * dP;
                tmp += -transH * mux / mu * dP;
                dFdXsU[(nc + 1) * ncol2 + np + j * nc + k] += tmp;
                dFdXsD[(nc + 1) * ncol2 + np + j * nc + k] +=
                    -transH * rhoWghtD * bvs.rhox[dId_np_j * nc + k] * dGamma;
            }
        }
    }
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           May/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/

