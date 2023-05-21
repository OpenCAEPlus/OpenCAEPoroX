/*! \file    OCPFlux.cpp
 *  \brief   OCPFlux class declaration
 *  \author  Shizhe Li
 *  \date    May/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


// OpenCAEPoroX header files
#include "OCPFlux.hpp"


////////////////////////////////////////////
// OCPFlux_IsoT
////////////////////////////////////////////


void OCPFlux_IsoT::CalFlux(const BulkPair& bp, const Bulk& bk)
{
	// Calculte upblock, rho, flux_vj, flux_ni

    fill(flux_ni.begin(), flux_ni.end(), 0.0);
    const USI& np     = numPhase;
    const USI& nc     = numCom;

    const OCP_USI bId = bp.BId();
    const OCP_USI eId = bp.EId();
    const OCP_DBL Akd = CONV1 * CONV2 * bp.Area(); 
    OCP_USI       bId_np_j, eId_np_j, uId_np_j;
    OCP_BOOL      exbegin, exend;
    OCP_DBL       dP;

    for (USI j = 0; j < np; j++) {
        bId_np_j = bId * np + j;
        eId_np_j = eId * np + j;

        exbegin  = bk.phaseExist[bId_np_j];
        exend    = bk.phaseExist[eId_np_j];

        if ((exbegin) && (exend)) {
            rho[j] = (bk.rho[bId_np_j] + bk.rho[eId_np_j]) / 2;
        }
        else if (exbegin && (!exend)) {
            rho[j] = bk.rho[bId_np_j];
        }
        else if ((!exbegin) && (exend)) {
            rho[j] = bk.rho[eId_np_j];
        }
        else {
            upblock[j] = bId;
            rho[j]     = 0;
            continue;
        }

        upblock[j] = bId;
        dP = (bk.Pj[bId_np_j] - GRAVITY_FACTOR * rho[j] * bk.depth[bId]) -
            (bk.Pj[eId_np_j] - GRAVITY_FACTOR * rho[j] * bk.depth[eId]);

        if (dP < 0)  upblock[j] = eId;

        uId_np_j = upblock[j] * np + j;

        if (!bk.phaseExist[uId_np_j]) continue;

        flux_vj[j] = Akd * bk.kr[uId_np_j] / bk.mu[uId_np_j] * dP;
       
        for (USI i = 0; i < nc; i++) {
            flux_ni[i] += flux_vj[j] * bk.xi[uId_np_j] * bk.xij[uId_np_j * nc + i];
        }
    }
}



void OCPFlux_IsoT::AssembleMatFIM(const BulkPair& bp, const OCP_USI& c, const BulkConnVal& bcv, const Bulk& bk)
{
    const USI& np     = numPhase;
    const USI& nc     = numCom;
    const USI  ncol   = nc + 1;
    const USI  ncol2  = np * nc + np;

    fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
    fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
    fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
    fill(dFdXsE.begin(), dFdXsE.end(), 0.0);

    const OCP_USI bId    = bp.BId();
    const OCP_USI eId    = bp.EId();
    const OCP_DBL Akd    = CONV1 * CONV2 * bp.Area();
    const OCP_DBL dGamma = GRAVITY_FACTOR * (bk.depth[bId] - bk.depth[eId]);

    OCP_USI  bId_np_j, eId_np_j, uId_np_j, dId_np_j;
    OCP_BOOL phaseExistDj;
    OCP_DBL  rhoWghtU, rhoWghtD;
    OCP_DBL  dP, transJ, transIJ, kr, mu, xi, xij, xiP, muP, rhox, xix, mux, tmp;

    OCP_DBL* dFdXpU;     // up    bulk: dF / dXp
    OCP_DBL* dFdXpD;     // down  bulk: dF / dXp
    OCP_DBL* dFdXsU;     // up    bulk: dF / dXs
    OCP_DBL* dFdXsD;     // down  bulk: dF / dXs

    for (USI j = 0; j < np; j++) {
        uId_np_j = bcv.upblock[c * np + j] * np + j;
        if (!bk.phaseExist[uId_np_j]) continue;
        bId_np_j     = bId * np + j;
        eId_np_j     = eId * np + j;

        if (bId_np_j == uId_np_j) {
            dId_np_j     = eId_np_j;
            phaseExistDj = bk.phaseExist[dId_np_j];                      
            dFdXpU       = &dFdXpB[0];
            dFdXpD       = &dFdXpE[0];
            dFdXsU       = &dFdXsB[0];
            dFdXsD       = &dFdXsE[0];            
        }
        else {
            dId_np_j     = bId_np_j;
            phaseExistDj = bk.phaseExist[dId_np_j];
                      
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

        dP      = bk.Pj[bId_np_j] - bk.Pj[eId_np_j] - bcv.rho[c * np + j] * dGamma;
        xi      = bk.xi[uId_np_j];
        kr      = bk.kr[uId_np_j];
        mu      = bk.mu[uId_np_j];
        muP     = bk.muP[uId_np_j];
        xiP     = bk.xiP[uId_np_j];
        transJ  = Akd * kr / mu;

        for (USI i = 0; i < nc; i++) {
            xij    = bk.xij[uId_np_j * nc + i];
            transIJ = xij * xi * transJ;

            // dP
            dFdXpB[(i + 1) * ncol] += transIJ;
            dFdXpE[(i + 1) * ncol] -= transIJ;

            tmp = transJ * xiP * xij * dP;
            tmp += -transIJ * muP / mu * dP;
            dFdXpU[(i + 1) * ncol] +=
                (tmp - transIJ * rhoWghtU * bk.rhoP[uId_np_j] * dGamma);
            dFdXpD[(i + 1) * ncol] +=
                -transIJ * rhoWghtD * bk.rhoP[dId_np_j] * dGamma;

            // dS
            for (USI k = 0; k < np; k++) {
                dFdXsB[(i + 1) * ncol2 + k] +=
                    transIJ * bk.dPcdS[bId_np_j * np + k];
                dFdXsE[(i + 1) * ncol2 + k] -=
                    transIJ * bk.dPcdS[eId_np_j * np + k];
                dFdXsU[(i + 1) * ncol2 + k] +=
                    Akd * bk.dKrdS[uId_np_j * np + k] / mu * xi * xij * dP;
            }
            // dxij
            for (USI k = 0; k < nc; k++) {
                rhox = bk.rhox[uId_np_j * nc + k];
                xix  = bk.xix[uId_np_j * nc + k];
                mux  = bk.mux[uId_np_j * nc + k];
                tmp  = -transIJ * rhoWghtU * rhox * dGamma;
                tmp += transJ * xix * xij * dP;
                tmp += -transIJ * mux / mu * dP;
                dFdXsU[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                dFdXsD[(i + 1) * ncol2 + np + j * nc + k] +=
                    -transIJ * rhoWghtD * bk.rhox[dId_np_j * nc + k] * dGamma;
            }
            dFdXsU[(i + 1) * ncol2 + np + j * nc + i] += transJ * xi * dP;
        }
    }
}


void OCPFlux_IsoT::AssembleMatAIM(const BulkPair& bp, const OCP_USI& c, const BulkConnVal& bcv, const Bulk& bk)
{

    const USI& np     = numPhase;
    const USI& nc     = numCom;
    const USI  ncol   = nc + 1;
    const USI  ncol2  = np * nc + np;

    const OCP_USI& bId    = bp.BId();
    const OCP_USI& eId    = bp.EId();
    const OCP_DBL  Akd    = CONV1 * CONV2 * bp.Area();
    const OCP_DBL  dGamma = GRAVITY_FACTOR * (bk.depth[bId] - bk.depth[eId]);

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

        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);

        for (USI j = 0; j < np; j++) {
            uId_np_j = bcv.upblock[c * np + j] * np + j;
            if (!bk.phaseExist[uId_np_j]) continue;
            xi = bk.xi[uId_np_j];
            kr = bk.kr[uId_np_j];
            mu = bk.mu[uId_np_j];
            transJ = Akd * xi * kr / mu;
            for (USI i = 0; i < nc; i++) {
                xij = bk.xij[uId_np_j * nc + i];
                transIJ = xij * transJ;
                // Pressure -- Primary var
                dFdXpB[(i + 1) * ncol] += transIJ;
                dFdXpE[(i + 1) * ncol] -= transIJ;

                // maybe more derivatives should be considered  -- xiP, rhoP, muP
            }
        }
    }
    else if ((bIdFIM && !eIdFIM) || (!bIdFIM && eIdFIM)) {
        // one is explicit, one is implicit

        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
        fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
        fill(dFdXsE.begin(), dFdXsE.end(), 0.0);

        OCP_BOOL uIdFIM;
        OCP_DBL  flag_be;

        for (USI j = 0; j < np; j++) {
            uId_np_j = bcv.upblock[c * np + j] * np + j;
            if (!bk.phaseExist[uId_np_j]) continue;
            bId_np_j     = bId * np + j;
            eId_np_j     = eId * np + j;

            if (bId_np_j == uId_np_j) {
                dId_np_j     = eId_np_j;
                phaseExistDj = bk.phaseExist[dId_np_j];                
                uIdFIM       = bIdFIM;               
                flag_be      = 1;
                dFdXpU       = &dFdXpB[0];
                dFdXpD       = &dFdXpE[0];
                dFdXsU       = &dFdXsB[0];
                dFdXsD       = &dFdXsE[0];  
            }
            else {
                dId_np_j     = bId_np_j;
                phaseExistDj = bk.phaseExist[dId_np_j];               
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

            dP     = bk.Pj[bId_np_j] - bk.Pj[eId_np_j] - bcv.rho[c * np + j] * dGamma;
            xi     = bk.xi[uId_np_j];
            kr     = bk.kr[uId_np_j];
            mu     = bk.mu[uId_np_j];
            transJ = Akd * kr / mu;            

            if (uIdFIM) {
                // upblock is implicit
                muP  = bk.muP[uId_np_j];
                xiP  = bk.xiP[uId_np_j];
                rhoP = bk.rhoP[uId_np_j];
                for (USI i = 0; i < nc; i++) {
                    xij     = bk.xij[uId_np_j * nc + i];
                    transIJ = xij * xi * transJ;

                    // Pressure -- Primary var
                    dFdXpB[(i + 1) * ncol] += transIJ;
                    dFdXpE[(i + 1) * ncol] -= transIJ;

                    tmp = transIJ * (-rhoP * dGamma) * rhoWghtU;
                    tmp += xij * transJ * xiP * dP;
                    tmp += -transIJ * muP / mu * dP;
                    dFdXpU[(i + 1) * ncol] += tmp;

                    // maybe more derivatives should be considered  -- xiP, rhoP, muP

                    // Second var
                    // dS
                    for (USI k = 0; k < np; k++) {
                        tmp = transIJ * bk.dPcdS[uId_np_j * np + k] * flag_be;
                        tmp += Akd * xij * xi / mu * bk.dKrdS[uId_np_j * np + k] * dP;
						dFdXsU[(i + 1) * ncol2 + k] += tmp;
                    }
                    // d xij
                    for (USI k = 0; k < nc; k++) {
                        rhox = bk.rhox[uId_np_j * nc + k] * rhoWghtU;
                        xix  = bk.xix[uId_np_j * nc + k];
                        mux  = bk.mux[uId_np_j * nc + k];
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
                rhoP = bk.rhoP[dId_np_j];

                for (USI i = 0; i < nc; i++) {
                    xij     = bk.xij[uId_np_j * nc + i];
                    transIJ = xij * xi * transJ;

                    // Pressure -- Primary var
                    dFdXpB[(i + 1) * ncol] += transIJ;
                    dFdXpE[(i + 1) * ncol] -= transIJ;

                    dFdXpD[(i + 1) * ncol] -= transIJ * rhoP * dGamma * rhoWghtD;

                    // maybe more derivatives should be considered  -- xiP, rhoP, muP

                    // Second var
                    // dS
                    for (USI k = 0; k < np; k++) {
                        dFdXsD[(i + 1) * ncol2 + k] +=
                            transIJ * bk.dPcdS[dId_np_j * np + k] * (-flag_be);
                    }
                    // d xij
                    for (USI k = 0; k < nc; k++) {
                        rhox = bk.rhox[dId_np_j * nc + k] * rhoWghtD;
                        tmp  = -transIJ * rhox * dGamma;
                        dFdXsD[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                    }
                }               
            }
        }
    }
    else {

        // both are implicit
        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
        fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
        fill(dFdXsE.begin(), dFdXsE.end(), 0.0);

        for (USI j = 0; j < np; j++) {
            uId_np_j     = bcv.upblock[c * np + j] * np + j;
            if (!bk.phaseExist[uId_np_j]) continue;
            bId_np_j     = bId * np + j;
            eId_np_j     = eId * np + j;

            if (bId_np_j == uId_np_j) {
                dId_np_j     = eId_np_j;
                phaseExistDj = bk.phaseExist[dId_np_j];
                dFdXpU       = &dFdXpB[0];
                dFdXpD       = &dFdXpE[0];
                dFdXsU       = &dFdXsB[0];
                dFdXsD       = &dFdXsE[0];
            }
            else {
                dId_np_j     = bId_np_j;
                phaseExistDj = bk.phaseExist[dId_np_j];
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

            dP     = bk.Pj[bId_np_j] - bk.Pj[eId_np_j] - bcv.rho[c * np + j] * dGamma;
            xi     = bk.xi[uId_np_j];
            kr     = bk.kr[uId_np_j];
            mu     = bk.mu[uId_np_j];
            muP    = bk.muP[uId_np_j];
            xiP    = bk.xiP[uId_np_j];
            transJ = Akd * kr / mu;

            for (USI i = 0; i < nc; i++) {
                xij     = bk.xij[uId_np_j * nc + i];
                transIJ = xij * xi * transJ;

                // dP
                dFdXpB[(i + 1) * ncol] += transIJ;
                dFdXpE[(i + 1) * ncol] -= transIJ;

                tmp = transJ * xiP * xij * dP;
                tmp += -transIJ * muP / mu * dP;
                dFdXpU[(i + 1) * ncol] +=
                    (tmp - transIJ * rhoWghtU * bk.rhoP[uId_np_j] * dGamma);
                dFdXpD[(i + 1) * ncol] +=
                    -transIJ * rhoWghtD * bk.rhoP[dId_np_j] * dGamma;

                // dS
                for (USI k = 0; k < np; k++) {
                    dFdXsB[(i + 1) * ncol2 + k] +=
                        transIJ * bk.dPcdS[bId_np_j * np + k];
                    dFdXsE[(i + 1) * ncol2 + k] -=
                        transIJ * bk.dPcdS[eId_np_j * np + k];
                    dFdXsU[(i + 1) * ncol2 + k] +=
                        Akd * bk.dKrdS[uId_np_j * np + k] / mu * xi * xij * dP;
                }
                // dxij
                for (USI k = 0; k < nc; k++) {
                    rhox = bk.rhox[uId_np_j * nc + k];
                    xix  = bk.xix[uId_np_j * nc + k];
                    mux  = bk.mux[uId_np_j * nc + k];
                    tmp  = -transIJ * rhoWghtU * rhox * dGamma;
                    tmp  += transJ * xix * xij * dP;
                    tmp  += -transIJ * mux / mu * dP;
                    dFdXsU[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                    dFdXsD[(i + 1) * ncol2 + np + j * nc + k] +=
                        -transIJ * rhoWghtD * bk.rhox[dId_np_j * nc + k] * dGamma;
                }
                dFdXsU[(i + 1) * ncol2 + np + j * nc + i] += transJ * xi * dP;
            }
        }
    }
}


void OCPFlux_IsoT::AssembleMatIMPEC(const BulkPair& bp, const OCP_USI& c, const BulkConnVal& bcv, const Bulk& bk)
{

    const USI&    np  = numPhase;
    const USI&    nc  = numCom;
    
    const OCP_USI bId = bp.BId();
    const OCP_USI eId = bp.EId();
    const OCP_DBL Akd = CONV1 * CONV2 * bp.Area();
    const OCP_DBL dD  = GRAVITY_FACTOR * (bk.depth[bId] - bk.depth[eId]);    
     
    valbb = 0;   valee = 0;
    rhsb  = 0;   rhse  = 0;

    for (USI j = 0; j < np; j++) {
        const OCP_USI uId_np_j = bcv.upblock[c * np + j] * np + j;        
        if (!bk.phaseExist[uId_np_j]) continue;

        const OCP_DBL rho = bcv.rho[c * np + j];

        OCP_DBL valbi = 0;
        OCP_DBL valei = 0;

        for (USI i = 0; i < nc; i++) {
            valbi += bk.vfi[bId * nc + i] * bk.xij[uId_np_j * nc + i];
            valei += bk.vfi[eId * nc + i] * bk.xij[uId_np_j * nc + i];
        }

        OCP_DBL tmp = bk.xi[uId_np_j] * Akd * bk.kr[uId_np_j] / bk.mu[uId_np_j];
        valbb += tmp * valbi;
        valee += tmp * valei;
        tmp   *= rho * dD - (bk.Pc[bId * np + j] - bk.Pc[eId * np + j]);
        rhsb  += tmp * valbi;
        rhse  -= tmp * valei;
    }
}


////////////////////////////////////////////
// OCPFlux_T
////////////////////////////////////////////


void OCPFlux_T::CalFlux(const BulkPair& bp, const Bulk& bk)
{
    // Calculte upblock, rho, flux_vj, flux_ni, Adkt

    fill(flux_ni.begin(), flux_ni.end(), 0.0);

    const USI&    np  = numPhase;
    const USI&    nc  = numCom;

    const OCP_USI bId = bp.BId();
    const OCP_USI eId = bp.EId();
    const OCP_DBL T1  = bk.kt[bId] * bp.AreaB();
    const OCP_DBL T2  = bk.kt[eId] * bp.AreaE();

    Adkt              = 1 / (1 / T1 + 1 / T2);

    if (bp.Type() == 0) {
        OCP_USI       bId_np_j, eId_np_j, uId_np_j;
        OCP_BOOL      exbegin, exend;
        OCP_DBL       dP;

        const OCP_DBL Akd = CONV1 * CONV2 * bp.Area();

        for (USI j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            exbegin = bk.phaseExist[bId_np_j];
            exend = bk.phaseExist[eId_np_j];

            if ((exbegin) && (exend)) {
                rho[j] = (bk.rho[bId_np_j] + bk.rho[eId_np_j]) / 2;
            }
            else if (exbegin && (!exend)) {
                rho[j] = bk.rho[bId_np_j];
            }
            else if ((!exbegin) && (exend)) {
                rho[j] = bk.rho[eId_np_j];
            }
            else {
                upblock[j] = bId;
                rho[j] = 0;
                continue;
            }

            upblock[j] = bId;
            dP = (bk.Pj[bId_np_j] - GRAVITY_FACTOR * rho[j] * bk.depth[bId]) -
                (bk.Pj[eId_np_j] - GRAVITY_FACTOR * rho[j] * bk.depth[eId]);

            if (dP < 0)  upblock[j] = eId;

            uId_np_j = upblock[j] * np + j;

            if (!bk.phaseExist[uId_np_j]) continue;

            flux_vj[j] = Akd * bk.kr[uId_np_j] / bk.mu[uId_np_j] * dP;

            for (USI i = 0; i < nc; i++) {
                flux_ni[i] += flux_vj[j] * bk.xi[uId_np_j] * bk.xij[uId_np_j * nc + i];
            }
        }
    }
}


void OCPFlux_T::AssembleMatFIM(const BulkPair& bp, const OCP_USI& c, const BulkConnVal& bcv, const Bulk& bk)
{
    fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
    fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
    fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
    fill(dFdXsE.begin(), dFdXsE.end(), 0.0);

    const USI& np    = numPhase;
    const USI& nc    = numCom;
    const USI  ncol  = nc + 2;
    const USI  ncol2 = np * nc + np;

    const OCP_USI bId   = bp.BId();
    const OCP_USI eId   = bp.EId();
    const OCP_DBL areaB = bp.AreaB();
    const OCP_DBL areaE = bp.AreaE();
    const OCP_DBL T1    = bk.kt[bId] * areaB;
    const OCP_DBL T2    = bk.kt[eId] * areaE;
    Adkt                = 1 / (1 / T1 + 1 / T2);
                        
    const OCP_DBL tmpB  = pow(Adkt, 2) / pow(T1, 2) * areaB;
    const OCP_DBL tmpE  = pow(Adkt, 2) / pow(T2, 2) * areaE;

    const OCP_DBL dT    = bk.T[bId] - bk.T[eId];
    // Thermal Conduction always exist
    // dP
    dFdXpB[(ncol - 1) * ncol] += tmpB * bk.ktP[bId] * dT;
    dFdXpE[(ncol - 1) * ncol] += tmpE * bk.ktP[eId] * dT;
    // dT
    dFdXpB[ncol * ncol - 1]   += Adkt + tmpB * bk.ktT[bId] * dT;
    dFdXpE[ncol * ncol - 1]   += -Adkt + tmpE * bk.ktT[eId] * dT;
    // dS
    for (OCP_USI j = 0; j < np; j++) {
        dFdXsB[(nc + 1) * ncol2 + j] += tmpB * bk.ktS[bId * np + j] * dT;
        dFdXsE[(nc + 1) * ncol2 + j] += tmpE * bk.ktS[eId * np + j] * dT;
    }

    if (bp.Type() == 0) {
        // Fluid Bulk Connection
        const OCP_DBL Akd = CONV1 * CONV2 * bp.Area();
        const OCP_DBL dGamma = GRAVITY_FACTOR * (bk.depth[bId] - bk.depth[eId]);

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
            uId_np_j = bcv.upblock[c * np + j] * np + j;
            if (!bk.phaseExist[uId_np_j]) continue;
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            if (bId_np_j == uId_np_j) {
                dId_np_j     = eId_np_j;
                phaseExistDj = bk.phaseExist[dId_np_j];
                dFdXpU       = &dFdXpB[0];
                dFdXpD       = &dFdXpE[0];
                dFdXsU       = &dFdXsB[0];
                dFdXsD       = &dFdXsE[0];
            }
            else {
                dId_np_j     = bId_np_j;
                phaseExistDj = bk.phaseExist[dId_np_j];
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

            dP     = bk.Pj[bId_np_j] - bk.Pj[eId_np_j] - bcv.rho[c * np + j] * dGamma;
            xi     = bk.xi[uId_np_j];
            kr     = bk.kr[uId_np_j];
            mu     = bk.mu[uId_np_j];
            xiP    = bk.xiP[uId_np_j];
            xiT    = bk.xiT[uId_np_j];
            muP    = bk.muP[uId_np_j];
            muT    = bk.muT[uId_np_j];
            H      = bk.H[uId_np_j];
            HT     = bk.HT[uId_np_j];
            transJ = Akd * kr / mu;

            // Mass Conservation
            for (USI i = 0; i < nc; i++) {
                xij     = bk.xij[uId_np_j * nc + i];
                transIJ = transJ * xi * xij;

                // dP
                dFdXpB[(i + 1) * ncol] += transIJ;
                dFdXpE[(i + 1) * ncol] -= transIJ;

                tmp = transJ * xiP * xij * dP;
                tmp += -transIJ * muP / mu * dP;
                dFdXpU[(i + 1) * ncol] +=
                    (tmp - transIJ * rhoWghtU * bk.rhoP[uId_np_j] * dGamma);
                dFdXpD[(i + 1) * ncol] +=
                    -transIJ * rhoWghtD * bk.rhoP[dId_np_j] * dGamma;

                // dT
                tmp = transJ * xiT * xij * dP;
                tmp += -transIJ * muT / mu * dP;
                dFdXpU[(i + 2) * ncol - 1] +=
                    (tmp - transIJ * rhoWghtU * bk.rhoT[uId_np_j] * dGamma);
                dFdXpD[(i + 2) * ncol - 1] +=
                    -transIJ * rhoWghtD * bk.rhoT[dId_np_j] * dGamma;

                // dS
                for (USI k = 0; k < np; k++) {
                    dFdXsB[(i + 1) * ncol2 + k] +=
                        transIJ * bk.dPcdS[bId_np_j * np + k];
                    dFdXsE[(i + 1) * ncol2 + k] -=
                        transIJ * bk.dPcdS[eId_np_j * np + k];
                    dFdXsU[(i + 1) * ncol2 + k] +=
                        Akd * bk.dKrdS[uId_np_j * np + k] / mu * xi * xij * dP;
                }
                // dxij
                for (USI k = 0; k < nc; k++) {
                    rhox = bk.rhox[uId_np_j * nc + k];
                    xix = bk.xix[uId_np_j * nc + k];
                    mux = bk.mux[uId_np_j * nc + k];
                    tmp = -transIJ * rhoWghtU * rhox * dGamma;
                    tmp += transJ * xix * xij * dP;
                    tmp += -transIJ * mux / mu * dP;
                    dFdXsU[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                    dFdXsD[(i + 1) * ncol2 + np + j * nc + k] +=
                        -transIJ * rhoWghtD * bk.rhox[dId_np_j * nc + k] * dGamma;
                }
                dFdXsU[(i + 1) * ncol2 + np + j * nc + i] += transJ * xi * dP;
            }

            // Energy Conservation
            transH = transJ * xi * H;
            // dP
            dFdXpB[(ncol - 1) * ncol] += transH;
            dFdXpE[(ncol - 1) * ncol] -= transH;

            tmp = transJ * xiP * H * dP;
            tmp += -transJ * xi * muP / mu * dP * H;
            dFdXpU[(ncol - 1) * ncol] +=
                (tmp - transH * rhoWghtU * bk.rhoP[uId_np_j] * dGamma);
            dFdXpD[(ncol - 1) * ncol] +=
                -transH * rhoWghtD * bk.rhoP[dId_np_j] * dGamma;

            // dT
            tmp = transJ * xiT * H * dP;
            tmp += transJ * xi * HT * dP;
            tmp += -transH * muT / mu * dP;
            dFdXpU[ncol * ncol - 1] +=
                (tmp - transH * rhoWghtU * bk.rhoT[uId_np_j] * dGamma);
            dFdXpD[ncol * ncol - 1] +=
                -transH * rhoWghtD * bk.rhoT[dId_np_j] * dGamma;

            // dS
            for (USI k = 0; k < np; k++) {
                dFdXsB[(nc + 1) * ncol2 + k] +=
                    transH * bk.dPcdS[bId_np_j * np + k];
                dFdXsE[(nc + 1) * ncol2 + k] -=
                    transH * bk.dPcdS[eId_np_j * np + k];
                dFdXsU[(nc + 1) * ncol2 + k] +=
                    Akd * bk.dKrdS[uId_np_j * np + k] / mu * xi * H * dP;
            }
            // dxij
            for (USI k = 0; k < nc; k++) {
                rhox = bk.rhox[uId_np_j * nc + k];
                xix = bk.xix[uId_np_j * nc + k];
                mux = bk.mux[uId_np_j * nc + k];
                Hx = bk.Hx[uId_np_j * nc + k];
                tmp = -transH * rhoWghtU * rhox * dGamma;
                tmp += transJ * xix * H * dP;
                tmp += transJ * xi * Hx * dP;
                tmp += -transH * mux / mu * dP;
                dFdXsU[(nc + 1) * ncol2 + np + j * nc + k] += tmp;
                dFdXsD[(nc + 1) * ncol2 + np + j * nc + k] +=
                    -transH * rhoWghtD * bk.rhox[dId_np_j * nc + k] * dGamma;
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

