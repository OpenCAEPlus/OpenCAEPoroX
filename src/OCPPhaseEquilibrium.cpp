/*! \file    OCPPhaseEquilibrium.cpp
 *  \brief   OCPPhaseEquilibrium class declaration
 *  \author  Shizhe Li
 *  \date    Jul/28/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPPhaseEquilibrium.hpp"


void OCPPhaseEquilibrium::Setup(const ComponentParam& param, const USI& tarId, EoSCalculation* eosin)
{

    NC    = param.numCom;
    NPmax = param.numPhase;

    Cname = param.Cname;

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
    if (param.Acf.activity)  Acf = param.Acf.data[tarId];
    else                     OCP_ABORT("ACF hasn't been input!");

    eos = eosin;

    flashCtrl.SSMsta.maxIt = stoi(param.SSMparamSTA[0]);
    flashCtrl.SSMsta.tol   = stod(param.SSMparamSTA[1]);
    flashCtrl.SSMsta.eYt   = stod(param.SSMparamSTA[2]);
    flashCtrl.SSMsta.tol2  = flashCtrl.SSMsta.tol * flashCtrl.SSMsta.tol;

    flashCtrl.NRsta.maxIt  = stoi(param.NRparamSTA[0]);
    flashCtrl.NRsta.tol    = stod(param.NRparamSTA[1]);
    flashCtrl.NRsta.tol2   = flashCtrl.NRsta.tol * flashCtrl.NRsta.tol;

    flashCtrl.SSMsp.maxIt  = stoi(param.SSMparamSP[0]);
    flashCtrl.SSMsp.tol    = stod(param.SSMparamSP[1]);
    flashCtrl.SSMsp.tol2   = flashCtrl.SSMsp.tol * flashCtrl.SSMsp.tol;

    flashCtrl.NRsp.maxIt   = stoi(param.NRparamSP[0]);
    flashCtrl.NRsp.tol     = stod(param.NRparamSP[1]);
    flashCtrl.NRsp.tol2    = flashCtrl.NRsp.tol * flashCtrl.NRsp.tol;

    flashCtrl.RR.maxIt     = stoi(param.RRparam[0]);
    flashCtrl.RR.tol       = stod(param.RRparam[1]);
    flashCtrl.RR.tol2      = flashCtrl.RR.tol * flashCtrl.RR.tol;

    AllocateBasicVars();
    AllocateMethodVars();
}


void OCPPhaseEquilibrium::PhaseEquilibrium(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin,
    const USI& ftypein, const OCP_DBL& lNPin, const OCP_DBL* lKsin)
{

    SetInitalValue(Pin, Tin, Niin, ftypein, lNPin, lKsin);

    // Attention: sum of components' moles equals 1
    switch (ftype) {
    case 0:
        // flash from single phase
        flagSkip = OCP_TRUE;
        NP = 1;
        nu[0] = 1;
        x[0]  = zi;
        CalKwilson();
        while (!PhaseStable()) {
            NP++;
            PhaseSplit();
            if (NP == NPmax || NP == 1) break;
        }
        if (NP > 1) {
            flagSkip = OCP_FALSE;
        }

        break;
    case 1:
        // Skip Phase Stability analysis, only single phase exists
        flagSkip = OCP_TRUE;
        NP = 1;
        nu[0] = 1;
        x[0] = zi;
        break;

    case 2:
        // Skip Phase Stability analysis, two phases exist
        flagSkip = OCP_FALSE;
        NP = 2;
        Yt = 1.01;
        CalKwilson();
        PhaseSplit();
        break;

    default:
        OCP_ABORT("Wrong flash type!");
        break;
    }
}


void OCPPhaseEquilibrium::SetInitalValue(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* ziin,
    const USI& ftypein, const OCP_DBL& lNPin, const OCP_DBL* lKsin)
{
    P = Pin;
    T = Tin;
    copy(ziin, ziin + NC, zi.begin());
    ftype = ftypein;
    lNP   = lNPin;
    if (lNP == 2)
        copy(lKsin, lKsin + NC, lKs.begin());
}


///////////////////////////////////////////////
// Basic variables
///////////////////////////////////////////////


void OCPPhaseEquilibrium::AllocateBasicVars()
{
    // Allocate Memoery for Phase variables
    zi.resize(NC);
    nu.resize(NPmax);
    MW.resize(NPmax);
    x.resize(NPmax);
    n.resize(NPmax);
    phi.resize(NPmax);
    fug.resize(NPmax);
    for (USI j = 0; j < NPmax; j++) {
        x[j].resize(NC);
        n[j].resize(NC);
        phi[j].resize(NC);
        fug[j].resize(NC);
    }
    ln = n;
}

void OCPPhaseEquilibrium::x2n()
{
    // Total moles are 1
    for (USI j = 0; j < NP; j++) {
        for (USI i = 0; i < NC; i++) {
            n[j][i] = nu[j] * x[j][i];
        }
    }
}


void OCPPhaseEquilibrium::CalMW()
{
    // Calculate Molecular Weight of phase
    for (USI j = 0; j < NP; j++) {
        MW[j] = 0;
        for (USI i = 0; i < NC; i++) {
            MW[j] += x[j][i] * MWC[i];
        }
    }
}


USI OCPPhaseEquilibrium::FindMWmax()
{
    // find the phase id with the highest Molecular Weight
    USI     tmpId = 0;
    OCP_DBL tmpMW = MW[0];
    for (USI j = 1; j < NP; j++) {
        if (tmpMW < MW[j]) {
            tmpMW = MW[j];
            tmpId = j;
        }
    }
    return tmpId;
}


///////////////////////////////////////////////
// Method
///////////////////////////////////////////////


void OCPPhaseEquilibrium::AllocateMethodVars()
{
    Kw.resize(4);
    for (USI i = 0; i < 4; i++) {
        Kw[i].resize(NC);
    }
    Ks.resize(NPmax - 1);
    for (USI i = 0; i < NPmax - 1; i++) {
        Ks[i].resize(NC);
    }
    phiSta.resize(NC);
    fugSta.resize(NC);

    Y.resize(NC);
    di.resize(NC);
    resSTA.resize(NC);

    JmatSTA.resize(NC * NC);
    lKs.resize(NC);

    resRR.resize(NPmax - 1);
    resSP.resize(static_cast<size_t>(NC) * NPmax);
    JmatSP.resize(static_cast<size_t>(NC * NC) * NPmax * NPmax);
    lnfugX.resize(NPmax);
    lnfugN.resize(NPmax);
    for (USI j = 0; j < NPmax; j++) {
        lnfugX[j].resize(NC * NC);
        lnfugN[j].resize(NC * NC);
    }
    pivot.resize(static_cast<size_t>(NC + 1) * NPmax, 1);
    lenJmatWork = NC * (NPmax - 1);
    JmatWork.resize(lenJmatWork);

    pivot.resize(NPmax * static_cast<size_t>(NC), 1);
}

void OCPPhaseEquilibrium::CalKwilson()
{
    for (USI i = 0; i < NC; i++) {
        Kw[0][i] = (Pc[i] / P) * exp(5.373 * (1 + Acf[i]) * (1 - Tc[i] / T));
        Kw[1][i] = 1 / Kw[0][i];
        Kw[2][i] = pow(Kw[0][i], 1.0 / 3);
        Kw[3][i] = pow(Kw[1][i], 1.0 / 3);
    }
}


///////////////////////////////////////////////
// Phase-Stability Analysis
///////////////////////////////////////////////


OCP_BOOL OCPPhaseEquilibrium::PhaseStable()
{
    if (NP == 1) {
        testPId = 0;
    }
    else {
        CalMW();
        testPId = FindMWmax();
    }

    flashCtrl.SSMsta.conflag = OCP_FALSE;
    flashCtrl.SSMsta.curIt = 0;
    flashCtrl.NRsta.conflag = OCP_FALSE;
    flashCtrl.NRsta.curIt = 0;

    // Test if a phase is stable, if stable return OCP_TRUE, else return OCP_FALSE
    OCP_BOOL flag;
    USI      tmpNP = NP;

    if (lNP == 0) {
        // strict stability ananlysis
        flag = StableSSM(testPId);
    }
    else {
        flag = StableSSM01(testPId);
        if (!flag) tmpNP++;

        if (tmpNP != lNP) {
            flag = StableSSM(testPId);
            flagSkip = OCP_FALSE;
        }
    }
    itersSSMSTA += flashCtrl.SSMsta.curIt;
    itersNRSTA += flashCtrl.NRsta.curIt;
    countsSSMSTA++;
    countsNRSTA++;
    return flag;
}

OCP_BOOL OCPPhaseEquilibrium::StableSSM01(const USI& Id)
{
    // if unsatble, return OCP_FALSE
    // if stable, return OCP_TRUE

    OCP_DBL Stol = flashCtrl.SSMsta.tol2;
    USI     maxIt = flashCtrl.SSMsta.maxIt;
    OCP_DBL eYt = flashCtrl.SSMsta.eYt;
    OCP_DBL Ktol = flashCtrl.SSMsta.Ktol;
    OCP_DBL dYtol = flashCtrl.SSMsta.dYtol;
    // OCP_DBL& Sk = EoSctrl.SSMsta.curSk;
    OCP_DBL Se, Sk, dY;

    OCP_BOOL flag, Tsol; // Tsol, trivial solution
    USI      iter, k;

    const vector<OCP_DBL>& xj = x[Id];

    eos->CalFug(P, T, &x[Id][0], &fug[Id][0]);

    const vector<OCP_DBL>& fugId = fug[Id];

    vector<OCP_DBL>& ks = Ks[0];
    for (k = 0; k < 2; k++) {

        ks = Kw[k];
        iter = 0;
        flag = OCP_FALSE;
        Tsol = OCP_FALSE;
        while (OCP_TRUE) {
            Yt = 0;
            for (USI i = 0; i < NC; i++) {
                Y[i] = xj[i] * ks[i];
                Yt += Y[i];
            }
            Dscalar(NC, 1 / Yt, &Y[0]);

            if (iter > 0) {
                dY = 0;
                for (USI i = 0; i < NC; i++) {
                    dY += pow((Y[i] - di[i]), 2);
                }
                if (dY < dYtol) {
                    // converges
                    flag = OCP_TRUE;
                    break;
                }
            }
            eos->CalFug(P, T, &Y[0], &fugSta[0]);
            Se = 0;
            Sk = 0;
            for (USI i = 0; i < NC; i++) {
                di[i] = fugId[i] / (fugSta[i] * Yt);
                ks[i] *= di[i];
                Se += pow(di[i] - 1, 2);
                // Sk += pow(ks[i] - 1, 2);
                Sk += pow(log(ks[i]), 2);
            }

            iter++;
            if (Se < Stol) {
                flag = OCP_TRUE;
                break;
            }
            if (Sk < Ktol) {
                // Sk < Ktol -> trivial solution
                flag = OCP_TRUE;
                Tsol = OCP_TRUE;
                break;
            }
            if (iter > maxIt) {
                break;
            }

            // Record last Y with di
            di = Y;
        }

        if (!Tsol) {
            // flag = StableNR(Id);
        }

        flashCtrl.SSMsta.curIt += iter;

        if (flag && Yt > 1 - 0.1 && Sk > 1) {
            // close to phase boundary, or more than 1 phase, So don't skip at next step
            flagSkip = OCP_FALSE;
        }
        if (flag && Yt > 1 + eYt) {
            flashCtrl.SSMsta.conflag = OCP_TRUE;
            flashCtrl.SSMsta.res = sqrt(Se);
            return OCP_FALSE;
        }
    }
    flashCtrl.SSMsta.res = sqrt(Se);
    return OCP_TRUE;
}

OCP_BOOL OCPPhaseEquilibrium::StableSSM(const USI& Id)
{
    // if unsatble, return OCP_FALSE
    // if stable, return OCP_TRUE

    const vector<OCP_DBL>& xj = x[Id];
    eos->CalFugPhi(P, T, &x[Id][0], &fug[Id][0], &phi[Id][0]);
    const vector<OCP_DBL>& fugId = fug[Id];

    for (USI i = 0; i < NC; i++) {
        di[i] = phi[Id][i] * xj[i];
    }

    OCP_DBL  Stol = flashCtrl.SSMsta.tol2;
    USI      maxIt = flashCtrl.SSMsta.maxIt;
    OCP_DBL  eYt = flashCtrl.SSMsta.eYt;
    OCP_DBL  Se;
    OCP_BOOL flag;
    USI      iter;
    USI      k;

    for (k = 0; k < 4; k++) {

        Yt = 0;
        for (USI i = 0; i < NC; i++) {
            Y[i] = xj[i] / Kw[k][i];
            Yt += Y[i];
        }
        Dscalar(NC, 1 / Yt, &Y[0]);

        eos->CalFugPhi(P, T, &Y[0], &fugSta[0], &phiSta[0]);

        Se = 0;
        for (USI i = 0; i < NC; i++) {
            Se += pow(log(fugSta[i] / fugId[i] * Yt), 2);
        }

        flag = OCP_TRUE;
        iter = 0;

        while (Se > Stol) {

            Yt = 0;
            for (USI i = 0; i < NC; i++) {
                Y[i] = di[i] / phiSta[i];
                Yt += Y[i];
            }
            Dscalar(NC, 1 / Yt, &Y[0]);

            eos->CalFugPhi(P, T, &Y[0], &fugSta[0], &phiSta[0]);

            Se = 0;
            for (USI i = 0; i < NC; i++) {
                Se += pow(log(fugSta[i] / fugId[i] * Yt), 2);
            }

            iter++;
            if (iter > maxIt) {
                flag = OCP_FALSE;
                break;
            }
        }
        flashCtrl.SSMsta.curIt += iter;
        flag = StableNR(Id);
        // here, a relaxation is needed, on the one hand it can prevent the influence
        // of rounding error, on the other hand, if Yt is too close to 1, phase
        // splitting calculation may get into trouble and single phase is indentified
        // finally
        if (flag && Yt > 1 + eYt) {
            flashCtrl.SSMsta.conflag = OCP_TRUE;
            flashCtrl.SSMsta.res = sqrt(Se);
            return OCP_FALSE;
        }
    }
    flashCtrl.SSMsta.res = sqrt(Se);
    return OCP_TRUE;
}

OCP_BOOL OCPPhaseEquilibrium::StableNR(const USI& Id)
{

    for (USI i = 0; i < NC; i++) {
        resSTA[i] = log(fug[Id][i] / (fugSta[i] * Yt));
    }

    USI     maxIt = flashCtrl.NRsta.maxIt;
    OCP_DBL Stol = flashCtrl.NRsta.tol;
    OCP_DBL Se = Dnorm2(NC, &resSTA[0]);
    OCP_DBL alpha = 1;
    USI     iter = 0;

    while (Se > Stol) {

        eos->CalLnFugX(P, T, &Y[0], &lnfugX[0][0]);
        AssembleJmatSTA();
        // LUSolve(1, NC, &JmatSTA[0], &resSTA[0], &pivot[0]);
        SYSSolve(1, &uplo, NC, &JmatSTA[0], &resSTA[0], &pivot[0], &JmatWork[0],
            lenJmatWork);

        const OCP_DBL lYt = Yt;
        Yt = 0;
        for (USI i = 0; i < NC; i++) {
            Y[i] = Y[i] * lYt + alpha * resSTA[i];
            Yt += Y[i];
        }
        Dscalar(NC, 1 / Yt, &Y[0]);

        eos->CalFug(P, T, &Y[0], &fugSta[0]);
        for (USI i = 0; i < NC; i++) {
            resSTA[i] = log(fug[Id][i] / (fugSta[i] * Yt));
        }
        Se = Dnorm2(NC, &resSTA[0]);
        iter++;
        if (iter > maxIt) {
            flashCtrl.NRsta.curIt += iter;
            flashCtrl.NRsta.conflag = OCP_FALSE;
            flashCtrl.NRsta.res = Se;
            return OCP_FALSE;
        }
    }
    flashCtrl.NRsta.curIt += iter;
    flashCtrl.NRsta.conflag = OCP_TRUE;
    flashCtrl.NRsta.res = Se;
    return OCP_TRUE;
}

void OCPPhaseEquilibrium::AssembleJmatSTA()
{
    vector<OCP_DBL>& fugx = lnfugX[0];
    fill(JmatSTA.begin(), JmatSTA.end(), 0.0);
    OCP_DBL tmp;
    for (USI i = 0; i < NC; i++) {

        tmp = 0;
        for (USI k = 0; k < NC; k++) {
            tmp += Y[k] * fugx[i * NC + k];
        }

        for (USI j = 0; j < NC; j++) {
            // Symmetric
            // JmatSTA[i * NC + j] = (fugx[i * NC + j] - tmp + delta(i, j) / Y[i]) / Yt;
            JmatSTA[i * NC + j] = (fugx[i * NC + j] - tmp + 1) / Yt;
        }
    }
}


///////////////////////////////////////////////
// Phase-Splitting Analysis
///////////////////////////////////////////////


void OCPPhaseEquilibrium::PhaseSplit()
{
    flashCtrl.SSMsp.conflag = OCP_FALSE;
    flashCtrl.SSMsp.curIt = 0;
    flashCtrl.NRsp.conflag = OCP_FALSE;
    flashCtrl.NRsp.curIt = 0;
    flashCtrl.RR.curIt = 0;

    // cout << "begin" << endl;
    SplitSSM(OCP_FALSE);
    SplitNR();
    while (!flashCtrl.NRsp.conflag) {
        SplitSSM(OCP_TRUE);
        SplitNR();
        if (!CheckSplit()) break;
        if (flashCtrl.SSMsp.conflag) break;
    }
    CheckSplit();

    itersSSMSP += flashCtrl.SSMsp.curIt;
    itersNRSP += flashCtrl.NRsp.curIt;
    itersRR += flashCtrl.RR.curIt;
    countsSSMSP++;
    countsNRSP++;
}


void OCPPhaseEquilibrium::SplitSSM(const OCP_BOOL& flag)
{
    if (NP == 2) {
        SplitSSM2(flag);
    }
    else {
        SplitSSM3(flag);
    }
}

void OCPPhaseEquilibrium::SplitSSM2(const OCP_BOOL& flag)
{
    // NP = 2 in this case
    // Ks is very IMPORTANT!
    // flag = OCP_TRUE : Restart SSM
    // flag = OCP_FALSE : New SSM
    flashCtrl.SSMsp.conflag = OCP_TRUE;
    OCP_DBL Se = 1;
    OCP_DBL Stol = flashCtrl.SSMsp.tol2;
    USI     maxIt = flashCtrl.SSMsp.maxIt;

    if (!flag) {
        if (lNP == 2) {
            Ks[NP - 2] = lKs;
        }
        else {
            if (Yt < 1.1 || OCP_TRUE) {
                Ks[NP - 2] = Kw[0];
            }
            else {
                for (USI i = 0; i < NC; i++) {
                    // Ks[NP - 2][i] = phi[testPId][i] / phiSta[i];
                    Ks[NP - 2][i] = Y[i] / x[testPId][i];
                }
            }
        }
    }
    else {
        Stol *= 1E-1;
        maxIt *= 2;
    }

    if (Yt < 1.1) {
        Stol *= 1E-1;
        maxIt *= 2;
    }

    USI iter = 0;
    while (Se > Stol) {

        RachfordRice2();
        UpdateXRR();
        for (USI j = 0; j < NP; j++)
            eos->CalFugPhi(P, T, &x[j][0], &fug[j][0], &phi[j][0]);

        Se = 0;
        for (USI i = 0; i < NC; i++) {
            Se += pow(fug[1][i] / fug[0][i] - 1, 2);
            Ks[0][i] = phi[1][i] / phi[0][i];
        }

        iter++;
        if (iter > maxIt) {
            // OCP_WARNING("SSM not converged in Phase Spliting!");
            flashCtrl.SSMsp.conflag = OCP_FALSE;
            break;
        }
    }

    flashCtrl.SSMsp.res = sqrt(Se);
    flashCtrl.SSMsp.curIt += iter;
}

void OCPPhaseEquilibrium::SplitSSM3(const OCP_BOOL& flag) {}

void OCPPhaseEquilibrium::RachfordRice2() ///< Used when NP = 2
{
    const vector<OCP_DBL>& Ktmp = Ks[0];
    OCP_DBL                Kmin = Ktmp[0];
    OCP_DBL                Kmax = Ktmp[0];

    for (USI i = 1; i < NC; i++) {
        if (Ktmp[i] < Kmin) Kmin = Ktmp[i];
        if (Ktmp[i] > Kmax) Kmax = Ktmp[i];
    }

    const OCP_DBL numin = 1 / (1 - Kmax);
    const OCP_DBL numax = 1 / (1 - Kmin);

    nu[0] = 0.5 * (numin + numax);

    // Solve RR with NR
    OCP_DBL tmp, rj, J, dnuj, tmpnu;

    USI           iter = 0;
    const OCP_DBL tol = flashCtrl.RR.tol;
    const OCP_DBL maxIt = flashCtrl.RR.maxIt;
    while (OCP_TRUE) {

        rj = 0;
        J = 0;
        for (USI i = 0; i < NC; i++) {
            tmp = 1 + nu[0] * (Ktmp[i] - 1);
            rj += zi[i] * (Ktmp[i] - 1) / tmp;
            J -= zi[i] * (Ktmp[i] - 1) * (Ktmp[i] - 1) / (tmp * tmp);
        }

        if (fabs(rj) < tol || iter > maxIt) break;

        dnuj = -rj / J;
        tmpnu = nu[0] + dnuj;
        if (tmpnu < numax && tmpnu > numin) {
            nu[0] = tmpnu;
        }
        else {
            if (dnuj > 0) {
                nu[0] = (nu[0] + numax) / 2;
            }
            else {
                nu[0] = (nu[0] + numin) / 2;
            }
        }
        iter++;
    }

    flashCtrl.RR.curIt += iter;
    nu[1] = 1 - nu[0];

    // cout << scientific << setprecision(6) << nu[0] << "   " << nu[1] << endl;
}

void OCPPhaseEquilibrium::RachfordRice2P()
{
    // modified RachfordRice equations
    // less iterations but more divergence --- unstable!

    const vector<OCP_DBL>& Ktmp = Ks[0];
    OCP_DBL                Kmin = Ktmp[0];
    OCP_DBL                Kmax = Ktmp[0];

    for (USI i = 1; i < NC; i++) {
        if (Ktmp[i] < Kmin) Kmin = Ktmp[i];
        if (Ktmp[i] > Kmax) Kmax = Ktmp[i];
    }

    const OCP_DBL numin = 1 / (1 - Kmax);
    const OCP_DBL numax = 1 / (1 - Kmin);

    nu[0] = 0.5 * (numin + numax);

    // Solve RR with NR
    OCP_DBL tmp, rj, J, dnuj, tmpnu;
    OCP_DBL f, df;

    USI           iter = 0;
    const OCP_DBL tol = flashCtrl.RR.tol;
    const OCP_DBL maxIt = flashCtrl.RR.maxIt;
    while (OCP_TRUE) {

        rj = 0;
        J = 0;
        for (USI i = 0; i < NC; i++) {
            tmp = 1 + nu[0] * (Ktmp[i] - 1);
            rj += zi[i] * (Ktmp[i] - 1) / tmp;
            J -= zi[i] * (Ktmp[i] - 1) * (Ktmp[i] - 1) / (tmp * tmp);
        }
        f = (nu[0] - numin) * (numax - nu[0]);
        df = -2 * nu[0] + (numax + numin);
        J *= f;
        J += df * rj;
        rj *= f;

        if (fabs(rj) < tol || iter > maxIt) break;

        dnuj = -rj / J;
        tmpnu = nu[0] + dnuj;
        if (tmpnu < numax && tmpnu > numin) {
            nu[0] = tmpnu;
        }
        else {
            if (dnuj > 0) {
                nu[0] = (nu[0] + numax) / 2;
            }
            else {
                nu[0] = (nu[0] + numin) / 2;
            }
        }
        iter++;
    }

    flashCtrl.RR.curIt += iter;

    cout << iter << "      " << scientific << setprecision(3) << fabs(rj) << "      "
        << nu[0] << "      " << numin << "      " << numax << endl;

    nu[1] = 1 - nu[0];
}

void OCPPhaseEquilibrium::RachfordRice3() ///< Used when NP > 2
{
}

void OCPPhaseEquilibrium::UpdateXRR()
{
    OCP_DBL tmp = 0;
    for (USI i = 0; i < NC; i++) {
        tmp = 1;
        for (USI j = 0; j < NP - 1; j++) {
            tmp += nu[j] * (Ks[j][i] - 1);
        }
        x[NP - 1][i] = zi[i] / tmp;
        for (USI j = 0; j < NP - 1; j++) {
            x[j][i] = Ks[j][i] * x[NP - 1][i];
        }
    }
}

void OCPPhaseEquilibrium::SplitBFGS()
{
    // Use BFGS to calculate phase splitting
    // JmatSP will store the BFGS mat
    // resSP will store the resSP - lresSP if necessary
    // n will store the n - ln if necessary

    // get initial value, n and ln, resSP and lresSP, H0

    // begin BFGS
}

void OCPPhaseEquilibrium::SplitNR()
{
    flashCtrl.NRsp.conflag = OCP_FALSE;
    // for (USI j = 0; j < NP; j++) {
    //     nu[j] = fabs(nu[j]);
    // }

    USI len = NC * (NP - 1);
    x2n();
    CalResSP();
    OCP_DBL eNR0;
    OCP_DBL eNR = Dnorm2(len, &resSP[0]);
    OCP_DBL NRtol = flashCtrl.NRsp.tol;
    OCP_DBL alpha;

    OCP_DBL en;
    USI     iter = 0;
    eNR0 = eNR;
    while (eNR > NRtol) {

        // eNR0 = eNR;
        ln = n;
        for (USI j = 0; j < NP; j++)
            eos->CalLnFugN(P, T, &x[j][0], nu[j], &lnfugN[j][0]);

        AssembleJmatSP();

        // LUSolve(1, len, &JmatSP[0], &resSP[0], &pivot[0]);
        // PrintDX(NC, &resSP[0]);

        int info = SYSSolve(1, &uplo, len, &JmatSP[0], &resSP[0], &pivot[0],
            &JmatWork[0], lenJmatWork);
        if (info > 0 && false) {
            for (USI i = 0; i < len; i++) {
                for (USI j = 0; j < len; j++) {
                    cout << scientific << setprecision(9) << JmatSP[i * len + j]
                        << "   ";
                }
                cout << ";\n";
            }
        }
        // PrintDX(NC, &resSP[0]);

        alpha = CalStepNRsp();

        n[NP - 1] = zi;
        nu[NP - 1] = 1;
        for (USI j = 0; j < NP - 1; j++) {
            nu[j] = 0;
            for (USI i = 0; i < NC; i++) {
                n[j][i] += alpha * resSP[j * NC + i];
                n[NP - 1][i] -= n[j][i];
                nu[j] += n[j][i];
            }
            nu[NP - 1] -= nu[j];

            for (USI i = 0; i < NC; i++) {
                // x[j][i] = n[j][i] / nu[j];
                x[j][i] = fabs(n[j][i] / nu[j]);
            }
        }
        for (USI i = 0; i < NC; i++) {
            x[NP - 1][i] = fabs(n[NP - 1][i] / nu[NP - 1]);
        }

        for (USI j = 0; j < NP; j++)
            eos->CalFug(P, T, &x[j][0], &fug[j][0]);

        CalResSP();
        eNR = Dnorm2(len, &resSP[0]);
        iter++;
        if (eNR > eNR0 || iter > flashCtrl.NRsp.maxIt) {
            break;
        }

        // Maybe it should be execute before "eNR > eNR0 || iter > EoSctrl.NRsp.maxIt"
        en = 0;
        for (USI j = 0; j < NP; j++) {
            for (USI i = 0; i < NC; i++) {
                en += (ln[j][i] - n[j][i]) * (ln[j][i] - n[j][i]);
            }
        }
        if (en < 1E-16 * (NP * NC) * (NP * NC)) {
            flashCtrl.NRsp.conflag = OCP_TRUE;
            break;
        }
    }
    flashCtrl.NRsp.res = eNR;
    if (eNR < NRtol) flashCtrl.NRsp.conflag = OCP_TRUE;
    flashCtrl.NRsp.curIt += iter;

    // cout << iter << "   " << scientific << setprecision(3) << eNR << endl;
}

void OCPPhaseEquilibrium::CalResSP()
{
    // So it equals -res
    for (USI j = 0; j < NP - 1; j++) {
        for (USI i = 0; i < NC; i++) {
            resSP[j * NC + i] = log(fug[NP - 1][i] / fug[j][i]);
        }
    }
}



void OCPPhaseEquilibrium::AssembleJmatSP()
{
    // Dim: (NP-1)*NC
    // Attention that fugNj is sysmetric
    fill(JmatSP.begin(), JmatSP.end(), 0);

    OCP_DBL* Jtmp = &JmatSP[0];
    const OCP_DBL* fugNp;
    const OCP_DBL* fugNj;

    for (USI j = 0; j < NP - 1; j++) {
        fugNp = &lnfugN[NP - 1][0];
        fugNj = &lnfugN[j][0];

        for (USI i = 0; i < NC; i++) {
            // ith components
            for (USI k = 0; k < NC; k++) {
                // kth fugacity
                Jtmp[k] = fugNj[k] + fugNp[k];
            }
            Jtmp += NC * (NP - 1);
            fugNp += NC;
            fugNj += NC;
        }
        Jtmp += NC;
    }
}

OCP_DBL OCPPhaseEquilibrium::CalStepNRsp()
{
    OCP_DBL alpha = 1;
    OCP_DBL tmp;

    for (USI j = 0; j < NP - 1; j++) {

        const OCP_DBL* nj = &n[j][0];
        const OCP_DBL* r = &resSP[j * NC];

        for (USI i = 0; i < NC; i++) {
            tmp = nj[i] + alpha * r[i];
            if (tmp < 0) {
                alpha *= 0.9 * fabs(nj[i] / r[i]);
            }
        }
    }
    return alpha;
}


OCP_BOOL OCPPhaseEquilibrium::CheckSplit()
{
    if (NP == 2) {

        OCP_DBL eX = 0;
        OCP_DBL nuMax = max(nu[0], nu[1]);

        for (USI i = 0; i < NC; i++) {
            eX += (x[0][i] - x[1][i]) * (x[0][i] - x[1][i]);
        }

        if (OCP_TRUE) {
            // Calculate Gibbs Energy

            eos->CalFug(P, T, &zi[0], &fugSta[0]);
            GibbsEnergyB = 0;
            GibbsEnergyE = 0;
            for (USI i = 0; i < NC; i++) {
                GibbsEnergyB += zi[i] * log(fugSta[i]);
                GibbsEnergyE += (n[0][i] * log(fug[0][i]) + n[1][i] * log(fug[1][i]));
            }

            cout << scientific << setprecision(6);
            // cout << GibbsEnergyB << "   " << GibbsEnergyE << endl;
            if (GibbsEnergyE > GibbsEnergyB) {
                cout << ftype << "   ";
                cout << GibbsEnergyB << "   ";
                cout << GibbsEnergyE << "   ";
                cout << nuMax << "   ";
                cout << eX << "   ";
                cout << flashCtrl.NRsp.conflag << "   ";
                cout << bulkId << endl;
            }
        }

        if (nuMax < 1 && flashCtrl.NRsp.conflag && isfinite(eX)) {
            // accept this result
        }
        else {
            if (!isfinite(eX) || 1 - nuMax < 1E-3) {
                NP = 1;
                x[0] = zi;
                nu[0] = 1;

                flashCtrl.SSMsta.conflag = OCP_FALSE;
                flashCtrl.NRsta.conflag = OCP_FALSE;
                return OCP_FALSE;
            }
        }
    }
    return OCP_TRUE;
}


///////////////////////////////////////////////
// Statistcs
///////////////////////////////////////////////

void OCPPhaseEquilibrium::OutMixtureIters() const
{
    cout << "SSMSTA:     " << setw(12) << itersSSMSTA << setw(15)
        << itersSSMSTA * 1.0 / countsSSMSTA << endl;
    cout << "NRSTA:      " << setw(12) << itersNRSTA << setw(15)
        << itersNRSTA * 1.0 / countsNRSTA << endl;
    cout << "SSMSP:      " << setw(12) << itersSSMSP << setw(15)
        << itersSSMSP * 1.0 / countsSSMSP << endl;
    cout << "NRSP:       " << setw(12) << itersNRSP << setw(15)
        << itersNRSP * 1.0 / countsNRSP << endl;
    cout << "NRRR:       " << setw(12) << itersRR << setw(15)
        << itersRR * 1.0 / itersSSMSP << endl;
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/28/2023      Create file                          */
/*----------------------------------------------------------------------------*/