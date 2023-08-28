/*! \file    WellPeaceman.cpp
 *  \brief   Peaceman PeacemanWell class definition
 *  \author  Shizhe Li
 *  \date    Aug/17/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "WellPeaceman.hpp"

 /// It calculates pressure difference between perforations iteratively.
 /// This function can be used in both black oil model and compositional model.
 /// stability of this method shoule be tested.
void PeacemanWell::CaldG(const Bulk& bk)
{
    OCP_FUNCNAME;

    if (opt.type == WellType::injector)
        CalInjdG(bk);
    else
        CalProddG01(bk);
}

void PeacemanWell::CalInjdG(const Bulk& bk)
{
    OCP_FUNCNAME;

    const OCP_DBL   maxlen = 10;
    USI             seg_num = 0;
    OCP_DBL         seg_len = 0;
    vector<OCP_DBL> dGperf(numPerf, 0);

    if (depth <= perf.front().depth) {
        // PeacemanWell is higher
        for (OCP_INT p = numPerf - 1; p >= 0; p--) {
            if (p == 0) {
                seg_num = ceil(fabs((perf[0].depth - depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[0].depth - depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p].depth - perf[p - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p].depth - perf[p - 1].depth) / seg_num;
            }
            OCP_USI n = perf[p].location;
            perf[p].P = bhp + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;
            auto    PVT = bk.PVTm.GetPVT(n);

            for (USI i = 0; i < seg_num; i++) {
                Ptmp -= PVT->RhoPhase(Ptmp, 0, opt.injTemp, opt.injZi, opt.injPhase) * GRAVITY_FACTOR * seg_len;
            }
            dGperf[p] = Pperf - Ptmp;
        }
        dG[0] = dGperf[0];
        for (USI p = 1; p < numPerf; p++) {
            dG[p] = dG[p - 1] + dGperf[p];
        }
    }
    else if (depth >= perf[numPerf - 1].depth) {
        // PeacemanWell is lower
        for (USI p = 0; p < numPerf; p++) {
            if (p == numPerf - 1) {
                seg_num = ceil(fabs((depth - perf[numPerf - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (depth - perf[numPerf - 1].depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p + 1].depth - perf[p].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p + 1].depth - perf[p].depth) / seg_num;
            }
            OCP_USI n = perf[p].location;
            perf[p].P = bhp + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;

            auto PVT = bk.PVTm.GetPVT(n);
            for (USI i = 0; i < seg_num; i++) {
                Ptmp += PVT->RhoPhase(Ptmp, 0, opt.injTemp, opt.injZi, opt.injPhase) * GRAVITY_FACTOR * seg_len;
            }
            dGperf[p] = Ptmp - Pperf;
        }
        dG[numPerf - 1] = dGperf[numPerf - 1];
        for (OCP_INT p = numPerf - 2; p >= 0; p--) {
            dG[p] = dG[p + 1] + dGperf[p];
        }
    }
}

// Use transj
void PeacemanWell::CalProddG01(const Bulk& bk)
{
    OCP_FUNCNAME;

    const BulkVarSet& bvs = bk.vs;

    const OCP_DBL   maxlen = 10;
    USI             seg_num = 0;
    OCP_DBL         seg_len = 0;
    vector<OCP_DBL> dGperf(numPerf, 0);
    vector<OCP_DBL> tmpNi(nc, 0);
    OCP_DBL         rhotmp, qtacc, rhoacc;

    if (depth <= perf.front().depth) {
        // PeacemanWell is higher
        for (OCP_INT p = numPerf - 1; p >= 0; p--) {
            if (p == 0) {
                seg_num = ceil(fabs((perf[0].depth - depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[0].depth - depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p].depth - perf[p - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p].depth - perf[p - 1].depth) / seg_num;
            }
            OCP_USI n = perf[p].location;
            perf[p].P = bhp + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;

            // fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (USI j = 0; j < np; j++) {
                const OCP_USI n_np_j = n * np + j;
                if (!bvs.phaseExist[n_np_j]) continue;
                for (USI k = 0; k < nc; k++) {
                    tmpNi[k] += (bvs.P[n] - perf[p].P) * perf[p].transj[j] *
                        bvs.xi[n_np_j] * bvs.xij[n_np_j * nc + k];
                }
            }
            OCP_DBL tmpSum = Dnorm1(nc, &tmpNi[0]);
            if (tmpSum < TINY) {
                for (USI i = 0; i < nc; i++) {
                    tmpNi[i] = bvs.Ni[n * nc + i];
                }
            }

            auto PVT = bk.PVTm.GetPVT(n);
            for (USI i = 0; i < seg_num; i++) {
                qtacc = rhoacc = 0;
                PVT->Flash(Ptmp, bvs.T[n], tmpNi.data());
                for (USI j = 0; j < np; j++) {
                    if (PVT->GetPhaseExist(j)) {
                        rhotmp = PVT->GetRho(j);
                        qtacc += PVT->GetVj(j);
                        rhoacc += PVT->GetVj(j) * rhotmp;
                    }
                }
                Ptmp -= rhoacc / qtacc * seg_len * GRAVITY_FACTOR;
            }
            dGperf[p] = Pperf - Ptmp;
        }
        dG[0] = dGperf[0];
        for (USI p = 1; p < numPerf; p++) {
            dG[p] = dG[p - 1] + dGperf[p];
        }
    }
    else if (depth >= perf[numPerf - 1].depth) {
        // PeacemanWell is lower
        for (USI p = 0; p < numPerf; p++) {
            if (p == numPerf - 1) {
                seg_num = ceil(fabs((depth - perf[numPerf - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (depth - perf[numPerf - 1].depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p + 1].depth - perf[p].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p + 1].depth - perf[p].depth) / seg_num;
            }
            OCP_USI n = perf[p].location;
            perf[p].P = bhp + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;

            // fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (USI j = 0; j < np; j++) {
                const OCP_USI n_np_j = n * np + j;
                if (!bvs.phaseExist[n_np_j]) continue;
                for (USI k = 0; k < nc; k++) {
                    tmpNi[k] += (bvs.P[n] - perf[p].P) * perf[p].transj[j] *
                        bvs.xi[n_np_j] * bvs.xij[n_np_j * nc + k];
                }
            }
            OCP_DBL tmpSum = Dnorm1(nc, &tmpNi[0]);
            if (tmpSum < TINY) {
                for (USI i = 0; i < nc; i++) {
                    tmpNi[i] = bvs.Ni[n * nc + i];
                }
            }

            auto PVT = bk.PVTm.GetPVT(n);
            for (USI i = 0; i < seg_num; i++) {
                qtacc = rhoacc = 0;
                PVT->Flash(Ptmp, bvs.T[n], tmpNi.data());
                for (USI j = 0; j < np; j++) {
                    if (PVT->GetPhaseExist(j)) {
                        rhotmp = PVT->GetRho(j);
                        qtacc += PVT->GetVj(j);
                        rhoacc += PVT->GetVj(j) * rhotmp;
                    }
                }
                Ptmp += rhoacc / qtacc * seg_len * GRAVITY_FACTOR;
            }
            dGperf[p] = Ptmp - Pperf;
        }
        dG[numPerf - 1] = dGperf[numPerf - 1];
        for (OCP_INT p = numPerf - 2; p >= 0; p--) {
            dG[p] = dG[p + 1] + dGperf[p];
        }
    }
}

// Use bulk
void PeacemanWell::CalProddG02(const Bulk& bk)
{
    OCP_FUNCNAME;

    const BulkVarSet& bvs = bk.vs;

    const OCP_DBL   maxlen = 10;
    USI             seg_num = 0;
    OCP_DBL         seg_len = 0;
    vector<OCP_DBL> dGperf(numPerf, 0);
    vector<OCP_DBL> tmpNi(nc, 0);
    OCP_DBL         rhotmp, qtacc, rhoacc;

    if (depth <= perf.front().depth) {
        // PeacemanWell is higher
        for (OCP_INT p = numPerf - 1; p >= 0; p--) {
            if (p == 0) {
                seg_num = ceil(fabs((perf[0].depth - depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[0].depth - depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p].depth - perf[p - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p].depth - perf[p - 1].depth) / seg_num;
            }
            OCP_USI n = perf[p].location;
            perf[p].P = bhp + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;

            fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (USI j = 0; j < np; j++) {
                const OCP_USI n_np_j = n * np + j;
                if (!bvs.phaseExist[n_np_j]) continue;
                for (USI k = 0; k < nc; k++) {
                    tmpNi[k] += (perf[p].transj[j] > 0) * bvs.xi[n_np_j] *
                        bvs.xij[n_np_j * nc + k];
                }
            }
            OCP_DBL tmpSum = Dnorm1(nc, &tmpNi[0]);
            if (tmpSum < TINY) {
                for (USI i = 0; i < nc; i++) {
                    tmpNi[i] = bvs.Ni[n * nc + i];
                }
            }

            auto PVT = bk.PVTm.GetPVT(n);
            for (USI i = 0; i < seg_num; i++) {
                qtacc = rhoacc = 0;
                PVT->Flash(Ptmp, bvs.T[n], tmpNi.data());
                for (USI j = 0; j < np; j++) {
                    if (PVT->GetPhaseExist(j)) {
                        rhotmp = PVT->GetRho(j);
                        qtacc += PVT->GetVj(j);
                        rhoacc += PVT->GetVj(j) * rhotmp;
                    }
                }
                Ptmp -= rhoacc / qtacc * seg_len * GRAVITY_FACTOR;
            }
            dGperf[p] = Pperf - Ptmp;
        }
        dG[0] = dGperf[0];
        for (USI p = 1; p < numPerf; p++) {
            dG[p] = dG[p - 1] + dGperf[p];
        }
    }
    else if (depth >= perf[numPerf - 1].depth) {
        // PeacemanWell is lower
        for (USI p = 0; p < numPerf; p++) {
            if (p == numPerf - 1) {
                seg_num = ceil(fabs((depth - perf[numPerf - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (depth - perf[numPerf - 1].depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p + 1].depth - perf[p].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p + 1].depth - perf[p].depth) / seg_num;
            }
            OCP_USI n = perf[p].location;
            perf[p].P = bhp + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;

            fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (USI j = 0; j < np; j++) {
                const OCP_USI n_np_j = n * np + j;
                if (!bvs.phaseExist[n_np_j]) continue;
                for (USI k = 0; k < nc; k++) {
                    tmpNi[k] += (perf[p].transj[j] > 0) * bvs.xi[n_np_j] *
                        bvs.xij[n_np_j * nc + k];
                }
            }
            OCP_DBL tmpSum = Dnorm1(nc, &tmpNi[0]);
            if (tmpSum < TINY) {
                for (USI i = 0; i < nc; i++) {
                    tmpNi[i] = bvs.Ni[n * nc + i];
                }
            }

            auto PVT = bk.PVTm.GetPVT(n);
            for (USI i = 0; i < seg_num; i++) {
                qtacc = rhoacc = 0;
                PVT->Flash(Ptmp, bvs.T[n], tmpNi.data());
                for (USI j = 0; j < np; j++) {
                    if (PVT->GetPhaseExist(j)) {
                        rhotmp = PVT->GetRho(j);
                        qtacc += PVT->GetVj(j);
                        rhoacc += PVT->GetVj(j) * rhotmp;
                    }
                }
                Ptmp += rhoacc / qtacc * seg_len * GRAVITY_FACTOR;
            }
            dGperf[p] = Ptmp - Pperf;
        }
        dG[numPerf - 1] = dGperf[numPerf - 1];
        for (OCP_INT p = numPerf - 2; p >= 0; p--) {
            dG[p] = dG[p + 1] + dGperf[p];
        }
    }
}

// Use qi_lbmol
void PeacemanWell::CalProddG(const Bulk& bk)
{
    OCP_FUNCNAME;

    const BulkVarSet& bvs = bk.vs;

    const OCP_DBL   maxlen = 5;
    USI             seg_num = 0;
    OCP_DBL         seg_len = 0;
    vector<OCP_DBL> tmpNi(nc, 0);
    vector<OCP_DBL> dGperf(numPerf, 0);
    OCP_DBL         qtacc = 0;
    OCP_DBL         rhoacc = 0;
    OCP_DBL         rhotmp = 0;

    if (depth <= perf.front().depth) {
        // PeacemanWell is higher

        // check qi_lbmol   ----   test
        if (perf[numPerf - 1].state == WellState::close) {
            for (OCP_INT p = numPerf - 2; p >= 0; p--) {
                if (perf[p].state == WellState::open) {
                    for (USI i = 0; i < nc; i++) {
                        perf[numPerf - 1].qi_lbmol[i] = perf[p].qi_lbmol[i];
                    }
                    break;
                }
            }
        }

        for (OCP_INT p = numPerf - 1; p >= 0; p--) {

            if (p == 0) {
                seg_num = ceil(fabs((perf[0].depth - depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[0].depth - depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p].depth - perf[p - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p].depth - perf[p - 1].depth) / seg_num;
            }

            OCP_USI n = perf[p].location;
            perf[p].P = bhp + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;

            // tmpNi.assign(nc, 0);
            // for (OCP_INT p1 = numPerf - 1; p1 >= p; p1--) {
            //    for (USI i = 0; i < nc; i++) {
            //        tmpNi[i] += perf[p1].qi_lbmol[i];
            //    }
            //}

            for (USI i = 0; i < nc; i++) {
                tmpNi[i] += perf[p].qi_lbmol[i];
            }

            // check tmpNi
            for (auto& v : tmpNi) {
                v = fabs(v);
            }

            auto PVT = bk.PVTm.GetPVT(n);
            for (USI k = 0; k < seg_num; k++) {
                PVT->Flash(Ptmp, bvs.T[n], tmpNi.data());
                for (USI j = 0; j < np; j++) {
                    if (PVT->GetPhaseExist(j)) {
                        rhotmp = PVT->GetRho(j);
                        qtacc += PVT->GetVj(j) / seg_num;
                        rhoacc += PVT->GetVj(j) * rhotmp / seg_num;
#ifdef DEBUG
                        if (rhotmp <= 0 || !isfinite(rhotmp)) {
                            OCP_ABORT("Wrong rho " + to_string(rhotmp));
                        }
#endif // DEBUG
                    }
                }
                Ptmp -= rhoacc / qtacc * seg_len * GRAVITY_FACTOR;
            }
            dGperf[p] = Pperf - Ptmp;
        }
        dG[0] = dGperf[0];
        for (USI p = 1; p < numPerf; p++) {
            dG[p] = dG[p - 1] + dGperf[p];
        }
    }
    else if (depth >= perf.back().depth) {
        // PeacemanWell is lower

        // check qi_lbmol   ----   test
        if (perf[0].state == WellState::close) {
            for (USI p = 1; p <= numPerf; p++) {
                if (perf[p].state == WellState::open) {
                    for (USI i = 0; i < nc; i++) {
                        perf[numPerf - 1].qi_lbmol[i] = perf[p].qi_lbmol[i];
                    }
                    break;
                }
            }
        }

        for (USI p = 0; p < numPerf; p++) {
            if (p == numPerf - 1) {
                seg_num = ceil(fabs((depth - perf[numPerf - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (depth - perf[numPerf - 1].depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p + 1].depth - perf[p].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p + 1].depth - perf[p].depth) / seg_num;
            }

            OCP_USI n = perf[p].location;
            perf[p].P = bhp + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;


            fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (OCP_INT p1 = numPerf - 1; p1 - p >= 0; p1--) {
                for (USI i = 0; i < nc; i++) {
                    tmpNi[i] += perf[p1].qi_lbmol[i];
                }
            }

            // check tmpNi
            for (auto& v : tmpNi) {
                v = fabs(v);
            }

            auto PVT = bk.PVTm.GetPVT(n);
            for (USI k = 0; k < seg_num; k++) {
                PVT->Flash(Ptmp, bvs.T[n], tmpNi.data());
                for (USI j = 0; j < np; j++) {
                    if (PVT->GetPhaseExist(j)) {
                        rhotmp = PVT->GetRho(j);
                        qtacc += PVT->GetVj(j) / seg_num;
                        rhoacc += PVT->GetVj(j) * rhotmp / seg_num;
                    }
                }
                Ptmp += rhoacc / qtacc * seg_len * GRAVITY_FACTOR;
            }
            dGperf[p] = Ptmp - Pperf;
        }
        dG[numPerf - 1] = dGperf[numPerf - 1];
        for (OCP_INT p = numPerf - 2; p >= 0; p--) {
            dG[p] = dG[p + 1] + dGperf[p];
        }

    }
    else {
        OCP_ABORT("Wrong well position!");
    }
}


void PeacemanWellIsoT::AssembleMatFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{
    if (opt.state == WellState::open) {
        if (opt.type == WellType::injector)        AssembleMatInjFIM(ls, bk, dt);
        else if (opt.type == WellType::productor)  AssembleMatProdFIM(ls, bk, dt);
        else                                       OCP_ABORT("WRONG Well Type!");
    }
}


void PeacemanWellIsoT::CalResFIM(OCP_USI& wId, OCPRes& res, const Bulk& bk, const OCP_DBL& dt) const
{
    if (opt.state == WellState::open) {

        const USI len = nc + 1;
        // Well to Bulk
        for (USI p = 0; p < numPerf; p++) {
            const OCP_USI k = perf[p].location;
            for (USI i = 0; i < nc; i++) {
                res.resAbs[k * len + 1 + i] += perf[p].qi_lbmol[i] * dt;
            }
        }
        // Well Self
        switch (opt.mode)
        {
        case WellOptMode::bhp:
            res.resAbs[wId] = bhp - opt.tarBHP;
            break;
        case WellOptMode::irate:
        case WellOptMode::orate:
        case WellOptMode::grate:
        case WellOptMode::wrate:
        case WellOptMode::lrate:
            CalFactor(bk);
            if (opt.type == WellType::injector)  res.resAbs[wId] = opt.tarRate;
            else                                 res.resAbs[wId] = -opt.tarRate;
            for (USI i = 0; i < nc; i++) {
                res.resAbs[wId] += qi_lbmol[i] * factor[i];
            }
            res.maxWellRelRes_mol =
                max(res.maxWellRelRes_mol,
                    fabs(res.resAbs[wId] / opt.tarRate));
            break;
        default:
            OCP_ABORT("Wrong well opt mode!");
            break;
        }
        wId += len;
    }
}


void PeacemanWellIsoT::AssembleMatInjFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{
    const BulkVarSet& bvs = bk.vs;

    const USI ncol   = nc + 1;
    const USI ncol2  = np * nc + np;
    const USI bsize  = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
    OCP_USI   n_np_j;

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);

    OCP_DBL mu, muP, dP, transIJ;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, bmat);

    CalFactor(bk);

    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        dP = bvs.P[n] - bhp - dG[p];

        for (USI j = 0; j < np; j++) {
            n_np_j = n * np + j;
            if (!bvs.phaseExist[n_np_j]) continue;

            mu = bvs.mu[n_np_j];
            muP = bvs.muP[n_np_j];

            for (USI i = 0; i < nc; i++) {
                // dQ / dP
                transIJ = perf[p].transj[j] * perf[p].xi * opt.injZi[i];
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu);
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                for (USI k = 0; k < np; k++) {
                    dQdXsB[(i + 1) * ncol2 + k] +=
                        CONV1 * perf[p].WI * perf[p].multiplier * perf[p].xi *
                        opt.injZi[i] * bvs.dKrdS[n_np_j * np + k] * dP / mu;
                }
                // dQ / dxij
                for (USI k = 0; k < nc; k++) {
                    dQdXsB[(i + 1) * ncol2 + np + j * nc + k] +=
                        -transIJ * dP / mu * bvs.mux[n_np_j * nc + k];
                }
            }
        }

        // Bulk to Well
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bvs.dSec_dPri[n * bsize2], 1,
            bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Bulk - Bulk -- add
        ls.AddDiag(n, bmat);

        // Bulk - Well -- insert
        bmat = dQdXpW;
        Dscalar(bsize, dt, bmat.data());
        ls.NewOffDiag(n, wId, bmat);

        // Well
        switch (opt.mode) {
        case WellOptMode::irate:
        case WellOptMode::orate:
        case WellOptMode::grate:
        case WellOptMode::wrate:
        case WellOptMode::lrate:
            // Well - Well -- add
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < nc; i++) {
                bmat[0] += dQdXpW[(i + 1) * ncol] * factor[i];
                bmat[(i + 1) * ncol + i + 1] = 1;
            }
            ls.AddDiag(wId, bmat);

            // Well - Bulk -- insert
            bmat = dQdXpB;
            DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bvs.dSec_dPri[n * bsize2],
                1, bmat.data());
            fill(bmat2.begin(), bmat2.end(), 0.0);
            for (USI i = 0; i < nc; i++) {
                Daxpy(ncol, factor[i], bmat.data() + (i + 1) * ncol, bmat2.data());
            }
            ls.NewOffDiag(wId, n, bmat2);
            break;

        case WellOptMode::bhp:
            // Well - Well -- add
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < ncol; i++) {
                bmat[i * ncol + i] = 1;
            }
            ls.AddDiag(wId, bmat);

            // Well - Bulk -- insert
            fill(bmat.begin(), bmat.end(), 0.0);
            ls.NewOffDiag(wId, n, bmat);
            break;

        default:
            OCP_ABORT("Wrong Well Opt mode!");
            break;
        }
    }
}


void PeacemanWellIsoT::AssembleMatProdFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{
    const BulkVarSet& bvs = bk.vs;

    const USI ncol = nc + 1;
    const USI ncol2 = np * nc + np;
    const USI bsize = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
    OCP_USI   n_np_j;

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);

    OCP_DBL xij, xi, mu, muP, xiP, dP, transIJ, tmp;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, bmat);

    // Set Commponent Weight
    CalFactor(bk);

    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        for (USI j = 0; j < np; j++) {
            n_np_j = n * np + j;
            if (!bvs.phaseExist[n_np_j]) continue;

            dP = bvs.Pj[n_np_j] - bhp - dG[p];
            xi = bvs.xi[n_np_j];
            mu = bvs.mu[n_np_j];
            muP = bvs.muP[n_np_j];
            xiP = bvs.xiP[n_np_j];

            for (USI i = 0; i < nc; i++) {
                xij = bvs.xij[n_np_j * nc + i];
                // dQ / dP
                transIJ = perf[p].transj[j] * xi * xij;
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu) +
                    dP * perf[p].transj[j] * xij * xiP;
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                for (USI k = 0; k < np; k++) {
                    tmp = CONV1 * perf[p].WI * perf[p].multiplier * dP / mu * xi *
                        xij * bvs.dKrdS[n_np_j * np + k];
                    // capillary pressure
                    tmp += transIJ * bvs.dPcdS[n_np_j * np + k];
                    dQdXsB[(i + 1) * ncol2 + k] += tmp;
                }
                // dQ / dCij
                for (USI k = 0; k < nc; k++) {
                    tmp = dP * perf[p].transj[j] * xij *
                        (bvs.xix[n_np_j * nc + k] - xi / mu * bvs.mux[n_np_j * nc + k]);
                    dQdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                }
                dQdXsB[(i + 1) * ncol2 + np + j * nc + i] +=
                    perf[p].transj[j] * xi * dP;
            }
        }

        // Bulk - Bulk -- add
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bvs.dSec_dPri[n * bsize2], 1,
            bmat.data());

        Dscalar(bsize, dt, bmat.data());
        ls.AddDiag(n, bmat);

        // Bulk - Well -- insert
        bmat = dQdXpW;
        Dscalar(bsize, dt, bmat.data());
        ls.NewOffDiag(n, wId, bmat);

        // Well
        switch (opt.mode) {
        case WellOptMode::irate:
        case WellOptMode::orate:
        case WellOptMode::grate:
        case WellOptMode::wrate:
        case WellOptMode::lrate:
            // Well - Well -- add
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < nc; i++) {
                bmat[0] += dQdXpW[(i + 1) * ncol] * factor[i];
                bmat[(i + 1) * ncol + i + 1] = 1;
            }
            ls.AddDiag(wId, bmat);

            // Well - Bulk -- insert
            bmat = dQdXpB;
            DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bvs.dSec_dPri[n * bsize2],
                1, bmat.data());
            fill(bmat2.begin(), bmat2.end(), 0.0);
            for (USI i = 0; i < nc; i++) {
                Daxpy(ncol, factor[i], bmat.data() + (i + 1) * ncol,
                    bmat2.data());
            }
            ls.NewOffDiag(wId, n, bmat2);
            break;

        case WellOptMode::bhp:
            // Well - Well -- add
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < ncol; i++) {
                bmat[i * ncol + i] = 1;
            }
            ls.AddDiag(wId, bmat);

            // Well - Bulk -- insert
            fill(bmat.begin(), bmat.end(), 0.0);
            ls.NewOffDiag(wId, n, bmat);
            break;

        default:
            OCP_ABORT("Wrong Well Opt mode!");
            break;
        }
    }
}


void PeacemanWellIsoT::AssembleMatIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{
    if (opt.state == WellState::open) {
        if (opt.type == WellType::injector)        AssembleMatInjIMPEC(ls, bk, dt);
        else if (opt.type == WellType::productor)  AssembleMatProdIMPEC(ls, bk, dt);
        else                                       OCP_ABORT("WRONG Well Type!");
    }
}


void PeacemanWellIsoT::AssembleMatInjIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{
    const BulkVarSet& bvs = bk.vs;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, 0.0);
    CalFactor(bk);
    OCP_DBL Vfi_zi, valb, valw, bb, bw;

    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;

        Vfi_zi = 0;
        for (USI i = 0; i < nc; i++) {
            Vfi_zi += bvs.vfi[n * nc + i] * opt.injZi[i];
        }

        valw = dt * perf[p].xi * perf[p].transINJ;
        bw = valw * dG[p];
        valb = valw * Vfi_zi;
        bb = valb * dG[p];

        // Bulk to Well
        ls.AddDiag(n, valb);
        ls.NewOffDiag(n, wId, -valb);
        ls.AddRhs(n, bb);

        // Well to Bulk
        switch (opt.mode) {
        case WellOptMode::irate:
        case WellOptMode::orate:
        case WellOptMode::grate:
        case WellOptMode::wrate:
        case WellOptMode::lrate:
            ls.AddDiag(wId, valw * factor[0]);
            ls.NewOffDiag(wId, n, -valw * factor[0]);
            ls.AddRhs(wId, -bw * factor[0]);
            break;
        case WellOptMode::bhp:
            ls.NewOffDiag(wId, n, 0);
            break;
        default:
            OCP_ABORT("Wrong well option mode!");
        }
    }

    // Well Self
    switch (opt.mode) {
    case WellOptMode::irate:
    case WellOptMode::orate:
    case WellOptMode::grate:
    case WellOptMode::wrate:
    case WellOptMode::lrate:
        ls.AddRhs(wId, dt * opt.tarRate);
        break;
    case WellOptMode::bhp:
        ls.AddDiag(wId, dt);
        ls.AddRhs(wId, dt * opt.tarBHP);
        ls.AssignGuess(wId, opt.tarBHP);
        break;
    default:
        OCP_ABORT("Wrong well option mode!");
    }
}


void PeacemanWellIsoT::AssembleMatProdIMPEC(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{
    const BulkVarSet& bvs = bk.vs;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, 0.0);

    // Set Prod Weight
    CalFactor(bk);

    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;

        OCP_DBL valb = 0;
        OCP_DBL bb = 0;
        OCP_DBL valw = 0;
        OCP_DBL bw = 0;

        for (USI j = 0; j < np; j++) {
            const OCP_USI n_np_j = n * np + j;
            if (!bvs.phaseExist[n_np_j]) continue;

            OCP_DBL tempb = 0;
            OCP_DBL tempw = 0;

            for (USI i = 0; i < nc; i++) {
                tempb += bvs.vfi[n * nc + i] * bvs.xij[n_np_j * nc + i];
                tempw += factor[i] * bvs.xij[n_np_j * nc + i];
            }
            OCP_DBL trans = dt * perf[p].transj[j] * bvs.xi[n_np_j];
            valb += tempb * trans;
            valw += tempw * trans;

            OCP_DBL dP = dG[p] - bvs.Pc[n_np_j];
            bb += tempb * trans * dP;
            bw += tempw * trans * dP;
        }

        // Bulk to Well
        ls.AddDiag(n, valb);
        ls.NewOffDiag(n, wId, -valb);
        ls.AddRhs(n, bb);

        // Well to Bulk
        switch (opt.mode) {
        case WellOptMode::irate:
        case WellOptMode::orate:
        case WellOptMode::grate:
        case WellOptMode::wrate:
        case WellOptMode::lrate:
            ls.AddDiag(wId, -valw);
            ls.NewOffDiag(wId, n, valw);
            ls.AddRhs(wId, bw);
            break;
        case WellOptMode::bhp:
            ls.NewOffDiag(wId, n, 0.0);
            break;
        default:
            OCP_ABORT("Wrong well option mode!");
        }
    }

    // Well Self
    switch (opt.mode) {
    case WellOptMode::irate:
    case WellOptMode::orate:
    case WellOptMode::grate:
    case WellOptMode::wrate:
    case WellOptMode::lrate:
        ls.AddRhs(wId, dt * opt.tarRate);
        break;
    case WellOptMode::bhp:
        ls.AddDiag(wId, dt);
        ls.AddRhs(wId, dt * opt.tarBHP);
        ls.AssignGuess(wId, opt.tarBHP);
        break;
    default:
        OCP_ABORT("Wrong well option mode!");
    }
}


void PeacemanWellT::CalResFIM(OCP_USI& wId, OCPRes& res, const Bulk& bk, const OCP_DBL& dt) const
{
    if (opt.state == WellState::open) {

        const BulkVarSet& bvs = bk.vs;

        const USI len = nc + 2;

        // Well to Bulk
        for (USI p = 0; p < numPerf; p++) {
            const OCP_USI k = perf[p].location;
            // Mass Conservation
            for (USI i = 0; i < nc; i++) {
                res.resAbs[k * len + 1 + i] += perf[p].qi_lbmol[i] * dt;
            }
        }

        // Energy Conservation
        if (opt.type == WellType::injector) {
            const OCP_DBL Hw = bk.PVTm.GetPVT(0)->CalInjWellEnthalpy(opt.injTemp, &opt.injZi[0]);
            for (USI p = 0; p < numPerf; p++) {
                const OCP_USI k = perf[p].location;
                res.resAbs[k * len + 1 + nc] += perf[p].qt_ft3 * perf[p].xi * Hw * dt;
            }
        }
        else {
            for (USI p = 0; p < numPerf; p++) {
                const OCP_USI k = perf[p].location;
                for (USI j = 0; j < np; j++) {
                    res.resAbs[k * len + 1 + nc] += perf[p].qj_ft3[j] *
                        bvs.xi[k * np + j] * bvs.H[k * np + j] * dt;
                }
            }
        }

        switch (opt.mode)
        {
        case WellOptMode::bhp:
            res.resAbs[wId] = bhp - opt.tarBHP;
            break;
        case WellOptMode::irate:
        case WellOptMode::orate:
        case WellOptMode::grate:
        case WellOptMode::wrate:
        case WellOptMode::lrate:
            CalFactor(bk);
            if (opt.type == WellType::injector)  res.resAbs[wId] = opt.tarRate;
            else                                 res.resAbs[wId] = -opt.tarRate;
            for (USI i = 0; i < nc; i++) {
                res.resAbs[wId] += qi_lbmol[i] * factor[i];
            }
            res.maxWellRelRes_mol =
                max(res.maxWellRelRes_mol,
                    fabs(res.resAbs[wId] / opt.tarRate));
            break;
        default:
            OCP_ABORT("Wrong well opt mode!");
            break;
        }
        wId += len;
    }
}


void PeacemanWellT::AssembleMatFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{
    if (opt.state == WellState::open) {
        if (opt.type == WellType::injector)        AssembleMatInjFIM(ls, bk, dt);
        else if (opt.type == WellType::productor)  AssembleMatProdFIM(ls, bk, dt);
        else                                       OCP_ABORT("WRONG Well Type!");
    }
}


void PeacemanWellT::AssembleMatInjFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{

    const BulkVarSet& bvs = bk.vs;

    const USI ncol = nc + 2;
    const USI ncol2 = np * nc + np;
    const USI bsize = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
    OCP_USI   n_np_j;

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);

    OCP_DBL mu, muP, muT, dP, Hw;
    OCP_DBL transJ, transIJ;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, bmat);

    CalFactor(bk);

    Hw = bk.PVTm.GetPVT(0)->CalInjWellEnthalpy(opt.injTemp, &opt.injZi[0]);

    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        dP = bvs.P[n] - bhp - dG[p];

        for (USI j = 0; j < np; j++) {
            n_np_j = n * np + j;
            if (!bvs.phaseExist[n_np_j]) continue;

            mu = bvs.mu[n_np_j];
            muP = bvs.muP[n_np_j];
            muT = bvs.muT[n_np_j];

            for (USI i = 0; i < nc; i++) {

                // Mass Conservation
                if (!ifUseUnweight) {
                    transIJ = perf[p].transj[j] * perf[p].xi * opt.injZi[i];
                    // dQ / dP
                    dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu);
                    dQdXpW[(i + 1) * ncol] += -transIJ;

                    // dQ / dT
                    dQdXpB[(i + 2) * ncol - 1] += transIJ * (-dP * muT / mu);
                    dQdXpW[(i + 2) * ncol - 1] += 0;
                }
                else {
                    // dQ / dP
                    transIJ = perf[p].transINJ * perf[p].xi * opt.injZi[i];
                    dQdXpB[(i + 1) * ncol] += transIJ;
                    dQdXpW[(i + 1) * ncol] += -transIJ;
                }

                if (!ifUseUnweight) {
                    // dQ / dS
                    for (USI k = 0; k < np; k++) {
                        dQdXsB[(i + 1) * ncol2 + k] +=
                            CONV1 * perf[p].WI * perf[p].multiplier * perf[p].xi *
                            opt.injZi[i] * bvs.dKrdS[n_np_j * np + k] * dP / mu;
                    }
                    // dQ / dxij
                    for (USI k = 0; k < nc; k++) {
                        dQdXsB[(i + 1) * ncol2 + np + j * nc + k] +=
                            -transIJ * dP / mu * bvs.mux[n_np_j * nc + k];
                    }
                }
            }

            // Energy Conservation
            if (!ifUseUnweight) {
                transJ = perf[p].transj[j] * perf[p].xi;
                // dQ / dP
                dQdXpB[(nc + 1) * ncol] += transJ * Hw * (1 - dP * muP / mu);
                dQdXpW[(nc + 1) * ncol] += -transJ * Hw;

                // dQ / dT
                dQdXpB[(nc + 2) * ncol - 1] += transJ * Hw * (-dP * muT / mu);
                dQdXpW[(nc + 2) * ncol - 1] += 0;

                // dQ / dS
                for (USI k = 0; k < np; k++) {
                    dQdXsB[(nc + 1) * ncol2 + k] +=
                        CONV1 * perf[p].WI * perf[p].multiplier * perf[p].xi *
                        bvs.dKrdS[n_np_j * np + k] * dP / mu * Hw;
                }
                // dQ / dxij
                for (USI k = 0; k < nc; k++) {
                    dQdXsB[(nc + 1) * ncol2 + np + j * nc + k] +=
                        -transJ * dP / mu * bvs.mux[n_np_j * nc + k] * Hw;
                }
            }
            else {
                transJ = perf[p].transINJ * perf[p].xi;
                // dQ / dP
                dQdXpB[(nc + 1) * ncol] += transJ * Hw;
                dQdXpW[(nc + 1) * ncol] += -transJ * Hw;
            }

            if (ifUseUnweight) break;
        }

        // Bulk to Well
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bvs.dSec_dPri[n * bsize2], 1,
            bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Add
        ls.AddDiag(n, bmat);

        // Insert
        bmat = dQdXpW;
        Dscalar(bsize, dt, bmat.data());
        ls.NewOffDiag(n, wId, bmat);

        // Well
        switch (opt.mode) {
        case WellOptMode::irate:
        case WellOptMode::orate:
        case WellOptMode::grate:
        case WellOptMode::wrate:
        case WellOptMode::lrate:
            // Diag
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < nc; i++) {
                bmat[0] += dQdXpW[(i + 1) * ncol] * factor[i];
                bmat[(i + 1) * ncol + i + 1] = 1;
            }
            bmat[ncol * ncol - 1] = 1;
            ls.AddDiag(wId, bmat);

            // OffDiag
            bmat = dQdXpB;
            DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bvs.dSec_dPri[n * bsize2],
                1, bmat.data());
            fill(bmat2.begin(), bmat2.end(), 0.0);
            for (USI i = 0; i < nc; i++) {
                Daxpy(ncol, factor[i], bmat.data() + (i + 1) * ncol, bmat2.data());
            }
            ls.NewOffDiag(wId, n, bmat2);
            break;

        case WellOptMode::bhp:
            // Diag
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < ncol; i++) {
                bmat[i * ncol + i] = 1;
            }
            // Add
            ls.AddDiag(wId, bmat);
            // OffDiag
            fill(bmat.begin(), bmat.end(), 0.0);
            // Insert
            ls.NewOffDiag(wId, n, bmat);
            break;

        default:
            OCP_ABORT("Wrong Well Opt mode!");
            break;
        }
    }
}


void PeacemanWellT::AssembleMatProdFIM(LinearSystem& ls, const Bulk& bk, const OCP_DBL& dt) const
{
    const BulkVarSet& bvs = bk.vs;

    const USI ncol = nc + 2;
    const USI ncol2 = np * nc + np;
    const USI bsize = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
    OCP_USI   n_np_j;

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);

    OCP_DBL xij, xi, xiP, xiT, mu, muP, muT, dP, transIJ, transJ, H, HT, Hx, tmp;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, bmat);

    // Set Prod Weight
    CalFactor(bk);

    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        for (USI j = 0; j < np; j++) {
            n_np_j = n * np + j;
            if (!bvs.phaseExist[n_np_j]) continue;

            dP = bvs.Pj[n_np_j] - bhp - dG[p];
            xi = bvs.xi[n_np_j];
            mu = bvs.mu[n_np_j];
            xiP = bvs.xiP[n_np_j];
            xiT = bvs.xiT[n_np_j];
            muP = bvs.muP[n_np_j];
            muT = bvs.muT[n_np_j];
            H = bvs.H[n_np_j];
            HT = bvs.HT[n_np_j];

            // Mass Conservation
            for (USI i = 0; i < nc; i++) {
                xij = bvs.xij[n_np_j * nc + i];
                Hx = bvs.Hx[n_np_j * nc + i];
                // dQ / dP
                transIJ = perf[p].transj[j] * xi * xij;
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu) +
                    dP * perf[p].transj[j] * xij * xiP;
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dT
                dQdXpB[(i + 2) * ncol - 1] +=
                    transIJ * (-dP * muT / mu) + dP * perf[p].transj[j] * xij * xiT;
                dQdXpW[(i + 2) * ncol - 1] += 0;

                // dQ / dS
                for (USI k = 0; k < np; k++) {
                    tmp = CONV1 * perf[p].WI * perf[p].multiplier * dP / mu * xi *
                        xij * bvs.dKrdS[n_np_j * np + k];
                    // capillary pressure
                    tmp += transIJ * bvs.dPcdS[n_np_j * np + k];
                    dQdXsB[(i + 1) * ncol2 + k] += tmp;
                }
                // dQ / dxij
                for (USI k = 0; k < nc; k++) {
                    tmp = dP * perf[p].transj[j] * xij *
                        (bvs.xix[n_np_j * nc + k] - xi / mu * bvs.mux[n_np_j * nc + k]);
                    dQdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                }
                dQdXsB[(i + 1) * ncol2 + np + j * nc + i] +=
                    perf[p].transj[j] * xi * dP;
            }

            // Energy Conservation
            transJ = perf[p].transj[j] * xi;
            // dQ / dP
            dQdXpB[(nc + 1) * ncol] +=
                transJ * (1 - dP * muP / mu) * H + dP * perf[p].transj[j] * xiP * H;
            dQdXpW[(nc + 1) * ncol] += -transJ * H;

            // dQ / dT
            dQdXpB[(nc + 2) * ncol - 1] += transJ * (-dP * muT / mu) * H +
                dP * perf[p].transj[j] * xiT * H +
                transJ * dP * HT;
            dQdXpW[(nc + 2) * ncol - 1] += 0;

            // dQ / dS
            for (USI k = 0; k < np; k++) {
                tmp = CONV1 * perf[p].WI * perf[p].multiplier * dP / mu * xi *
                    bvs.dKrdS[n_np_j * np + k] * H;
                // capillary pressure
                tmp += transJ * bvs.dPcdS[n_np_j * np + k] * H;
                dQdXsB[(nc + 1) * ncol2 + k] += tmp;
            }

            // dQ / dxij
            for (USI k = 0; k < nc; k++) {
                tmp = dP * perf[p].transj[j] *
                    (bvs.xix[n_np_j * nc + k] - xi / mu * bvs.mux[n_np_j * nc + k]) *
                    H +
                    transJ * dP * Hx;
                dQdXsB[(nc + 1) * ncol2 + np + j * nc + k] += tmp;
            }
        }

        // Bulk to Well
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bvs.dSec_dPri[n * bsize2], 1,
            bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Add
        ls.AddDiag(n, bmat);
        // Insert
        bmat = dQdXpW;
        Dscalar(bsize, dt, bmat.data());
        ls.NewOffDiag(n, wId, bmat);

        // Well
        switch (opt.mode) {
        case WellOptMode::irate:
        case WellOptMode::orate:
        case WellOptMode::grate:
        case WellOptMode::wrate:
        case WellOptMode::lrate:
            // Diag
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < nc; i++) {
                bmat[0] += dQdXpW[(i + 1) * ncol] * factor[i];
                bmat[(i + 1) * ncol + i + 1] = 1;
            }
            bmat[ncol * ncol - 1] = 1;
            ls.AddDiag(wId, bmat);

            // OffDiag
            bmat = dQdXpB;
            DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bvs.dSec_dPri[n * bsize2],
                1, bmat.data());
            fill(bmat2.begin(), bmat2.end(), 0.0);
            for (USI i = 0; i < nc; i++) {
                Daxpy(ncol, factor[i], bmat.data() + (i + 1) * ncol,
                    bmat2.data());
            }
            ls.NewOffDiag(wId, n, bmat2);
            break;

        case WellOptMode::bhp:
            // Diag
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < ncol; i++) {
                bmat[i * ncol + i] = 1;
            }
            // Add
            ls.AddDiag(wId, bmat);
            // OffDiag
            fill(bmat.begin(), bmat.end(), 0.0);
            // Insert
            ls.NewOffDiag(wId, n, bmat);
            break;

        default:
            OCP_ABORT("Wrong Well Opt mode!");
            break;
        }
    }
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/17/2023      Create file                          */
/*----------------------------------------------------------------------------*/