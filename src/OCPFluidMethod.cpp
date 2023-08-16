/*! \file    OCPFluidMethod.cpp
 *  \brief   Definition of solution methods for fluid part in OpenCAEPoroX
 *  \author  Shizhe Li
 *  \date    Nov/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPFluidMethod.hpp"

////////////////////////////////////////////
// IsothermalMethod
////////////////////////////////////////////


void IsothermalMethod::CalRock(Bulk& bk) const
{
    for (OCP_USI n = 0; n < bk.numBulk; n++) {
        bk.rock[bk.ROCKNUM[n]]->CalPoro(bk.P[n], bk.T[n], bk.poroInit[n], 0);
        bk.poro[n]   = bk.rock[bk.ROCKNUM[n]]->GetPoro();
        bk.poroP[n]  = bk.rock[bk.ROCKNUM[n]]->GetdPorodP();
        bk.rockVp[n] = bk.v[n] * bk.poro[n];
    }
}

////////////////////////////////////////////
// IsoT_IMPEC
////////////////////////////////////////////

void IsoT_IMPEC::Setup(Reservoir& rs, LinearSystem& ls, const OCPControl& ctrl)
{
    // Allocate Memory of auxiliary variables for IMPEC
    AllocateReservoir(rs);
    // Allocate Memory of Matrix for IMPEC
    AllocateLinearSystem(ls, rs, ctrl);
}

/// Initialize reservoir
void IsoT_IMPEC::InitReservoir(Reservoir& rs) const
{
    rs.bulk.InitPTSw(50);

    CalRock(rs.bulk);

    InitFlash(rs.bulk);
    CalKrPc(rs.bulk);

    CalBulkFlux(rs);

    rs.allWells.InitBHP(rs.bulk);

    UpdateLastTimeStep(rs);
}

void IsoT_IMPEC::Prepare(Reservoir& rs, OCPControl& ctrl)
{
    rs.allWells.PrepareWell(rs.bulk);
    rs.CalCFL(ctrl.GetCurDt(), OCP_TRUE);
    ctrl.Check(rs, {"CFL"});
}

void IsoT_IMPEC::AssembleMat(LinearSystem&    ls,
                             const Reservoir& rs,
                             const OCP_DBL&   dt) const
{
    AssembleMatBulks(ls, rs, dt);
    AssembleMatWells(ls, rs, dt);
}

void IsoT_IMPEC::SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl)
{
#ifdef DEBUG
    ls.CheckEquation();
#endif // DEBUG

    GetWallTime timer;

    timer.Start();
    ls.CalCommTerm(rs.GetNumOpenWell());
    ls.AssembleMatLinearSolver();
    OCPTIME_ASSEMBLE_MAT_FOR_LS += timer.Stop() / 1000;
    
    timer.Start();
    int status = ls.Solve();
    if (status < 0) {
        status = ls.GetNumIters();
    }

#ifdef DEBUG
    //OCP_INT myrank = rs.domain.myrank;
    //ls.OutputLinearSystem("proc" + to_string(myrank) + "_testA_IMPEC.out", 
    //                      "proc" + to_string(myrank) + "_testb_IMPEC.out");
    //ls.OutputSolution("proc" + to_string(myrank) + "_testx_IMPEC.out");
    //MPI_Barrier(rs.domain.myComm);
    //OCP_ABORT("Stop");
#endif // DEBUG

    OCPTIME_LSOLVER += timer.Stop() / 1000;
    ctrl.UpdateIterLS(status);
    ctrl.UpdateIterNR();

#ifdef DEBUG
    // ls.OutputSolution("testx_IMPEC.out");
#endif // DEBUG

    timer.Start();
    GetSolution(rs, ls.GetSolution());
    OCPTIME_NRSTEP += timer.Stop() / 1000;
    ls.ClearData();
}

OCP_BOOL IsoT_IMPEC::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    // First check : Pressure check
    if (!ctrl.Check(rs, {"BulkP", "WellP"})) {
        return OCP_FALSE;
    }

    // Calculate Flux between bulks and between bulks and wells
    CalFlux(rs);  
    MassConserve(rs, dt);

    // Second check : CFL check
    rs.CalCFL(dt, OCP_TRUE);
    // Third check: Ni check
    if (!ctrl.Check(rs, { "CFL","BulkNi"})) {
        ResetToLastTimeStep01(rs, ctrl);
        return OCP_FALSE;
    }

    CalRock(rs.bulk);
    CalFlash(rs.bulk);

    // Fouth check: Volume error check
    if (!ctrl.Check(rs, {"BulkVe"})) {
        ResetToLastTimeStep02(rs, ctrl);
        return OCP_FALSE;
    }

    CalKrPc(rs.bulk);
    CalBulkFlux(rs);

    return OCP_TRUE;
}

OCP_BOOL IsoT_IMPEC::FinishNR(const Reservoir& rs) { return OCP_TRUE; }

void IsoT_IMPEC::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.CalIPRT(ctrl.GetCurDt());
    rs.CalMaxChange();
    UpdateLastTimeStep(rs);
    // ctrl.CalNextTstepIMPEC(rs);
    ctrl.CalNextTimeStep(rs, {"dP", "dN", "dS", "eV"});
}

void IsoT_IMPEC::AllocateReservoir(Reservoir& rs)
{
    Bulk&         bk = rs.bulk;
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;
    const USI     nc = bk.numCom;

    // Rock
    bk.poro.resize(nb);
    bk.rockVp.resize(nb);

    bk.lporo.resize(nb);
    bk.lrockVp.resize(nb);

    // derivatives
    bk.poroP.resize(nb);
    bk.lporoP.resize(nb);

    // Fluid
    bk.phaseNum.resize(nb);
    bk.Nt.resize(nb);
    bk.Ni.resize(nb * nc);
    bk.vf.resize(nb);
    bk.T.resize(nb);
    bk.P.resize(nb);
    bk.Pb.resize(nb);
    bk.Pj.resize(nb * np);
    bk.Pc.resize(nb * np);
    bk.phaseExist.resize(nb * np);
    bk.S.resize(nb * np);
    bk.vj.resize(nb * np);
    bk.xij.resize(nb * np * nc);
    bk.rho.resize(nb * np);
    bk.xi.resize(nb * np);
    bk.mu.resize(nb * np);
    bk.kr.resize(nb * np);

    bk.lphaseNum.resize(nb);
    bk.lNt.resize(nb);
    bk.lNi.resize(nb * nc);
    bk.lvf.resize(nb);
    bk.lT.resize(nb);
    bk.lP.resize(nb);
    bk.lPj.resize(nb * np);
    bk.lPc.resize(nb * np);
    bk.lphaseExist.resize(nb * np);
    bk.lS.resize(nb * np);
    bk.vj.resize(nb * np);
    bk.lxij.resize(nb * np * nc);
    bk.lrho.resize(nb * np);
    bk.lxi.resize(nb * np);
    bk.lmu.resize(nb * np);
    bk.lkr.resize(nb * np);

    // derivatives
    bk.vfP.resize(nb);
    bk.vfi.resize(nb * nc);

    bk.lvfP.resize(nb);
    bk.lvfi.resize(nb * nc);

    // others
    bk.cfl.resize(nb * np);

    BulkConn& conn = rs.conn;

    conn.bcval.upblock.resize(conn.numConn * np);
    conn.bcval.rho.resize(conn.numConn * np);
    conn.bcval.velocity.resize(conn.numConn * np);
    conn.bcval.flux_ni.resize(conn.numConn * nc);

    conn.bcval.lupblock.resize(conn.numConn * np);
    conn.bcval.lrho.resize(conn.numConn * np);
    conn.bcval.lvelocity.resize(conn.numConn * np);
    conn.bcval.lflux_ni.resize(conn.numConn * nc);
}

void IsoT_IMPEC::AllocateLinearSystem(LinearSystem&     ls,
                                      const Reservoir&  rs,
                                      const OCPControl& ctrl)
{
    ls.SetupDomain(rs.domain);
    ls.AllocateRowMem(1);
    ls.AllocateColMem();
    ls.SetupLinearSolver(ISOTHERMALMODEL, ctrl.GetWorkDir(), ctrl.GetLsFile());
}

void IsoT_IMPEC::InitFlash(Bulk& bk) const
{
    for (OCP_USI n = 0; n < bk.numBulk; n++) {
        bk.flashCal[bk.PVTNUM[n]]->InitFlashIMPEC(bk.P[n], bk.Pb[n], bk.T[n],
                                                  &bk.S[n * bk.numPhase], bk.rockVp[n],
                                                  bk.Ni.data() + n * bk.numCom, n);
        for (USI i = 0; i < bk.numCom; i++) {
            bk.Ni[n * bk.numCom + i] = bk.flashCal[bk.PVTNUM[n]]->GetNi(i);
        }
        PassFlashValue(bk, n);
    }
}

void IsoT_IMPEC::CalFlash(Bulk& bk)
{
    for (OCP_USI n = 0; n < bk.numBulk; n++) {

        bk.flashCal[bk.PVTNUM[n]]->FlashIMPEC(bk.P[n], bk.T[n], &bk.Ni[n * bk.numCom],
                                              bk.phaseNum[n],
                                              &bk.xij[n * bk.numPhase * bk.numCom], n);
        PassFlashValue(bk, n);
    }
}

void IsoT_IMPEC::PassFlashValue(Bulk& bk, const OCP_USI& n) const
{
    const USI     np     = bk.numPhase;
    const USI     nc     = bk.numCom;
    const OCP_USI bIdp   = n * np;
    const USI     pvtnum = bk.PVTNUM[n];

    bk.phaseNum[n] = 0;
    bk.Nt[n]       = bk.flashCal[pvtnum]->GetNt();
    bk.vf[n]       = bk.flashCal[pvtnum]->GetVf();

    for (USI j = 0; j < np; j++) {
        // Important! Saturation must be passed no matter if the phase exists. This is
        // because it will be used to calculate relative permeability and capillary
        // pressure at each time step. Make sure that all saturations are updated at
        // each step!
        bk.phaseExist[bIdp + j] = bk.flashCal[pvtnum]->GetPhaseExist(j);
        bk.S[bIdp + j]          = bk.flashCal[pvtnum]->GetS(j);
        if (bk.phaseExist[bIdp + j]) {
            bk.phaseNum[n]++;
            for (USI i = 0; i < nc; i++) {
                bk.xij[bIdp * nc + j * nc + i] = bk.flashCal[pvtnum]->GetXij(j, i);
            }
            bk.vj[bIdp + j]  = bk.flashCal[pvtnum]->GetVj(j);
            bk.rho[bIdp + j] = bk.flashCal[pvtnum]->GetRho(j);
            bk.xi[bIdp + j]  = bk.flashCal[pvtnum]->GetXi(j);
            bk.mu[bIdp + j]  = bk.flashCal[pvtnum]->GetMu(j);
        }
    }

    bk.vfP[n] = bk.flashCal[pvtnum]->GetVfP();
    for (USI i = 0; i < nc; i++) {
        bk.vfi[n * nc + i] = bk.flashCal[pvtnum]->GetVfi(i);
    }
}

void IsoT_IMPEC::CalKrPc(Bulk& bk) const
{
    for (OCP_USI n = 0; n < bk.numBulk; n++) {
        OCP_USI bId = n * bk.numPhase;
        bk.flow[bk.SATNUM[n]]->CalKrPc(&bk.S[bId], n);
        copy(bk.flow[bk.SATNUM[n]]->GetKr().begin(), bk.flow[bk.SATNUM[n]]->GetKr().end(), &bk.kr[bId]);
        copy(bk.flow[bk.SATNUM[n]]->GetPc().begin(), bk.flow[bk.SATNUM[n]]->GetPc().end(), &bk.Pc[bId]);
        for (USI j = 0; j < bk.numPhase; j++)
            bk.Pj[n * bk.numPhase + j] = bk.P[n] + bk.Pc[n * bk.numPhase + j];
    }
}

void IsoT_IMPEC::CalFlux(Reservoir& rs) const
{
    CalBulkFlux(rs);
    rs.allWells.CalFlux(rs.bulk);
}

void IsoT_IMPEC::CalBulkFlux(Reservoir& rs) const
{
    const Bulk& bk   = rs.bulk;
    BulkConn&   conn = rs.conn;
    const USI   np   = bk.numPhase;
    const USI   nc   = bk.numCom;

    // calculate a step flux using iteratorConn

    vector<OCPFlux*>& flux  = conn.flux;
    BulkConnVal&      bcval = conn.bcval;

    for (OCP_USI c = 0; c < conn.numConn; c++) {

        const USI cType = conn.iteratorConn[c].Type();
        flux[cType]->CalFlux(conn.iteratorConn[c], bk);
        copy(flux[cType]->GetUpblock().begin(), flux[cType]->GetUpblock().end(), &bcval.upblock[c * np]);
        copy(flux[cType]->GetRho().begin(), flux[cType]->GetRho().end(), &bcval.rho[c * np]);
        copy(flux[cType]->GetFluxVj().begin(), flux[cType]->GetFluxVj().end(), &bcval.velocity[c * np]);
        copy(flux[cType]->GetFluxNi().begin(), flux[cType]->GetFluxNi().end(), &bcval.flux_ni[c * nc]);
    }
}

void IsoT_IMPEC::MassConserve(Reservoir& rs, const OCP_DBL& dt) const
{

    // Bulk to Bulk
    Bulk&           bk   = rs.bulk;
    const USI       nc   = bk.numCom;
    const BulkConn& conn = rs.conn;
    
    OCP_USI bId, eId;

    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId = conn.iteratorConn[c].BId();
        eId = conn.iteratorConn[c].EId();

        for (USI i = 0; i < nc; i++) {
            bk.Ni[eId * nc + i] += dt * conn.bcval.flux_ni[c * nc + i];
            bk.Ni[bId * nc + i] -= dt * conn.bcval.flux_ni[c * nc + i];
        }
    }

    // Well to Bulk
    for (auto& wl : rs.allWells.wells) {
        if (wl.IsOpen()) {
            for (USI p = 0; p < wl.PerfNum(); p++) {
                OCP_USI k = wl.PerfLocation(p);
                for (USI i = 0; i < nc; i++) {
                    bk.Ni[k * nc + i] -= wl.PerfQi_lbmol(p, i) * dt;
                }
            }
        }
    }

    // Exchange Ghost Ni
    const Domain& domain = rs.domain;

    for (USI i = 0; i < domain.numRecvProc; i++) {
        const vector<OCP_USI>& rel = domain.recv_element_loc[i];
        MPI_Irecv(&bk.Ni[rel[1] * nc], (rel[2] - rel[1]) * nc, MPI_DOUBLE, rel[0], 0, domain.myComm, &domain.recv_request[i]);
    }

    vector<vector<OCP_DBL>> send_buffer(domain.numSendProc);
    for (USI i = 0; i < domain.numSendProc; i++) {
        const vector<OCP_USI>& sel = domain.send_element_loc[i];
        vector<OCP_DBL>&       s   = send_buffer[i];
        s.resize(1 + (sel.size() - 1) * nc);
        s[0] = sel[0];
        for (USI j = 1; j < sel.size(); j++) {
            const OCP_DBL* bId = &bk.Ni[0] + sel[j] * nc;
            copy(bId, bId + nc, &s[1 + (j - 1) * nc]);
        }
        MPI_Isend(s.data() + 1, s.size() - 1, MPI_DOUBLE, s[0], 0, domain.myComm, &domain.send_request[i]);
    }

    MPI_Waitall(domain.numSendProc, domain.send_request.data(), MPI_STATUS_IGNORE);
    MPI_Waitall(domain.numRecvProc, domain.recv_request.data(), MPI_STATUS_IGNORE);
}

void IsoT_IMPEC::AssembleMatBulks(LinearSystem&    ls,
                                  const Reservoir& rs,
                                  const OCP_DBL&   dt) const
{
    const USI numWell = rs.GetNumOpenWell();

    const Bulk&     bk   = rs.bulk;
    const BulkConn& conn = rs.conn;
    const OCP_USI   nb   = bk.numBulkInterior;

    ls.AddDim(nb);

    // accumulate term
    OCP_DBL Vpp, Vp, vf, vfP, P;
    for (OCP_USI n = 0; n < nb; n++) {
        vf  = bk.vf[n];
        vfP = bk.vfP[n];
        P   = bk.lP[n];
        Vpp = bk.v[n] * bk.poroP[n];
        Vp  = bk.rockVp[n];

        ls.NewDiag(n, Vpp - vfP);
        ls.AddRhs(n, (Vpp - vfP) * P + dt * (vf - Vp));
    }


    // flux term
    OCP_USI bId, eId;
    USI     cType;
    OCP_DBL valbb, rhsb, valee, rhse;

    // Be careful when first bulk has no neighbors!
    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId   = conn.iteratorConn[c].BId();
        eId   = conn.iteratorConn[c].EId();
        cType = conn.iteratorConn[c].Type();

        conn.flux[cType]->AssembleMatIMPEC(conn.iteratorConn[c], c, conn.bcval, bk);
        valbb  = dt * conn.flux[cType]->GetValbb();
        valee  = dt * conn.flux[cType]->GetValee();
        rhsb   = dt * conn.flux[cType]->GetRhsb();
        rhse   = dt * conn.flux[cType]->GetRhse();


        if (eId < nb) {
            // interior grid
            ls.AddDiag(bId, valbb);
            ls.AddDiag(eId, valee);
            ls.NewOffDiag(bId, eId, -valbb);
            ls.NewOffDiag(eId, bId, -valee);
            ls.AddRhs(bId, rhsb);
            ls.AddRhs(eId, rhse);
        }
        else {
            // ghost grid
            ls.AddDiag(bId, valbb);
            ls.NewOffDiag(bId, eId + numWell, -valbb);
            ls.AddRhs(bId, rhsb);
        }
    }
}

void IsoT_IMPEC::AssembleMatWells(LinearSystem&    ls,
                                  const Reservoir& rs,
                                  const OCP_DBL&   dt) const
{
    for (auto& wl : rs.allWells.wells) {
        wl.AssembleMatIMPEC(ls, rs.bulk, dt);
    }

    // for Reinjection
    // for (auto& wG : wellGroup) {
    //    if (wG.reInj) {
    //        for (auto& prod : wellGroup[wG.prodGroup].wIdPROD) {
    //            if (wells[prod].IsOpen()) {
    //                wells[prod].AssembleMatReinjection_IMPEC(myBulk, myLS, dt, wells,
    //                    wG.wIdINJ);
    //            }
    //        }
    //    }
    //}
}


void IsoT_IMPEC::GetSolution(Reservoir& rs, vector<OCP_DBL>& u)
{
    Bulk&         bk     = rs.bulk;
    const OCP_USI nb     = bk.numBulk;
    const USI     np     = bk.numPhase;
    const Domain& domain = rs.domain;

    // Well first
    USI wId = bk.numBulkInterior;
    for (auto& wl : rs.allWells.wells) {
        if (wl.IsOpen()) {
            wl.SetBHP(u[wId]);
            wl.CalPerfP();
            wId++;
        }
    }

    // Exchange Solution

    for (USI i = 0; i < domain.numRecvProc; i++) {
        const vector<OCP_USI>& rel = domain.recv_element_loc[i];
        MPI_Irecv(&u[rel[1]], rel[2] - rel[1], MPI_DOUBLE, rel[0], 0, domain.myComm, &domain.recv_request[i]);
    }

    vector<vector<OCP_DBL>> send_buffer(domain.numSendProc);
    for (USI i = 0; i < domain.numSendProc; i++) {
        const vector<OCP_USI>& sel = domain.send_element_loc[i];
        vector<OCP_DBL>&       s   = send_buffer[i];
        s.resize(sel.size());
        s[0] = sel[0];
        for (USI j = 1; j < sel.size(); j++) {
            s[j] = u[sel[j]];
        }
        MPI_Isend(s.data() + 1, s.size() - 1, MPI_DOUBLE, s[0], 0, domain.myComm, &domain.send_request[i]);
    }

    // Bulk
    // interior first, ghost second
    OCP_USI bId = 0;
    OCP_USI eId = bk.GetInteriorBulkNum();
    for (USI p = 0; p < 2; p++) {

        for (OCP_USI n = bId; n < eId; n++) {
            bk.P[n] = u[n];
            for (USI j = 0; j < np; j++) {
                bk.Pj[n * np + j] = bk.P[n] + bk.Pc[n * np + j];
            }
        }

        if (p == 0) {
            bId = eId;
            eId = nb;
            MPI_Waitall(domain.numRecvProc, domain.recv_request.data(), MPI_STATUS_IGNORE);
        }
        else {
            break;
        }
    }
    MPI_Waitall(domain.numSendProc, domain.send_request.data(), MPI_STATUS_IGNORE);
}


void IsoT_IMPEC::ResetToLastTimeStep01(Reservoir& rs, OCPControl& ctrl)
{
    // Bulk
    rs.bulk.Ni = rs.bulk.lNi;
    rs.bulk.Pj = rs.bulk.lPj;
    // Bulk Conn

    rs.conn.bcval.upblock      = rs.conn.bcval.lupblock;
    rs.conn.bcval.rho          = rs.conn.bcval.lrho;
    rs.conn.bcval.velocity     = rs.conn.bcval.lvelocity;
    rs.conn.bcval.flux_ni      = rs.conn.bcval.lflux_ni;

    // Iters
    ctrl.ResetIterNRLS();
}

void IsoT_IMPEC::ResetToLastTimeStep02(Reservoir& rs, OCPControl& ctrl)
{
    Bulk& bk = rs.bulk;
    // Rock
    bk.rockVp = bk.lrockVp;
    bk.poro   = bk.lporo;
    bk.poroP  = bk.lporoP;

    // Fluid
    bk.phaseNum   = bk.lphaseNum;
    bk.Nt         = bk.lNt;
    bk.Ni         = bk.lNi;
    bk.vf         = bk.lvf;
    bk.Pj         = bk.lPj;
    bk.phaseExist = bk.lphaseExist;
    bk.S          = bk.lS;
    bk.vj         = bk.lvj;
    bk.xij        = bk.lxij;
    bk.rho        = bk.lrho;
    bk.xi         = bk.lxi;
    bk.mu         = bk.lmu;

    // derivatives
    bk.vfP = bk.lvfP;
    bk.vfi = bk.lvfi;

    // Bulk Conn
    rs.conn.bcval.upblock      = rs.conn.bcval.lupblock;
    rs.conn.bcval.rho          = rs.conn.bcval.lrho;
    rs.conn.bcval.velocity     = rs.conn.bcval.lvelocity;
    rs.conn.bcval.flux_ni      = rs.conn.bcval.lflux_ni;

    // Optional Features
    rs.optFeatures.ResetToLastTimeStep();

    // Iters
    ctrl.ResetIterNRLS();
}

void IsoT_IMPEC::UpdateLastTimeStep(Reservoir& rs) const
{

    Bulk& bk = rs.bulk;

    // Rock
    bk.lporo   = bk.poro;
    bk.lporoP  = bk.poroP;
    bk.lrockVp = bk.rockVp;

    // Fluid
    bk.lphaseNum   = bk.phaseNum;
    bk.lNt         = bk.Nt;
    bk.lNi         = bk.Ni;
    bk.lvf         = bk.vf;
    bk.lP          = bk.P;
    bk.lPj         = bk.Pj;
    bk.lPc         = bk.Pc;
    bk.lphaseExist = bk.phaseExist;
    bk.lS          = bk.S;
    bk.lvj         = bk.vj;
    bk.lxij        = bk.xij;
    bk.lrho        = bk.rho;
    bk.lxi         = bk.xi;
    bk.lmu         = bk.mu;
    bk.lkr         = bk.kr;

    // derivatives
    bk.lvfP = bk.vfP;
    bk.lvfi = bk.vfi;

    BulkConn& conn = rs.conn;

    conn.bcval.lupblock    = conn.bcval.upblock;
    conn.bcval.lrho        = conn.bcval.rho;
    conn.bcval.lvelocity   = conn.bcval.velocity;
    conn.bcval.lflux_ni    = conn.bcval.flux_ni;

    rs.allWells.UpdateLastTimeStepBHP();
    rs.optFeatures.UpdateLastTimeStep();
}

////////////////////////////////////////////
// IsoT_FIM
////////////////////////////////////////////

void IsoT_FIM::Setup(Reservoir& rs, LinearSystem& ls, const OCPControl& ctrl)
{
    // Allocate memory for reservoir
    AllocateReservoir(rs);
    // Allocate memory for linear system
    AllocateLinearSystem(ls, rs, ctrl);
}

void IsoT_FIM::InitReservoir(Reservoir& rs) const
{
    // Calculate initial bulk pressure and temperature and water saturation
    rs.bulk.InitPTSw(50);
    // Initialize rock property
    CalRock(rs.bulk);
    // Initialize fluid properties
    InitFlash(rs.bulk);
    CalKrPc(rs.bulk);
    // Initialize well pressure
    rs.allWells.InitBHP(rs.bulk);
    // Update variables at last time step
    UpdateLastTimeStep(rs);
}

void IsoT_FIM::Prepare(Reservoir& rs, const OCP_DBL& dt)
{
    // Calculate well property at the beginning of next time step
    rs.allWells.PrepareWell(rs.bulk);
    // Calculate initial residual
    CalRes(rs, dt, OCP_TRUE);
}

void IsoT_FIM::AssembleMat(LinearSystem&    ls,
                           const Reservoir& rs,
                           const OCP_DBL&   dt) const
{
    // Assemble matrix
    AssembleMatBulks(ls, rs, dt);
    AssembleMatWells(ls, rs, dt);
    // Assemble rhs -- from residual
    ls.AssembleRhsCopy(rs.bulk.res.resAbs);
}

void IsoT_FIM::SolveLinearSystem(LinearSystem& ls,
                                 Reservoir&    rs,
                                 OCPControl&   ctrl) const
{
#ifdef DEBUG
    // Check if inf or nan occurs in A and b
    ls.CheckEquation();
#endif // DEBUG

    GetWallTime timer;
    timer.Start();
    ls.CalCommTerm(rs.GetNumOpenWell());
    ls.AssembleMatLinearSolver();
    OCPTIME_ASSEMBLE_MAT_FOR_LS += timer.Stop() / 1000;

    // Solve linear system
    
    timer.Start();
    int status = ls.Solve();
    if (status < 0) {
        status = ls.GetNumIters();
    }
    // Record time, iterations
    OCPTIME_LSOLVER += timer.Stop() / 1000;
    ctrl.UpdateIterLS(status);
    ctrl.UpdateIterNR();

     
#ifdef DEBUG
    // Output A, b, x
    
    //ls.OutputLinearSystem("proc" + to_string(CURRENT_RANK) + "_testA_FIM.out",
    //                      "proc" + to_string(CURRENT_RANK) + "_testb_FIM.out");
    //MPI_Barrier(rs.domain.myComm);
    //OCP_ABORT("Stop");
    //
    //ls.OutputSolution("proc" + to_string(CURRENT_RANK) + "_testx_FIM.out");
    // Check if inf or nan occurs in solution
    // ls.CheckSolution();
#endif // DEBUG

    // Get solution from linear system to Reservoir
    timer.Start();
    GetSolution(rs, ls.GetSolution(), ctrl);
    OCPTIME_NRSTEP += timer.Stop() / 1000;   
    // rs.PrintSolFIM(ctrl.workDir + "testPNi.out");
    ls.ClearData();
}

OCP_BOOL IsoT_FIM::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    if (!ctrl.Check(rs, {"BulkNi", "BulkP"})) {
        ResetToLastTimeStep(rs, ctrl);
        cout << "Cut time step size and repeat! current dt = " << fixed
             << setprecision(3) << dt << " days\n";
        return OCP_FALSE;
    }

    // Update fluid property
    CalFlash(rs.bulk);
    CalKrPc(rs.bulk);
    // Update rock property
    CalRock(rs.bulk);
    // Update well property
    rs.allWells.CalTrans(rs.bulk);
    rs.allWells.CalFlux(rs.bulk);
    // Update residual
    CalRes(rs, dt, OCP_FALSE);

    return OCP_TRUE;
}

OCP_BOOL IsoT_FIM::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    OCP_USI dSn;

    const OCP_DBL NRdSmax = rs.GetNRdSmax(dSn);
    const OCP_DBL NRdPmax = rs.GetNRdPmax();
    // const OCP_DBL NRdNmax = rs.GetNRdNmax();
   
    OCP_INT conflag_loc = -1;
    if (((rs.bulk.res.maxRelRes_V <= rs.bulk.res.maxRelRes0_V * ctrl.ctrlNR.NRtol ||
        rs.bulk.res.maxRelRes_V <= ctrl.ctrlNR.NRtol ||
        rs.bulk.res.maxRelRes_N <= ctrl.ctrlNR.NRtol) &&
        rs.bulk.res.maxWellRelRes_mol <= ctrl.ctrlNR.NRtol) ||
        (fabs(NRdPmax) <= ctrl.ctrlNR.NRdPmin &&
            fabs(NRdSmax) <= ctrl.ctrlNR.NRdSmin)) {
        conflag_loc = 0;
    }

    GetWallTime timer;
    timer.Start();

    OCP_INT conflag;
    MPI_Allreduce(&conflag_loc, &conflag, 1, MPI_INT, MPI_MIN, rs.domain.myComm);

    OCPTIME_COMM_COLLECTIVE += timer.Stop() / 1000;


    if (conflag == 0) {

        if (!ctrl.Check(rs, {"WellP"})) {
            ResetToLastTimeStep(rs, ctrl);
            return OCP_FALSE;
        } else {
            return OCP_TRUE;
        }
    } else if (ctrl.iterNR >= ctrl.ctrlNR.maxNRiter) {
        ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
        ResetToLastTimeStep(rs, ctrl);
        cout << "### WARNING: NR not fully converged! Cut time step size and repeat!  "
                "current dt = "
             << fixed << setprecision(3) << ctrl.current_dt << " days\n";
        return OCP_FALSE;
    } else {
        return OCP_FALSE;
    }
}

void IsoT_FIM::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.CalIPRT(ctrl.GetCurDt());
    rs.CalMaxChange();
    UpdateLastTimeStep(rs);
    ctrl.CalNextTimeStep(rs, {"dP", "dS", "iter"});
}

void IsoT_FIM::AllocateReservoir(Reservoir& rs)
{
    Bulk&         bk = rs.bulk;
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;
    const USI     nc = bk.numCom;

    // Rock
    bk.poro.resize(nb);
    bk.rockVp.resize(nb);

    bk.lporo.resize(nb);
    bk.lrockVp.resize(nb);

    // derivatives
    bk.poroP.resize(nb);
    bk.lporoP.resize(nb);

    // Fluid
    bk.phaseNum.resize(nb);
    bk.Nt.resize(nb);
    bk.Ni.resize(nb * nc);
    bk.vf.resize(nb);
    bk.T.resize(nb);
    bk.P.resize(nb);
    bk.Pb.resize(nb);
    bk.Pj.resize(nb * np);
    bk.Pc.resize(nb * np);
    bk.phaseExist.resize(nb * np);
    bk.S.resize(nb * np);
    bk.xij.resize(nb * np * nc);
    bk.rho.resize(nb * np);
    bk.xi.resize(nb * np);
    bk.mu.resize(nb * np);
    bk.kr.resize(nb * np);

    bk.lphaseNum.resize(nb);
    bk.lNt.resize(nb);
    bk.lNi.resize(nb * nc);
    bk.lvf.resize(nb);
    bk.lT.resize(nb);
    bk.lP.resize(nb);
    bk.lPj.resize(nb * np);
    bk.lPc.resize(nb * np);
    bk.lphaseExist.resize(nb * np);
    bk.lS.resize(nb * np);
    bk.lxij.resize(nb * np * nc);
    bk.lrho.resize(nb * np);
    bk.lxi.resize(nb * np);
    bk.lmu.resize(nb * np);
    bk.lkr.resize(nb * np);

    // derivatives
    bk.vfP.resize(nb);
    bk.vfi.resize(nb * nc);
    bk.rhoP.resize(nb * np);
    bk.rhox.resize(nb * nc * np);
    bk.xiP.resize(nb * np);
    bk.xix.resize(nb * nc * np);
    bk.muP.resize(nb * np);
    bk.mux.resize(nb * nc * np);
    bk.dPcdS.resize(nb * np * np);
    bk.dKrdS.resize(nb * np * np);

    bk.lvfP.resize(nb);
    bk.lvfi.resize(nb * nc);
    bk.lrhoP.resize(nb * np);
    bk.lrhox.resize(nb * nc * np);
    bk.lxiP.resize(nb * np);
    bk.lxix.resize(nb * nc * np);
    bk.lmuP.resize(nb * np);
    bk.lmux.resize(nb * nc * np);
    bk.ldPcdS.resize(nb * np * np);
    bk.ldKrdS.resize(nb * np * np);

    // FIM-Specified
    bk.maxLendSdP = (nc + 1) * (nc + 1) * np;
    bk.dSec_dPri.resize(nb * bk.maxLendSdP);

    bk.ldSec_dPri.resize(nb * bk.maxLendSdP);

    // Allocate Residual
    bk.res.Setup_IsoT(bk.numBulkInterior, rs.allWells.numWell, nc);

    // NR
    bk.NRstep.resize(nb);
    bk.NRphaseNum.resize(nb);
    bk.dSNR.resize(nb * np);
    bk.dSNRP.resize(nb * np);
    bk.dNNR.resize(nb * nc);
    bk.dPNR.resize(nb);

    // BulkConn
    BulkConn& conn = rs.conn;

    conn.bcval.upblock.resize(conn.numConn* np);
    conn.bcval.rho.resize(conn.numConn* np);
    conn.bcval.velocity.resize(conn.numConn* np);
}

void IsoT_FIM::AllocateLinearSystem(LinearSystem&     ls,
                                    const Reservoir&  rs,
                                    const OCPControl& ctrl)
{
    ls.SetupDomain(rs.domain);
    ls.AllocateRowMem(rs.GetComNum() + 1);
    ls.AllocateColMem();
    ls.SetupLinearSolver(ISOTHERMALMODEL, ctrl.GetWorkDir(), ctrl.GetLsFile());
}

void IsoT_FIM::InitFlash(Bulk& bk) const
{
    for (OCP_USI n = 0; n < bk.numBulk; n++) {
        bk.flashCal[bk.PVTNUM[n]]->InitFlashFIM(bk.P[n], bk.Pb[n], bk.T[n],
                                                &bk.S[n * bk.numPhase], bk.rockVp[n],
                                                bk.Ni.data() + n * bk.numCom, n);
        for (USI i = 0; i < bk.numCom; i++) {
            bk.Ni[n * bk.numCom + i] = bk.flashCal[bk.PVTNUM[n]]->GetNi(i);
        }
        PassFlashValue(bk, n);
    }
}

void IsoT_FIM::CalFlash(Bulk& bk)
{
    bk.maxNRdSSP       = 0;
    bk.index_maxNRdSSP = 0;

    for (OCP_USI n = 0; n < bk.numBulk; n++) {

        bk.flashCal[bk.PVTNUM[n]]->FlashFIM(bk.P[n], bk.T[n], &bk.Ni[n * bk.numCom],
                                            &bk.S[n * bk.numPhase], bk.phaseNum[n],
                                            &bk.xij[n * bk.numPhase * bk.numCom], n);
        PassFlashValue(bk, n);
    }
}

void IsoT_FIM::PassFlashValue(Bulk& bk, const OCP_USI& n) const
{
    const USI     np     = bk.numPhase;
    const USI     nc     = bk.numCom;
    const OCP_USI bIdp   = n * np;
    const USI     pvtnum = bk.PVTNUM[n];

    bk.phaseNum[n] = 0;
    bk.Nt[n]       = bk.flashCal[pvtnum]->GetNt();
    bk.vf[n]       = bk.flashCal[pvtnum]->GetVf();

    for (USI j = 0; j < np; j++) {
        // Important! Saturation must be passed no matter if the phase exists. This is
        // because it will be used to calculate relative permeability and capillary
        // pressure at each time step. Make sure that all saturations are updated at
        // each step!
        bk.S[bIdp + j]    = bk.flashCal[pvtnum]->GetS(j);
        bk.dSNR[bIdp + j] = bk.S[bIdp + j] - bk.dSNR[bIdp + j];
        if (bk.phaseExist[bIdp + j]) {
            if (fabs(bk.maxNRdSSP) < fabs(bk.dSNR[bIdp + j] - bk.dSNRP[bIdp + j])) {
                bk.maxNRdSSP       = bk.dSNR[bIdp + j] - bk.dSNRP[bIdp + j];
                bk.index_maxNRdSSP = n;
            }
        }

        bk.phaseExist[bIdp + j] = bk.flashCal[pvtnum]->GetPhaseExist(j);
        if (bk.phaseExist[bIdp + j]) {
            bk.phaseNum[n]++;
            bk.rho[bIdp + j] = bk.flashCal[pvtnum]->GetRho(j);
            bk.xi[bIdp + j]  = bk.flashCal[pvtnum]->GetXi(j);
            bk.mu[bIdp + j]  = bk.flashCal[pvtnum]->GetMu(j);

            // Derivatives
            bk.rhoP[bIdp + j] = bk.flashCal[pvtnum]->GetRhoP(j);
            bk.xiP[bIdp + j]  = bk.flashCal[pvtnum]->GetXiP(j);
            bk.muP[bIdp + j]  = bk.flashCal[pvtnum]->GetMuP(j);

            for (USI i = 0; i < nc; i++) {
                bk.xij[bIdp * nc + j * nc + i]  = bk.flashCal[pvtnum]->GetXij(j, i);
                bk.rhox[bIdp * nc + j * nc + i] = bk.flashCal[pvtnum]->GetRhoX(j, i);
                bk.xix[bIdp * nc + j * nc + i]  = bk.flashCal[pvtnum]->GetXiX(j, i);
                bk.mux[bIdp * nc + j * nc + i]  = bk.flashCal[pvtnum]->GetMuX(j, i);
            }
        }
    }

    bk.vfP[n] = bk.flashCal[pvtnum]->GetVfP();
    for (USI i = 0; i < nc; i++) {
        bk.vfi[n * nc + i] = bk.flashCal[pvtnum]->GetVfi(i);
    }

    copy(bk.flashCal[pvtnum]->GetDXsDXp().begin(), bk.flashCal[pvtnum]->GetDXsDXp().end(), &bk.dSec_dPri[n * bk.maxLendSdP]);
}

void IsoT_FIM::CalKrPc(Bulk& bk) const
{
    const USI& np = bk.numPhase;
    for (OCP_USI n = 0; n < bk.numBulk; n++) {
        const OCP_USI bId = n * np;
        bk.flow[bk.SATNUM[n]]->CalKrPcFIM(&bk.S[bId], n);
        copy(bk.flow[bk.SATNUM[n]]->GetKr().begin(), bk.flow[bk.SATNUM[n]]->GetKr().end(), &bk.kr[bId]);
        copy(bk.flow[bk.SATNUM[n]]->GetPc().begin(), bk.flow[bk.SATNUM[n]]->GetPc().end(), &bk.Pc[bId]);
        copy(bk.flow[bk.SATNUM[n]]->GetdKrdS().begin(), bk.flow[bk.SATNUM[n]]->GetdKrdS().end(), &bk.dKrdS[bId * np]);
        copy(bk.flow[bk.SATNUM[n]]->GetdPcdS().begin(), bk.flow[bk.SATNUM[n]]->GetdPcdS().end(), &bk.dPcdS[bId * np]);
        for (USI j = 0; j < np; j++) bk.Pj[bId + j] = bk.P[n] + bk.Pc[bId + j];
    }
}

void IsoT_FIM::CalRes(Reservoir& rs, const OCP_DBL& dt, const OCP_BOOL& resetRes0) const
{
    const Bulk&            bk    = rs.bulk;
    const USI              nb    = bk.numBulkInterior;
    const USI              np    = bk.numPhase;
    const USI              nc    = bk.numCom;
    const USI              len   = nc + 1;
    OCPRes&                Res   = bk.res;
    
    Res.SetZero();

    // Bulk to Bulk
    OCP_USI bId, eId, bIdb;    
    // Accumalation Term
    for (OCP_USI n = 0; n < nb; n++) {
        bId             = n * len;
        bIdb            = n * nc;
        Res.resAbs[bId] = bk.rockVp[n] - bk.vf[n];
        for (USI i = 0; i < nc; i++) {
            Res.resAbs[bId + 1 + i] = bk.Ni[bIdb + i] - bk.lNi[bIdb + i];
        }
    }

    // Flux Term
    BulkConn&              conn  = rs.conn;
    BulkConnVal&           bcval = conn.bcval;
    vector<OCPFlux*>&      flux  = conn.flux;
    USI     cType;
    for (OCP_USI c = 0; c < conn.numConn; c++) {

        bId   = conn.iteratorConn[c].BId();
        eId   = conn.iteratorConn[c].EId();
        cType = conn.iteratorConn[c].Type();

        flux[cType]->CalFlux(conn.iteratorConn[c], bk);
        copy(flux[cType]->GetUpblock().begin(), flux[cType]->GetUpblock().end(), &bcval.upblock[c * np]);
        copy(flux[cType]->GetRho().begin(), flux[cType]->GetRho().end(), &bcval.rho[c * np]);
        copy(flux[cType]->GetFluxVj().begin(), flux[cType]->GetFluxVj().end(), &bcval.velocity[c * np]);
               
        if (eId < nb) {
            for (USI i = 0; i < nc; i++) {               
                Res.resAbs[bId * len + 1 + i] += dt * flux[cType]->GetFluxNi()[i];
                Res.resAbs[eId * len + 1 + i] -= dt * flux[cType]->GetFluxNi()[i];
            }
        }
        else {
            for (USI i = 0; i < nc; i++) {
                Res.resAbs[bId * len + 1 + i] += dt * flux[cType]->GetFluxNi()[i];
            }
        }
    }

    // Well to Bulk, Well
    USI wId = nb * len;
    for (const auto& wl : rs.allWells.wells) {
        wl.CalResFIM(wId, Res, bk, dt);
    }

    // Calculate RelRes
    OCP_DBL tmp;
    for (OCP_USI n = 0; n < nb; n++) {

        for (USI i = 0; i < len; i++) {
            tmp = fabs(Res.resAbs[n * len + i] / bk.rockVp[n]);
            if (Res.maxRelRes_V < tmp) {
                Res.maxRelRes_V = tmp;
                Res.maxId_V     = n;
            }
            Res.resRelV[n] += tmp * tmp;
        }
        Res.resRelV[n] = sqrt(Res.resRelV[n]);

        for (USI i = 1; i < len; i++) {
            tmp = fabs(Res.resAbs[n * len + i] / bk.Nt[n]);
            if (Res.maxRelRes_N < tmp) {
                Res.maxRelRes_N = tmp;
                Res.maxId_N     = n;
            }
            Res.resRelN[n] += tmp * tmp;
        }
        Res.resRelN[n] = sqrt(Res.resRelN[n]);
    }

    Dscalar(Res.resAbs.size(), -1.0, Res.resAbs.data());
    if (resetRes0) {
        Res.SetInitRes();

        GetWallTime timer;
        timer.Start();

        OCP_DBL tmploc = Res.maxRelRes0_V;
        MPI_Allreduce(&tmploc, &Res.maxRelRes0_V, 1, MPI_DOUBLE, MPI_MIN, rs.domain.myComm);

        OCPTIME_COMM_COLLECTIVE += timer.Stop() / 1000;
    }
}

void IsoT_FIM::AssembleMatBulks(LinearSystem&    ls,
                                const Reservoir& rs,
                                const OCP_DBL&   dt) const
{
    const USI numWell = rs.GetNumOpenWell();

    const Bulk&     bk     = rs.bulk;
    const BulkConn& conn   = rs.conn;
    const OCP_USI   nb     = bk.numBulkInterior;
    const USI       np     = bk.numPhase;
    const USI       nc     = bk.numCom;
    const USI       ncol   = nc + 1;
    const USI       ncol2  = np * nc + np;
    const USI       bsize  = ncol * ncol;
    const USI       bsize2 = ncol * ncol2;

    ls.AddDim(nb);

    vector<OCP_DBL> bmat(bsize, 0);
    // Accumulation term
    for (USI i = 1; i < ncol; i++) {
        bmat[i * ncol + i] = 1;
    }
    for (OCP_USI n = 0; n < nb; n++) {
        bmat[0] = bk.v[n] * bk.poroP[n] - bk.vfP[n];
        for (USI i = 0; i < nc; i++) {
            bmat[i + 1] = -bk.vfi[n * nc + i];
        }
        ls.NewDiag(n, bmat);
    }

    // flux term
    vector<OCPFlux*>& flux = conn.flux;
    OCP_USI  bId, eId;
    USI      cType;
    for (OCP_USI c = 0; c < conn.numConn; c++) {

        bId   = conn.iteratorConn[c].BId();
        eId   = conn.iteratorConn[c].EId();
        cType = conn.iteratorConn[c].Type();

        flux[cType]->AssembleMatFIM(conn.iteratorConn[c], c, conn.bcval, bk);
        
        bmat = flux[cType]->GetdFdXpB();
        DaABpbC(ncol, ncol, ncol2, 1, flux[cType]->GetdFdXsB().data(), &bk.dSec_dPri[bId * bsize2], 1,
            bmat.data());       
        Dscalar(bsize, dt, bmat.data());
        // Assemble
        // Begin - Begin -- add
        ls.AddDiag(bId, bmat);
        // End - Begin -- insert
        if (eId < nb) {
            // Interior grid
            Dscalar(bsize, -1, bmat.data());
            ls.NewOffDiag(eId, bId, bmat);
        }

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or INF in bmat !");
        }
#endif

        // End
        bmat = flux[cType]->GetdFdXpE();
        DaABpbC(ncol, ncol, ncol2, 1, flux[cType]->GetdFdXsE().data(), &bk.dSec_dPri[eId * bsize2], 1,
            bmat.data());
        Dscalar(bsize, dt, bmat.data());
        
        if (eId < nb) {
            // Interior grid
            // Begin - End -- insert
            ls.NewOffDiag(bId, eId, bmat);
            // End - End -- add
            Dscalar(bsize, -1, bmat.data());
            ls.AddDiag(eId, bmat);
        }
        else {
            // ghost grid
            // Begin - End -- insert
            ls.NewOffDiag(bId, eId + numWell, bmat);
        }

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or INF in bmat !");
        }
#endif
    }
}


void IsoT_FIM::AssembleMatWells(LinearSystem&    ls,
                                const Reservoir& rs,
                                const OCP_DBL&   dt) const
{
    for (auto& wl : rs.allWells.wells) {
        wl.AssembleMatFIM(ls, rs.bulk, dt);
    }

    //// for Reinjection
    // for (auto& wG : rs.allWells.wellGroup) {
    //     if (wG.IfReInj()) {
    //         for (auto& prod : rs.allWells.wellGroup[wG.prodGroup].wIdPROD) {
    //             if (rs.allWells.wells[prod].IsOpen()) {
    //                 rs.allWells.wells[prod].AssembleMatReinjection_FIM(rs.bulk, ls,
    //                 dt, rs.allWells.wells,
    //                     wG.wIdINJ);
    //             }
    //         }
    //     }
    // }
}


void IsoT_FIM::GetSolution(Reservoir&             rs,
                           vector<OCP_DBL>& u,
                           const OCPControl&      ctrl) const
{
    const Domain&   domain = rs.domain;
    Bulk&           bk     = rs.bulk;
    const OCP_USI   nb     = bk.numBulk;
    const USI       np     = bk.numPhase;
    const USI       nc     = bk.numCom;
    const USI       row    = np * (nc + 1);
    const USI       col    = nc + 1;

    // Well first
    USI wId = bk.numBulkInterior * col;
    for (auto& wl : rs.allWells.wells) {
        if (wl.IsOpen()) {
            wl.SetBHP(wl.BHP() + u[wId]);
            wId += col;
        }
    }

    GetWallTime timerT;         ///< total timer
    GetWallTime timerC;         ///< calculation timer
    OCP_DBL     time_cal = 0;   ///< calculation time
    timerT.Start();
    

    // Exchange Solution for ghost grid
    for (USI i = 0; i < domain.numRecvProc; i++) {
        const vector<OCP_USI>& rel = domain.recv_element_loc[i];
        MPI_Irecv(&u[rel[1] * col], (rel[2] - rel[1]) * col, MPI_DOUBLE, rel[0], 0, domain.myComm, &domain.recv_request[i]);
    }
   
    vector<vector<OCP_DBL>> send_buffer(domain.numSendProc);
    for (USI i = 0; i < domain.numSendProc; i++) {
        const vector<OCP_USI>& sel = domain.send_element_loc[i];
        vector<OCP_DBL>&       s   = send_buffer[i];
        s.resize(1 + (sel.size() - 1) * col);
        s[0] = sel[0];
        for (USI j = 1; j < sel.size(); j++) {
             const OCP_DBL* bId = u.data() + sel[j] * col;
             copy(bId, bId + col, &s[1 + (j - 1) * col]);
        }
        MPI_Isend(s.data() + 1, s.size() - 1, MPI_DOUBLE, s[0], 0, domain.myComm, &domain.send_request[i]);
    }

    // Bulk
    const OCP_DBL dSmaxlim = ctrl.ctrlNR.NRdSmax;
    // const OCP_DBL dPmaxlim = ctrl.ctrlNR.NRdPmax;

    vector<OCP_DBL> dtmp(row, 0);
    OCP_DBL         chopmin = 1;
    OCP_DBL         choptmp = 0;

    bk.dSNR       = bk.S;
    bk.NRphaseNum = bk.phaseNum;
    bk.NRdPmax    = 0;
    bk.NRdNmax    = 0;

    OCP_USI bId = 0;
    OCP_USI eId = bk.GetInteriorBulkNum();

    // interior first, ghost second
    for (USI p = 0; p < 2; p++) {
  
        timerC.Start();

		for (OCP_USI n = bId; n < eId; n++) {
			// const vector<OCP_DBL>& scm = satcm[SATNUM[n]];

			chopmin = 1;
			// compute the chop
			fill(dtmp.begin(), dtmp.end(), 0.0);

			DaAxpby(row, col, 1, &bk.dSec_dPri[n * bk.maxLendSdP], u.data() + n * col, 1,
				dtmp.data());

			for (USI j = 0; j < np; j++) {
				choptmp = 1;
				if (fabs(dtmp[j]) > dSmaxlim) {
					choptmp = dSmaxlim / fabs(dtmp[j]);
				}
				else if (bk.S[n * np + j] + dtmp[j] < 0.0) {
					choptmp = 0.9 * bk.S[n * np + j] / fabs(dtmp[j]);
				}
				// if (fabs(S[n_np_j] - scm[j]) > TINY &&
				//     (S[n_np_j] - scm[j]) / (choptmp * dtmp[js]) < 0)
				//     choptmp *= min(1.0, -((S[n_np_j] - scm[j]) / (choptmp * dtmp[js])));
				chopmin = min(chopmin, choptmp);
			}

			// dS
			for (USI j = 0; j < np; j++) {
				bk.dSNRP[n * np + j] = chopmin * dtmp[j];
				bk.S[n * np + j] += bk.dSNRP[n * np + j];
			}

			// dxij   ---- Compositional model only

			if (bk.IfUseEoS()) {
				USI js = np;
				if (bk.phaseNum[n] >= 3) {
					// num of Hydroncarbon phase >= 2
					OCP_USI bId = 0;
					for (USI j = 0; j < 2; j++) {
						bId = n * np * nc + j * nc;
						for (USI i = 0; i < bk.numComH; i++) {
							bk.xij[bId + i] += chopmin * dtmp[js];
							js++;
						}
						js++;
					}
				}
			}

			// dP
			OCP_DBL dP = u[n * col];
			// choptmp = dPmaxlim / fabs(dP);
			// chopmin = min(chopmin, choptmp);
			if (fabs(bk.NRdPmax) < fabs(dP)) bk.NRdPmax = dP;
			bk.P[n] += dP; // seems better
			bk.dPNR[n] = dP;

			// dNi
			bk.NRstep[n] = chopmin;
			for (USI i = 0; i < nc; i++) {
				bk.dNNR[n * nc + i] = u[n * col + 1 + i] * chopmin;
				if (fabs(bk.NRdNmax) < fabs(bk.dNNR[n * nc + i]) / bk.Nt[n])
					bk.NRdNmax = bk.dNNR[n * nc + i] / bk.Nt[n];

				bk.Ni[n * nc + i] += bk.dNNR[n * nc + i];

				// if (bk.Ni[n * nc + i] < 0 && bk.Ni[n * nc + i] > -1E-3) {
				//     bk.Ni[n * nc + i] = 1E-20;
				// }
			}
		}

        time_cal += timerC.Stop();

        if (p == 0) {
            bId = eId;
            eId = nb;
            MPI_Waitall(domain.numRecvProc, domain.recv_request.data(), MPI_STATUS_IGNORE);
        }
        else {
            break;
        }        
    }

    MPI_Waitall(domain.numSendProc, domain.send_request.data(), MPI_STATUS_IGNORE);

    OCPTIME_COMM_P2P += (timerT.Stop() - time_cal) / 1000;
    OCPTIME_NRSTEPC  += time_cal / 1000;
}

void IsoT_FIM::ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl)
{
    Bulk& bk = rs.bulk;

    // Rock
    bk.poro   = bk.lporo;
    bk.poroP  = bk.lporoP;
    bk.rockVp = bk.lrockVp;
    // Fluid
    bk.phaseNum   = bk.lphaseNum;
    bk.Nt         = bk.lNt;
    bk.Ni         = bk.lNi;
    bk.vf         = bk.lvf;
    bk.P          = bk.lP;
    bk.Pj         = bk.lPj;
    bk.Pc         = bk.lPc;
    bk.phaseExist = bk.lphaseExist;
    bk.S          = bk.lS;
    bk.xij        = bk.lxij;
    bk.rho        = bk.lrho;
    bk.xi         = bk.lxi;
    bk.mu         = bk.lmu;
    bk.kr         = bk.lkr;
    // derivatives
    bk.vfP     = bk.lvfP;
    bk.vfi     = bk.lvfi;
    bk.rhoP    = bk.lrhoP;
    bk.rhox    = bk.lrhox;
    bk.xiP     = bk.lxiP;
    bk.xix     = bk.lxix;
    bk.muP     = bk.lmuP;
    bk.mux     = bk.lmux;
    bk.dPcdS   = bk.ldPcdS;
    bk.dKrdS   = bk.ldKrdS;
    // FIM-Specified
    bk.dSec_dPri    = bk.ldSec_dPri;

    // Wells
    rs.allWells.ResetBHP();
    rs.allWells.CalTrans(bk);
    rs.allWells.CaldG(bk);
    rs.allWells.CalFlux(bk);

    // Optional Features
    rs.optFeatures.ResetToLastTimeStep();

    // Iters
    ctrl.ResetIterNRLS();

    // Residual
    CalRes(rs, ctrl.GetCurDt(), OCP_TRUE);
}

void IsoT_FIM::UpdateLastTimeStep(Reservoir& rs) const
{
    Bulk& bk = rs.bulk;

    // Rock
    bk.lporo   = bk.poro;
    bk.lporoP  = bk.poroP;
    bk.lrockVp = bk.rockVp;

    // Fluid
    bk.lphaseNum   = bk.phaseNum;
    bk.lNt         = bk.Nt;
    bk.lNi         = bk.Ni;
    bk.lvf         = bk.vf;
    bk.lP          = bk.P;
    bk.lPj         = bk.Pj;
    bk.lPc         = bk.Pc;
    bk.lphaseExist = bk.phaseExist;
    bk.lS          = bk.S;
    bk.lxij        = bk.xij;
    bk.lrho        = bk.rho;
    bk.lxi         = bk.xi;
    bk.lmu         = bk.mu;
    bk.lkr         = bk.kr;

    // derivatives
    bk.lvfP     = bk.vfP;
    bk.lvfi     = bk.vfi;
    bk.lrhoP    = bk.rhoP;
    bk.lrhox    = bk.rhox;
    bk.lxiP     = bk.xiP;
    bk.lxix     = bk.xix;
    bk.lmuP     = bk.muP;
    bk.lmux     = bk.mux;
    bk.ldPcdS = bk.dPcdS;
    bk.ldKrdS  = bk.dKrdS;

    // FIM-Specified
    bk.ldSec_dPri    = bk.dSec_dPri;

    rs.allWells.UpdateLastTimeStepBHP();
    rs.optFeatures.UpdateLastTimeStep();
}

////////////////////////////////////////////
// IsoT_AIMc
////////////////////////////////////////////

void IsoT_AIMc::Setup(Reservoir& rs, LinearSystem& ls, const OCPControl& ctrl)
{
    // Allocate Bulk and BulkConn Memory
    AllocateReservoir(rs);
    // Setup neighbor
    SetupNeighbor(rs);
    // Allocate memory for internal matrix structure
    IsoT_FIM::AllocateLinearSystem(ls, rs, ctrl);
}


void IsoT_AIMc::SetupNeighbor(Reservoir& rs)
{
    // Note: 
    // for interior bulk: neighbor stores their all 1-neighbors
    // for ghost    bulk: neighbor stores their 1-neighbors belonging to current process

    const OCP_USI nb   = rs.GetBulkNum();
    OCP_USI       bId, eId;
    BulkConn&     conn = rs.conn;
    
    conn.neighbor.resize(nb);
    for (const auto& c : conn.iteratorConn) {
        bId = c.BId();
        eId = c.EId();

        conn.neighbor[bId].push_back(eId);
        conn.neighbor[eId].push_back(bId);
    }
}


void IsoT_AIMc::InitReservoir(Reservoir& rs) const
{
    rs.bulk.InitPTSw(50);

    CalRock(rs.bulk);

    IsoT_IMPEC::InitFlash(rs.bulk);
    IsoT_IMPEC::CalKrPc(rs.bulk);

    rs.allWells.InitBHP(rs.bulk);

    UpdateLastTimeStep(rs);
}

void IsoT_AIMc::Prepare(Reservoir& rs, const OCP_DBL& dt)
{
    rs.allWells.PrepareWell(rs.bulk);
    CalRes(rs, dt, OCP_TRUE);

    // Set FIM Bulk
    rs.CalCFL(dt, OCP_FALSE);
    rs.allWells.SetupWellBulk(rs.bulk);
    SetFIMBulk(rs);
    //  Calculate FIM Bulk properties
    CalFlashI(rs.bulk);
    CalKrPcI(rs.bulk);

    UpdateLastTimeStep(rs);
}

void IsoT_AIMc::AssembleMat(LinearSystem&    ls,
                            const Reservoir& rs,
                            const OCP_DBL&   dt) const
{
    AssembleMatBulks(ls, rs, dt);
    IsoT_FIM::AssembleMatWells(ls, rs, dt);
    ls.AssembleRhsCopy(rs.bulk.res.resAbs);
}

void IsoT_AIMc::SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl)
{
#ifdef DEBUG
    ls.CheckEquation();
#endif // DEBUG

    GetWallTime timer;

    timer.Start();
    ls.CalCommTerm(rs.GetNumOpenWell());
    ls.AssembleMatLinearSolver();
    OCPTIME_ASSEMBLE_MAT_FOR_LS += timer.Stop() / 1000;

    timer.Start();
    int status = ls.Solve();
    if (status < 0) {
        status = ls.GetNumIters();
    }

#ifdef DEBUG
    //OCP_INT myrank = rs.domain.myrank;
    //ls.OutputLinearSystem("proc" + to_string(CURRENT_RANK) + "_testA_AIMc.out",
    //                      "proc" + to_string(CURRENT_RANK) + "_testb_AIMc.out");
    //MPI_Barrier(rs.domain.myComm);
    //ls.OutputSolution("proc" + to_string(CURRENT_RANK) + "_testx_AIMc.out");
    // ls.CheckSolution();
#endif // DEBUG

    OCPTIME_LSOLVER += timer.Stop() / 1000;
    ctrl.UpdateIterLS(status);
    ctrl.UpdateIterNR();

    timer.Start();
    GetSolution(rs, ls.GetSolution(), ctrl);
    OCPTIME_NRSTEP += timer.Stop() / 1000;
    ls.ClearData();
}

OCP_BOOL IsoT_AIMc::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    const OCP_DBL& dt = ctrl.current_dt;

    // First check: Ni check and bulk Pressure check
    if (!ctrl.Check(rs, {"BulkNi", "BulkP"})) {
        ResetToLastTimeStep(rs, ctrl);
        cout << "Cut time step size and repeat! current dt = " << fixed
             << setprecision(3) << dt << " days\n";
        return OCP_FALSE;
    }

    CalFlashI(rs.bulk);
    CalFlashEp(rs.bulk);
    CalKrPcI(rs.bulk);

    CalRock(rs.bulk);

    rs.allWells.CalTrans(rs.bulk);
    rs.allWells.CalFlux(rs.bulk);

    CalRes(rs, dt, OCP_FALSE);
    return OCP_TRUE;
}

OCP_BOOL IsoT_AIMc::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    OCP_USI       dSn;
    const OCP_DBL NRdSmax = rs.GetNRdSmax(dSn);
    const OCP_DBL NRdPmax = rs.GetNRdPmax();
    // const OCP_DBL NRdNmax = rs.GetNRdNmax();

#ifdef DEBUG
    // cout << "### DEBUG: Residuals = " << setprecision(3) << scientific
    //      << resAIMc.maxRelRes0_V << "  " << resAIMc.maxRelRes_V << "  "
    //      << resAIMc.maxRelRes_N << "  " << NRdPmax << "  " << NRdSmax << endl;
#endif

    OCP_INT conflag_loc = -1;
    if (((rs.bulk.res.maxRelRes_V <= rs.bulk.res.maxRelRes0_V * ctrl.ctrlNR.NRtol ||
        rs.bulk.res.maxRelRes_V <= ctrl.ctrlNR.NRtol ||
        rs.bulk.res.maxRelRes_N <= ctrl.ctrlNR.NRtol) &&
        rs.bulk.res.maxWellRelRes_mol <= ctrl.ctrlNR.NRtol) ||
        (fabs(NRdPmax) <= ctrl.ctrlNR.NRdPmin &&
            fabs(NRdSmax) <= ctrl.ctrlNR.NRdSmin)) {
        conflag_loc = 0;
    }

    GetWallTime timer;
    timer.Start();

    OCP_INT conflag;
    MPI_Allreduce(&conflag_loc, &conflag, 1, MPI_INT, MPI_MIN, rs.domain.myComm);

    OCPTIME_COMM_COLLECTIVE += timer.Stop() / 1000;

    if (conflag == 0) {

        if (!ctrl.Check(rs, {"WellP"})) {
            ResetToLastTimeStep(rs, ctrl);
            return OCP_FALSE;
        } else {
            CalFlashEa(rs.bulk);
            CalKrPcE(rs.bulk);
            return OCP_TRUE;
        }

    } else if (ctrl.iterNR > ctrl.ctrlNR.maxNRiter) {
        ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
        ResetToLastTimeStep(rs, ctrl);
        cout << "### WARNING: NR not fully converged! Cut time step size and repeat!  "
                "current dt = "
             << fixed << setprecision(3) << ctrl.current_dt << " days\n";
        return OCP_FALSE;
    } else {
        return OCP_FALSE;
    }
}

/// Finish a time step.
void IsoT_AIMc::FinishStep(Reservoir& rs, OCPControl& ctrl) const
{
    rs.CalIPRT(ctrl.GetCurDt());
    rs.CalMaxChange();
    UpdateLastTimeStep(rs);
    ctrl.CalNextTimeStep(rs, {"dP", "dS", "iter"});
}

/// Allocate memory for reservoir
void IsoT_AIMc::AllocateReservoir(Reservoir& rs)
{
    IsoT_FIM::AllocateReservoir(rs);

    Bulk&         bk = rs.bulk;
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;
    const USI     nc = bk.numCom;

    bk.vj.resize(nb * np);
    bk.lvj.resize(nb * np);

    bk.xijNR.resize(nb * np * nc);
    bk.cfl.resize(nb * np);
    bk.bulkTypeAIM.Setup(nb);
}

void IsoT_AIMc::SetFIMBulk(Reservoir& rs)
{
    // IMPORTANT: implicity of the same grid in different processes should be consistent
    
    const OCP_INT nlayers = 2;

    // We just consider at most 1 layer neighbor now

    Bulk&           bk   = rs.bulk;
    const BulkConn& conn = rs.conn;
    const OCP_USI   nb   = bk.numBulkInterior;
    const USI       np   = bk.numPhase;
    const USI       nc   = bk.numCom;

    // all impec
    bk.bulkTypeAIM.Init(-1);

    OCP_USI  bIdp, bIdc;
    OCP_BOOL flag;

    for (OCP_USI n = 0; n < nb; n++) {
        bIdp = n * np;
        bIdc = n * nc;
        flag = OCP_FALSE;
        // CFL
        for (USI j = 0; j < np; j++) {
            if (bk.cfl[bIdp + j] > 0.8) {
                flag = OCP_TRUE;
                break;
            }
        }
        // Volume error
        if (!flag) {
            if ((fabs(bk.vf[n] - bk.rockVp[n]) / bk.rockVp[n]) > 1E-3) {
                flag = OCP_TRUE;
            }
        }

        // NR Step
        if (!flag && OCP_FALSE) {
            // dP
            if (fabs(bk.dPNR[n] / bk.P[n]) > 1E-3) {
                flag = OCP_TRUE;
            }
            // dNi
            if (!flag) {
                for (USI i = 0; i < bk.numCom; i++) {
                    if (fabs(bk.dNNR[bIdc + i] / bk.Ni[bIdc + i]) > 1E-3) {
                        flag = OCP_TRUE;
                        break;
                    }
                }
            }
        }

        if (flag) {
            SetKNeighbor(conn.neighbor, n, bk.bulkTypeAIM, nlayers);
        }
    }

    // add WellBulk's 1-neighbor as Implicit bulk
    for (auto& p : bk.wellBulkId) {
        SetKNeighbor(conn.neighbor, p, bk.bulkTypeAIM, nlayers);
    }

    // exchange information of implicity of grid
    const Domain& domain = rs.domain;

    vector<vector<OCP_INT>> recv_buffer(domain.numRecvProc);
    for (USI i = 0; i < domain.numRecvProc; i++) {
        const vector<OCP_USI>& rel = domain.recv_element_loc[i];
        vector<OCP_INT>&       r   = recv_buffer[i];
        r.resize(rel[2] - rel[1]);
        MPI_Irecv(&r[0], rel[2] - rel[1], MPI_INT, rel[0], 0, domain.myComm, &domain.recv_request[i]);
    }

    vector<vector<OCP_INT>> send_buffer(domain.numSendProc);
    for (USI i = 0; i < domain.numSendProc; i++) {
        const vector<OCP_USI>& sel = domain.send_element_loc[i];
        vector<OCP_INT>&       s   = send_buffer[i];
        s.resize(sel.size());
        s[0] = sel[0];
        for (USI j = 1; j < sel.size(); j++) {
            s[j] = bk.bulkTypeAIM.GetBulkType(sel[j]);
        }
        MPI_Isend(s.data() + 1, s.size() - 1, MPI_INT, s[0], 0, domain.myComm, &domain.send_request[i]);
    }

    MPI_Waitall(domain.numRecvProc, domain.recv_request.data(), MPI_STATUS_IGNORE);
    for (USI i = 0; i < domain.numRecvProc; i++) {
        const vector<OCP_USI>& rel = domain.recv_element_loc[i];
        const vector<OCP_INT>& r   = recv_buffer[i];

        for (OCP_USI n = 0; n < rel[2] - rel[1]; n++) {
            SetKNeighbor(conn.neighbor, n + rel[1], bk.bulkTypeAIM, r[n]);
        }
    }

    MPI_Waitall(domain.numSendProc, domain.send_request.data(), MPI_STATUS_IGNORE);

    // Check Consistency
    for (USI i = 0; i < domain.numRecvProc; i++) {
        const vector<OCP_USI>& rel = domain.recv_element_loc[i];
        vector<OCP_INT>&       r   = recv_buffer[i];
        r.resize(rel[2] - rel[1]);
        MPI_Irecv(&r[0], rel[2] - rel[1], MPI_INT, rel[0], 0, domain.myComm, &domain.recv_request[i]);
    }

    for (USI i = 0; i < domain.numSendProc; i++) {
        const vector<OCP_USI>& sel = domain.send_element_loc[i];
        vector<OCP_INT>&       s   = send_buffer[i];
        s.resize(sel.size());
        s[0] = sel[0];
        for (USI j = 1; j < sel.size(); j++) {
            s[j] = bk.bulkTypeAIM.GetBulkType(sel[j]);
        }
        MPI_Isend(s.data() + 1, s.size() - 1, MPI_INT, s[0], 0, domain.myComm, &domain.send_request[i]);
    }

    MPI_Waitall(domain.numRecvProc, domain.recv_request.data(), MPI_STATUS_IGNORE);
     
    for (USI i = 0; i < domain.numRecvProc; i++) {
        const vector<OCP_USI>& rel = domain.recv_element_loc[i];
        vector<OCP_INT>&       r   = recv_buffer[i];
        for (OCP_USI n = 0; n < rel[2] - rel[1]; n++) {
            bk.bulkTypeAIM.SetBulkType(n + rel[1], r[n]);     // Maybe not a good idea
        }
    }

    MPI_Waitall(domain.numSendProc, domain.send_request.data(), MPI_STATUS_IGNORE);

    if (OCP_TRUE) {
        cout << fixed << setprecision(2) << "Rank " << CURRENT_RANK << "  " << bk.bulkTypeAIM.GetNumFIMBulk() * 1.0 / bk.numBulk * 100 << "% " << endl;
    }
}


void IsoT_AIMc::SetKNeighbor(const vector<vector<OCP_USI>>& neighbor, const OCP_USI& p, BulkTypeAIM& tar, OCP_INT k)
{
    tar.SetBulkType(p, max(k, tar.GetBulkType(p)));
    if (k > 0) {
        k--;
        for (const auto& v : neighbor[p]) {
            SetKNeighbor(neighbor, v, tar, k);
        }
    }
}


void IsoT_AIMc::CalFlashEp(Bulk& bk)
{
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;
    const USI     nc = bk.numCom;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bulkTypeAIM.IfIMPECbulk(n)) {
            // Explicit bulk

            bk.flashCal[bk.PVTNUM[n]]->FlashIMPEC(bk.P[n], bk.T[n], &bk.Ni[n * nc],
                                                  bk.phaseNum[n],
                                                  &bk.xijNR[n * np * nc], n);
            // bk.PassFlashValueAIMcEp(n);
            PassFlashValueEp(bk, n);
        }
    }
}

void IsoT_AIMc::CalFlashEa(Bulk& bk)
{
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;
    const USI     nc = bk.numCom;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bulkTypeAIM.IfIMPECbulk(n)) {
            // Explicit bulk

            bk.flashCal[bk.PVTNUM[n]]->FlashIMPEC(bk.P[n], bk.T[n], &bk.Ni[n * nc],
                                                  bk.phaseNum[n], &bk.xij[n * np * nc],
                                                  n);
            // bk.PassFlashValueAIMcEa(n);

            IsoT_IMPEC::PassFlashValue(bk, n);
        }
    }
}

void IsoT_AIMc::CalFlashI(Bulk& bk)
{
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;
    const USI     nc = bk.numCom;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bulkTypeAIM.IfFIMbulk(n)) {
            // Implicit bulk

            bk.flashCal[bk.PVTNUM[n]]->FlashFIM(bk.P[n], bk.T[n], &bk.Ni[n * nc],
                                                &bk.S[n * np], bk.phaseNum[n],
                                                &bk.xij[n * np * nc], n);
            IsoT_FIM::PassFlashValue(bk, n);
            for (USI j = 0; j < np; j++) {
                bk.vj[n * np + j] = bk.vf[n] * bk.S[n * np + j];
            }
        }
    }
}

void IsoT_AIMc::PassFlashValueEp(Bulk& bk, const OCP_USI& n)
{
    // only var about volume needs, some flash var also
    OCP_FUNCNAME;

    const USI     np     = bk.numPhase;
    const USI     nc     = bk.numCom;
    const OCP_USI bIdp   = n * np;
    const USI     pvtnum = bk.PVTNUM[n];

    bk.Nt[n]  = bk.flashCal[pvtnum]->GetNt();
    bk.vf[n]  = bk.flashCal[pvtnum]->GetVf();
    bk.vfP[n] = bk.flashCal[pvtnum]->GetVfP();
    for (USI i = 0; i < nc; i++) {
        bk.vfi[n * nc + i] = bk.flashCal[pvtnum]->GetVfi(i);
    }

    bk.phaseNum[n] = 0;
    for (USI j = 0; j < np; j++) {
        if (bk.flashCal[pvtnum]->GetPhaseExist(j)) {
            bk.phaseNum[n]++;

            // IMPORTANT -- need for next Flash
            // But xij in nonlinear equations has been modified
            for (USI i = 0; i < nc; i++) {
                bk.xijNR[bIdp * nc + j * nc + i] = bk.flashCal[pvtnum]->GetXij(j, i);
            }
        }
    }
}

void IsoT_AIMc::CalKrPcE(Bulk& bk)
{
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bulkTypeAIM.IfIMPECbulk(n)) {
            // Explicit bulk
            const OCP_USI bId = n * np;
            bk.flow[bk.SATNUM[n]]->CalKrPc(&bk.S[bId], n);
            copy(bk.flow[bk.SATNUM[n]]->GetKr().begin(), bk.flow[bk.SATNUM[n]]->GetKr().end(), &bk.kr[bId]);
            copy(bk.flow[bk.SATNUM[n]]->GetPc().begin(), bk.flow[bk.SATNUM[n]]->GetPc().end(), &bk.Pc[bId]);
            for (USI j = 0; j < np; j++) bk.Pj[bId + j] = bk.P[n] + bk.Pc[bId + j];
        }
    }
}

void IsoT_AIMc::CalKrPcI(Bulk& bk)
{
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bulkTypeAIM.IfFIMbulk(n)) {
            // Implicit bulk
            const OCP_USI bId = n * np;
            bk.flow[bk.SATNUM[n]]->CalKrPcFIM(&bk.S[bId], n);
            copy(bk.flow[bk.SATNUM[n]]->GetKr().begin(), bk.flow[bk.SATNUM[n]]->GetKr().end(), &bk.kr[bId]);
            copy(bk.flow[bk.SATNUM[n]]->GetPc().begin(), bk.flow[bk.SATNUM[n]]->GetPc().end(), &bk.Pc[bId]);
            copy(bk.flow[bk.SATNUM[n]]->GetdKrdS().begin(), bk.flow[bk.SATNUM[n]]->GetdKrdS().end(), &bk.dKrdS[bId * np]);
            copy(bk.flow[bk.SATNUM[n]]->GetdPcdS().begin(), bk.flow[bk.SATNUM[n]]->GetdPcdS().end(), &bk.dPcdS[bId * np]);
            for (USI j = 0; j < np; j++) bk.Pj[bId + j] = bk.P[n] + bk.Pc[bId + j];
        }
    }
}

void IsoT_AIMc::AssembleMatBulks(LinearSystem&    ls,
                                 const Reservoir& rs,
                                 const OCP_DBL&   dt) const
{
    const USI numWell = rs.GetNumOpenWell();

    const Bulk&     bk      = rs.bulk;
    const BulkConn& conn    = rs.conn;
    const OCP_USI   nb      = bk.numBulkInterior;
    const USI       np      = bk.numPhase;
    const USI       nc      = bk.numCom;
    const USI       ncol    = nc + 1;
    const USI       ncol2   = np * nc + np;
    const USI       bsize   = ncol * ncol;
    const USI       bsize2  = ncol * ncol2;

    ls.AddDim(nb);

    vector<OCP_DBL> bmat(bsize, 0);
    // Accumulation term
    for (USI i = 1; i < nc + 1; i++) {
        bmat[i * ncol + i] = 1;
    }
    for (OCP_USI n = 0; n < nb; n++) {
        bmat[0] = bk.v[n] * bk.poroP[n] - bk.vfP[n];
        for (USI i = 0; i < nc; i++) {
            bmat[i + 1] = -bk.vfi[n * nc + i];
        }
        ls.NewDiag(n, bmat);
    }

    // flux term
    vector<OCPFlux*>& flux = conn.flux;
    OCP_BOOL          bIdFIM, eIdFIM;
    OCP_USI           bId, eId;
    USI               cType;
    
   
    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId    = conn.iteratorConn[c].BId();
        eId    = conn.iteratorConn[c].EId();
        cType  = conn.iteratorConn[c].Type();

        if (bk.bulkTypeAIM.IfFIMbulk(bId))  bIdFIM = OCP_TRUE;
        else                                bIdFIM = OCP_FALSE;

        if (bk.bulkTypeAIM.IfFIMbulk(eId))  eIdFIM = OCP_TRUE;
        else                                eIdFIM = OCP_FALSE;

        flux[cType]->AssembleMatAIM(conn.iteratorConn[c], c, conn.bcval, bk);

        // Assemble
        bmat = flux[cType]->GetdFdXpB();
        if (bIdFIM) {
            DaABpbC(ncol, ncol, ncol2, 1, flux[cType]->GetdFdXsB().data(), &bk.dSec_dPri[bId * bsize2],
                    1, bmat.data());
        }
        Dscalar(bsize, dt, bmat.data());
        // Begin - Begin -- add
        ls.AddDiag(bId, bmat);
        // End - Begin -- insert
        if (eId < nb) {
            // Interior grid
            Dscalar(bsize, -1, bmat.data());
            ls.NewOffDiag(eId, bId, bmat);
        }    

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif

        // End
        bmat = flux[cType]->GetdFdXpE();
        if (eIdFIM) {
            DaABpbC(ncol, ncol, ncol2, 1, flux[cType]->GetdFdXsE().data(), &bk.dSec_dPri[eId * bsize2],
                    1, bmat.data());
        }
        Dscalar(bsize, dt, bmat.data());

        if (eId < nb) {
            // Interior grid
            // Begin - End -- insert
            ls.NewOffDiag(bId, eId, bmat);
            // End - End -- add
            Dscalar(bsize, -1, bmat.data());
            ls.AddDiag(eId, bmat);
        }
        else {
            // Ghost grid
            ls.NewOffDiag(bId, eId + numWell, bmat);
        }

   
#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif
    }
}

void IsoT_AIMc::GetSolution(Reservoir&             rs,
                            vector<OCP_DBL>& u,
                            const OCPControl&      ctrl) const
{
    const Domain&   domain = rs.domain;
    Bulk&           bk     = rs.bulk;
    const OCP_USI   nb     = bk.numBulk;
    const USI       np     = bk.numPhase;
    const USI       nc     = bk.numCom;
    const USI       row    = np * (nc + 1);
    const USI       col    = nc + 1;

    // Well first
    USI wId = bk.numBulkInterior * col;
    for (auto& wl : rs.allWells.wells) {
        if (wl.IsOpen()) {
            wl.SetBHP(wl.BHP() + u[wId]);
            wId += col;
        }
    }

    // Exchange Solution for ghost grid
    for (USI i = 0; i < domain.numRecvProc; i++) {
        const vector<OCP_USI>& rel = domain.recv_element_loc[i];
        MPI_Irecv(&u[rel[1] * col], (rel[2] - rel[1]) * col, MPI_DOUBLE, rel[0], 0, domain.myComm, &domain.recv_request[i]);
    }

    vector<vector<OCP_DBL>> send_buffer(domain.numSendProc);
    for (USI i = 0; i < domain.numSendProc; i++) {
        const vector<OCP_USI>& sel = domain.send_element_loc[i];
        vector<OCP_DBL>&       s   = send_buffer[i];
        s.resize(1 + (sel.size() - 1) * col);
        s[0] = sel[0];
        for (USI j = 1; j < sel.size(); j++) {
            const OCP_DBL* bId = u.data() + sel[j] * col;
            copy(bId, bId + col, &s[1 + (j - 1) * col]);
        }
        MPI_Isend(s.data() + 1, s.size() - 1, MPI_DOUBLE, s[0], 0, domain.myComm, &domain.send_request[i]);
    }

    // Bulk
    const OCP_DBL dSmaxlim = ctrl.ctrlNR.NRdSmax;
    // const OCP_DBL dPmaxlim = ctrl.ctrlNR.NRdPmax;

    vector<OCP_DBL> dtmp(row, 0);
    OCP_DBL         chopmin = 1;
    OCP_DBL         choptmp = 0;

    bk.dSNR       = bk.S;
    bk.NRphaseNum = bk.phaseNum;
    bk.NRdPmax    = 0;
    bk.NRdNmax    = 0;

    OCP_USI bId = 0;
    OCP_USI eId = bk.GetInteriorBulkNum();

    for (USI p = 0; p < 2; p++) {

        for (OCP_USI n = bId; n < eId; n++) {
            if (bk.bulkTypeAIM.IfIMPECbulk(n)) {
                // IMPEC Bulk
                // Pressure
                OCP_DBL dP   = u[n * col];
                bk.NRdPmax   = max(bk.NRdPmax, fabs(dP));
                bk.P[n]      += dP; // seems better
                bk.dPNR[n]   = dP;
                bk.NRstep[n] = 1;
                // Ni
                for (USI i = 0; i < nc; i++) {
                    bk.dNNR[n * nc + i] = u[n * col + 1 + i];
                    bk.Ni[n * nc + i] += bk.dNNR[n * nc + i];

                    // if (bk.Ni[n * nc + i] < 0 && bk.Ni[n * nc + i] > -1E-3) {
                    //     bk.Ni[n * nc + i] = 1E-20;
                    // }
                }
                // Pj
                for (USI j = 0; j < np; j++) {
                    bk.Pj[n * np + j] = bk.P[n] + bk.Pc[n * np + j];
                }
                continue;
            }

            chopmin = 1;
            // compute the chop
            fill(dtmp.begin(), dtmp.end(), 0.0);
            DaAxpby(row, col, 1, &bk.dSec_dPri[n * bk.maxLendSdP],
                u.data() + n * col, 1, dtmp.data());

            for (USI j = 0; j < np; j++) {
                choptmp = 1;
                if (fabs(dtmp[j]) > dSmaxlim) {
                    choptmp = dSmaxlim / fabs(dtmp[j]);
                }
                else if (bk.S[n * np + j] + dtmp[j] < 0.0) {
                    choptmp = 0.9 * bk.S[n * np + j] / fabs(dtmp[j]);
                }

                // if (fabs(S[n * np + j] - scm[j]) > TINY &&
                //     (S[n * np + j] - scm[j]) / (choptmp * dtmp[js]) < 0)
                //     choptmp *= min(1.0, -((S[n * np + j] - scm[j]) / (choptmp * dtmp[js])));

                chopmin = min(chopmin, choptmp);
            }

            // dS
            for (USI j = 0; j < np; j++) {
                bk.dSNRP[n * np + j] = chopmin * dtmp[j];
            }

            // dxij   ---- Compositional model only
            if (bk.IfUseEoS()) {
                USI js = np;
                if (bk.phaseNum[n] >= 3) {
                    OCP_USI bId = 0;
                    for (USI j = 0; j < 2; j++) {
                        bId = n * np * nc + j * nc;
                        for (USI i = 0; i < bk.numComH; i++) {
                            bk.xij[bId + i] += chopmin * dtmp[js];
                            js++;
                        }
                        js++;
                    }
                }
            }

            // dP
            OCP_DBL dP = u[n * col];
            if (fabs(bk.NRdPmax) < fabs(dP)) bk.NRdPmax = dP;
            bk.P[n] += dP; // seems better
            bk.dPNR[n] = dP;

            // dNi
            bk.NRstep[n] = chopmin;
            for (USI i = 0; i < nc; i++) {
                bk.dNNR[n * nc + i] = u[n * col + 1 + i] * chopmin;
                if (fabs(bk.NRdNmax) < fabs(bk.dNNR[n * nc + i]) / bk.Nt[n])
                    bk.NRdNmax = bk.dNNR[n * nc + i] / bk.Nt[n];

                bk.Ni[n * nc + i] += bk.dNNR[n * nc + i];

                // if (bk.Ni[n * nc + i] < 0 && bk.Ni[n * nc + i] > -1E-3) {
                //     bk.Ni[n * nc + i] = 1E-20;
                // }
            }
        }
        if (p == 0) {
            bId = eId;
            eId = nb;
            MPI_Waitall(domain.numRecvProc, domain.recv_request.data(), MPI_STATUS_IGNORE);
        }
        else {
            break;
        }
    }

    MPI_Waitall(domain.numSendProc, domain.send_request.data(), MPI_STATUS_IGNORE);
}

void IsoT_AIMc::ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.bulk.vj    = rs.bulk.lvj;
    rs.bulk.xijNR = rs.bulk.lxij;
    IsoT_FIM::ResetToLastTimeStep(rs, ctrl);

    // all FIM
    rs.bulk.bulkTypeAIM.Init(0);
    CalFlashI(rs.bulk);
    CalKrPcI(rs.bulk);

    //if (OCP_TRUE) {
    //    cout << "Rank " << CURRENT_RANK << "  " << rs.bulk.bulkTypeAIM.GetNumFIMBulk() * 1.0 / rs.bulk.numBulk * 100 << "% " << endl;
    //}
}

void IsoT_AIMc::UpdateLastTimeStep(Reservoir& rs) const
{
    IsoT_FIM::UpdateLastTimeStep(rs);

    rs.bulk.lvj   = rs.bulk.vj;
    rs.bulk.xijNR = rs.bulk.xij;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/01/2021      Create file                          */
/*  Chensong Zhang      Jan/08/2022      Update output                        */
/*----------------------------------------------------------------------------*/