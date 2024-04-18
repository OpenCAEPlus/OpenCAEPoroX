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

#include "IsoThermalMethod.hpp"

////////////////////////////////////////////
// IsothermalMethod
////////////////////////////////////////////


void IsothermalMethod::CalRock(Bulk& bk) const
{
    auto& bvs = bk.vs;
    for (OCP_USI n = 0; n < bvs.nb; n++) {
        auto ROCK = bk.ROCKm.GetROCK(n);

        ROCK->CalPoro(bvs.P[n], bvs.T[n], bvs.poroInit[n], BulkContent::rf);
        bvs.poro[n]   = ROCK->GetPoro();
        bvs.poroP[n]  = ROCK->GetdPorodP();
        bvs.rockVp[n] = bvs.v[n] * bvs.poro[n];
    }
}


void IsothermalMethod::ExchangeSolutionP(Reservoir& rs) const
{
    // Exchange Ghost P
    const Domain& domain = rs.domain;
    BulkVarSet&   bvs    = rs.bulk.vs;

    USI iter = 0;
    for (const auto& r : domain.recv_element_loc) {
        const auto& rv = r.second;
        MPI_Irecv(&bvs.P[rv[0]], rv[1] - rv[0], OCPMPI_DBL, r.first, 0, domain.global_comm, &domain.recv_request[iter]);
        iter++;
    }

    iter = 0;
    vector<vector<OCP_DBL>> send_buffer(domain.send_element_loc.size());
    for (const auto& s : domain.send_element_loc) {
        const auto& sv = s.second;
        auto&       sb = send_buffer[iter];
        sb.reserve(sv.size());
        for (const auto& sv1 : sv) {
            sb.push_back(bvs.P[sv1]);
        }
        MPI_Isend(sb.data(), sb.size(), OCPMPI_DBL, s.first, 0, domain.global_comm, &domain.send_request[iter]);
        iter++;
    }

    MPI_Waitall(iter, domain.send_request.data(), MPI_STATUS_IGNORE);
    MPI_Waitall(iter, domain.recv_request.data(), MPI_STATUS_IGNORE);
}


void IsothermalMethod::ExchangeSolutionNi(Reservoir& rs) const
{
    // Exchange Ghost Ni
    const Domain& domain = rs.domain;
    BulkVarSet&   bvs    = rs.bulk.vs;
    const USI     nc     = bvs.nc;

    USI iter = 0;
    for (const auto& r : domain.recv_element_loc) {
        const auto& rv = r.second;
        MPI_Irecv(&bvs.Ni[rv[0] * nc], (rv[1] - rv[0]) * nc, OCPMPI_DBL, r.first, 0, domain.global_comm, &domain.recv_request[iter]);
        iter++;
    }

    iter = 0;
    vector<vector<OCP_DBL>> send_buffer(domain.send_element_loc.size());
    for (const auto& s : domain.send_element_loc) {
        const auto& sv = s.second;
        auto&       sb = send_buffer[iter];
        sb.reserve(sv.size() * nc);
        for (const auto& sv1 : sv) {
            const OCP_DBL* bId = &bvs.Ni[0] + sv1 * nc;
            sb.insert(sb.end(), bId, bId + nc);
        }
        MPI_Isend(sb.data(), sb.size(), OCPMPI_DBL, s.first, 0, domain.global_comm, &domain.send_request[iter]);
        iter++;
    }

    MPI_Waitall(iter, domain.send_request.data(), MPI_STATUS_IGNORE);
    MPI_Waitall(iter, domain.recv_request.data(), MPI_STATUS_IGNORE);
}


void IsothermalMethod::SetWorkLS(const USI& w, const USI& i) 
{
    wls = w; 
    if (i > 0) preM = OCP_TRUE;
}



////////////////////////////////////////////
// IsoT_IMPEC
////////////////////////////////////////////

void IsoT_IMPEC::Setup(Reservoir& rs, const OCPControl& ctrl)
{
    // Allocate Memory of auxiliary variables for IMPEC
    AllocateReservoir(rs);
}

/// Initialize reservoir
void IsoT_IMPEC::InitReservoir(Reservoir& rs) const
{
    rs.bulk.Initialize(rs.domain);
    CalRock(rs.bulk);

    InitFlash(rs.bulk);
    CalKrPc(rs.bulk);

    rs.conn.CalFluxCoeff(rs.bulk);
    CalBulkFlux(rs);

    rs.allWells.InitBHP(rs.bulk);

    UpdateLastTimeStep(rs);
}

void IsoT_IMPEC::Prepare(Reservoir& rs, OCPControl& ctrl)
{
    rs.allWells.PrepareWell(rs.bulk);
    NR.CalCFL(rs, ctrl.time.GetCurrentDt(), OCP_TRUE);
    if (!NR.CheckPhysical(rs, { "CFL" }, ctrl.time.GetCurrentDt())) {
        ctrl.time.CutDt(NR);
    }
    NR.InitIter();
}

void IsoT_IMPEC::AssembleMat(LinearSystem&    ls,
                             const Reservoir& rs,
                             const OCP_DBL&   dt) const
{
    AssembleMatBulks(ls, rs, dt);
    AssembleMatWells(ls, rs, dt);

    rs.domain.SetNumActWellLocal(rs.GetNumOpenWell());
}

OCP_BOOL IsoT_IMPEC::SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl)
{
    GetWallTime timer;

    timer.Start();
    ls.AssembleMatLinearSolver();
    OCPTIME_CONVERT_MAT_FOR_LS_IF += timer.Stop();
    
    timer.Start();

    int status = ls.Solve();
    OCPTIME_LSOLVER += timer.Stop();

    if (status < 0) {
        ls.ClearData();
        ctrl.time.CutDt(-1.0);
        NR.ResetIter();
        return OCP_FALSE;
    }


    NR.UpdateIter(abs(status));

#ifdef DEBUG
    //OCP_INT myrank = rs.domain.global_rank;
    //ls.OutputLinearSystem("proc" + to_string(myrank) + "_testA_IMPEC.out", 
    //                      "proc" + to_string(myrank) + "_testb_IMPEC.out");
    //ls.OutputSolution("proc" + to_string(myrank) + "_testx_IMPEC.out");
    //MPI_Barrier(rs.domain.global_comm);
    //OCP_ABORT("Stop");
#endif // DEBUG

    timer.Start();
    GetSolution(rs, ls.GetSolution());
    OCPTIME_NRSTEP += timer.Stop();
    ls.ClearData();

    return OCP_TRUE;
}

OCP_BOOL IsoT_IMPEC::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{

    // First check : Pressure check
    if (!NR.CheckPhysical(rs, { "BulkP", "WellP" }, ctrl.time.GetCurrentDt())) {
        ctrl.time.CutDt(NR);
        NR.ResetIter();
        rs.bulk.vs.P = rs.bulk.vs.lP;
        return OCP_FALSE;
    }

    // Calculate Flux between bulks and between bulks and wells
    CalFlux(rs);  
    MassConserve(rs, ctrl.time.GetCurrentDt());

    // Second check : CFL check
    NR.CalCFL(rs, ctrl.time.GetCurrentDt(), OCP_TRUE);
    // Third check: Ni check

    if (!NR.CheckPhysical(rs, { "CFL","BulkNi" }, ctrl.time.GetCurrentDt())) {
        ctrl.time.CutDt(NR);
        ResetToLastTimeStep01(rs, ctrl);
        return OCP_FALSE;
    }

    CalRock(rs.bulk);
    CalFlash(rs.bulk);

    // Fouth check: Volume error check
    if (!NR.CheckPhysical(rs, { "BulkVe" }, ctrl.time.GetCurrentDt())) {
        ctrl.time.CutDt(NR);
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
    rs.CalIPRT(ctrl.time.GetCurrentDt());
    NR.CalMaxChangeTime(rs);
    ctrl.CalNextTimeStep(NR, {"dP", "dN", "dS", "eV"});
    UpdateLastTimeStep(rs);
}

void IsoT_IMPEC::AllocateReservoir(Reservoir& rs)
{
    Bulk&         bk  = rs.bulk;
    BulkVarSet&   bvs = bk.vs;
    const OCP_USI nb  = bvs.nb;
    const USI     np  = bvs.np;
    const USI     nc  = bvs.nc;

    // Rock
    bvs.poro.resize(nb);
    bvs.rockVp.resize(nb);

    bvs.lporo.resize(nb);
    bvs.lrockVp.resize(nb);

    // derivatives
    bvs.poroP.resize(nb);
    bvs.lporoP.resize(nb);

    // Fluid
    bvs.Nt.resize(nb);
    bvs.Ni.resize(nb * nc);
    bvs.vf.resize(nb);
    bvs.initT.resize(nb);
    bvs.T.resize(nb);
    bvs.P.resize(nb);
    bvs.Pb.resize(nb);
    bvs.Pj.resize(nb * np);
    bvs.Pc.resize(nb * np);
    bvs.phaseExist.resize(nb * np, OCP_FALSE);
    bvs.S.resize(nb * np);
    bvs.vj.resize(nb * np);
    bvs.xij.resize(nb * np * nc);
    bvs.rho.resize(nb * np);
    bvs.xi.resize(nb * np);
    bvs.mu.resize(nb * np);
    bvs.kr.resize(nb * np);

    bvs.lNt.resize(nb);
    bvs.lNi.resize(nb * nc);
    bvs.lvf.resize(nb);
    bvs.lT.resize(nb);
    bvs.lP.resize(nb);
    bvs.lPj.resize(nb * np);
    bvs.lPc.resize(nb * np);
    bvs.lphaseExist.resize(nb * np, OCP_FALSE);
    bvs.lS.resize(nb * np);
    bvs.vj.resize(nb * np);
    bvs.lxij.resize(nb * np * nc);
    bvs.lrho.resize(nb * np);
    bvs.lxi.resize(nb * np);
    bvs.lmu.resize(nb * np);
    bvs.lkr.resize(nb * np);

    // derivatives
    bvs.vfP.resize(nb);
    bvs.vfi.resize(nb * nc);

    bvs.lvfP.resize(nb);
    bvs.lvfi.resize(nb * nc);

    BulkConn& conn = rs.conn;

    conn.vs.upblock.resize(conn.numConn * np);
    conn.vs.dP.resize(conn.numConn * np);
    conn.vs.flux_vj.resize(conn.numConn * np);
    conn.vs.flux_ni.resize(conn.numConn * nc);

    conn.vs.lupblock.resize(conn.numConn * np);
    conn.vs.ldP.resize(conn.numConn * np);


    NR.Setup(bvs, rs.domain);
}


void IsoT_IMPEC::InitFlash(Bulk& bk) const
{
    BulkVarSet& bvs = bk.vs;

    for (OCP_USI n = 0; n < bvs.nb; n++) {

        auto PVT = bk.PVTm.GetPVT(n);
        PVT->InitFlashIMPEC(n, bvs);
        
        for (USI i = 0; i < bvs.nc; i++) {
            bvs.Ni[n * bvs.nc + i] = PVT->GetNi(i);
        }
        PassFlashValue(bk, n);
    }
}

void IsoT_IMPEC::CalFlash(Bulk& bk)
{
    const BulkVarSet& bvs = bk.vs;

    for (OCP_USI n = 0; n < bvs.nb; n++) {

        bk.PVTm.GetPVT(n)->FlashIMPEC(n, bvs);
        PassFlashValue(bk, n);
    }
}

void IsoT_IMPEC::PassFlashValue(Bulk& bk, const OCP_USI& n) const
{
    BulkVarSet& bvs = bk.vs;
    const auto  PVT = bk.PVTm.GetPVT(n);

    const USI     np     = bvs.np;
    const USI     nc     = bvs.nc;
    const OCP_USI bIdp   = n * np;

    bvs.Nt[n]       = PVT->GetNt();
    bvs.vf[n]       = PVT->GetVf();

    for (USI j = 0; j < np; j++) {
        // Important! Saturation must be passed no matter if the phase exists. This is
        // because it will be used to calculate relative permeability and capillary
        // pressure at each time step. Make sure that all saturations are updated at
        // each step!
        bvs.phaseExist[bIdp + j] = PVT->GetPhaseExist(j);
        bvs.S[bIdp + j]          = PVT->GetS(j);
        if (bvs.phaseExist[bIdp + j]) {
            for (USI i = 0; i < nc; i++) {
                bvs.xij[bIdp * nc + j * nc + i] = PVT->GetXij(j, i);
            }
            bvs.vj[bIdp + j]  = PVT->GetVj(j);
            bvs.rho[bIdp + j] = PVT->GetRho(j);
            bvs.xi[bIdp + j]  = PVT->GetXi(j);
            bvs.mu[bIdp + j]  = PVT->GetMu(j);
        }
    }

    bvs.vfP[n] = PVT->GetVfP();
    for (USI i = 0; i < nc; i++) {
        bvs.vfi[n * nc + i] = PVT->GetVfi(i);
    }
}

void IsoT_IMPEC::CalKrPc(Bulk& bk) const
{
    BulkVarSet& bvs = bk.vs;

    for (OCP_USI n = 0; n < bvs.nb; n++) {

        auto SAT = bk.SATm.GetSAT(n);

        OCP_USI bId = n * bvs.np;
        SAT->CalKrPc(n, &bvs.S[bId]);
        copy(SAT->GetKr().begin(), SAT->GetKr().end(), &bvs.kr[bId]);
        copy(SAT->GetPc().begin(), SAT->GetPc().end(), &bvs.Pc[bId]);
        for (USI j = 0; j < bvs.np; j++)
            bvs.Pj[n * bvs.np + j] = bvs.P[n] + bvs.Pc[n * bvs.np + j];
    }
}

void IsoT_IMPEC::CalFlux(Reservoir& rs) const
{
    CalBulkFlux(rs);
    rs.allWells.CalFlux(rs.bulk);
}

void IsoT_IMPEC::CalBulkFlux(Reservoir& rs) const
{
    const Bulk&       bk  = rs.bulk;
    const BulkVarSet& bvs = bk.vs;

    BulkConn&   conn = rs.conn;
    const USI   np   = bvs.np;
    const USI   nc   = bvs.nc;

    // calculate a step flux using iteratorConn
    BulkConnVarSet&   bcvs = conn.vs;

    for (OCP_USI c = 0; c < conn.numConn; c++) {

        auto Flux = conn.FLUXm.GetFlux(c);

        Flux->CalFlux(conn.iteratorConn[c], bk);
        copy(Flux->GetConvectUpblock().begin(), Flux->GetConvectUpblock().end(), &bcvs.upblock[c * np]);
        copy(Flux->GetConvectDP().begin(), Flux->GetConvectDP().end(), &bcvs.dP[c * np]);
        copy(Flux->GetConvectVj().begin(), Flux->GetConvectVj().end(), &bcvs.flux_vj[c * np]);
        copy(Flux->GetFluxNi().begin(), Flux->GetFluxNi().end(), &bcvs.flux_ni[c * nc]);
    }
}

void IsoT_IMPEC::MassConserve(Reservoir& rs, const OCP_DBL& dt) const
{

    BulkVarSet& bvs = rs.bulk.vs;

    // Bulk to Bulk
    const USI       nc   = bvs.nc;
    const BulkConn& conn = rs.conn;
    
    OCP_USI bId, eId;

    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId = conn.iteratorConn[c].BId();
        eId = conn.iteratorConn[c].EId();

        for (USI i = 0; i < nc; i++) {
            bvs.Ni[eId * nc + i] += dt * conn.vs.flux_ni[c * nc + i];
            bvs.Ni[bId * nc + i] -= dt * conn.vs.flux_ni[c * nc + i];
        }
    }

    // Well to Bulk
    for (auto& wl : rs.allWells.wells) {
        if (wl->IsOpen()) {
            for (USI p = 0; p < wl->PerfNum(); p++) {
                OCP_USI k = wl->PerfLocation(p);
                for (USI i = 0; i < nc; i++) {
                    bvs.Ni[k * nc + i] -= wl->PerfQi_lbmol(p, i) * dt;
                }
            }
        }
    }

    ExchangeSolutionNi(rs);
}

void IsoT_IMPEC::AssembleMatBulks(LinearSystem&    ls,
                                  const Reservoir& rs,
                                  const OCP_DBL&   dt) const
{
    const Bulk&       bk  = rs.bulk;
    const BulkVarSet& bvs = bk.vs;

    const USI       numWell = rs.GetNumOpenWell();
    const BulkConn& conn    = rs.conn;
    const OCP_USI   nbI     = bvs.nbI;

    ls.AddDim(nbI);

    // accumulate term
    OCP_DBL val, rhs;
    for (OCP_USI n = 0; n < nbI; n++) {
        bk.ACCm.GetAccumuTerm()->CalValRhsIMPEC(n, bvs, dt, val, rhs);
        ls.NewDiag(n, val);
        ls.AddRhs(n, rhs);
    }


    // flux term
    OCP_USI bId, eId;
    OCP_DBL valbb, rhsb, valee, rhse;

    // Be careful when first bulk has no neighbors!
    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId       = conn.iteratorConn[c].BId();
        eId       = conn.iteratorConn[c].EId();
        auto Flux = conn.FLUXm.GetFlux(c);

        Flux->AssembleMatIMPEC(conn.iteratorConn[c], c, conn.vs, bk);
        valbb  = dt * Flux->GetValbb();
        valee  = dt * Flux->GetValee();
        rhsb   = dt * Flux->GetRhsb();
        rhse   = dt * Flux->GetRhse();


        if (eId < nbI) {
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
        wl->AssembleMatIMPEC(ls, rs.bulk, dt);
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
    Bulk&       bk = rs.bulk;
    BulkVarSet& bvs = bk.vs;
    const OCP_USI nb     = bvs.nb;
    const USI     np     = bvs.np;
    const Domain& domain = rs.domain;

    // Well first
    USI wId = bvs.nbI;
    for (auto& wl : rs.allWells.wells) {
        wl->GetSolutionIMPEC(u, wId);
    }

    // Exchange Solution
    USI iter = 0;
    for (const auto& r : domain.recv_element_loc) {
        const auto& rv = r.second;
        MPI_Irecv(&u[rv[0]], rv[1] - rv[0], OCPMPI_DBL, r.first, 0, domain.global_comm, &domain.recv_request[iter]);
        iter++;
    }

    iter = 0;
    vector<vector<OCP_DBL>> send_buffer(domain.send_element_loc.size());
    for (const auto& s : domain.send_element_loc) {
        const auto& sv = s.second;
        auto&       sb = send_buffer[iter];
        sb.reserve(sv.size());
        for (const auto& sv1 : sv) {
            sb.push_back(u[sv1]);
        }
        MPI_Isend(sb.data(), sb.size(), OCPMPI_DBL, s.first, 0, domain.global_comm, &domain.send_request[iter]);
        iter++;
    }

    // Bulk
    // interior first, ghost second
    OCP_USI bId = 0;
    OCP_USI eId = bvs.nbI;
    for (USI p = 0; p < 2; p++) {

        for (OCP_USI n = bId; n < eId; n++) {
            bvs.P[n] = u[n];
            for (USI j = 0; j < np; j++) {
                bvs.Pj[n * np + j] = bvs.P[n] + bvs.Pc[n * np + j];
            }
        }

        if (p == 0) {
            bId = eId;
            eId = nb;
            MPI_Waitall(iter, domain.recv_request.data(), MPI_STATUS_IGNORE);
        }
        else {
            break;
        }
    }
    MPI_Waitall(iter, domain.send_request.data(), MPI_STATUS_IGNORE);
}


void IsoT_IMPEC::ResetToLastTimeStep01(Reservoir& rs, OCPControl& ctrl)
{

    // Bulk
    rs.bulk.vs.P  = rs.bulk.vs.lP;
    rs.bulk.vs.Ni = rs.bulk.vs.lNi;
    rs.bulk.vs.Pj = rs.bulk.vs.lPj;
    // Bulk Conn
    rs.conn.vs.upblock      = rs.conn.vs.lupblock;
    rs.conn.vs.dP           = rs.conn.vs.ldP;

    // Iters
    NR.ResetIter();
}

void IsoT_IMPEC::ResetToLastTimeStep02(Reservoir& rs, OCPControl& ctrl)
{
    Bulk&       bk  = rs.bulk;
    BulkVarSet& bvs = bk.vs;
    // Rock
    bvs.rockVp = bvs.lrockVp;
    bvs.poro   = bvs.lporo;
    bvs.poroP  = bvs.lporoP;

    // Fluid
    bvs.P          = bvs.lP;
    bvs.Nt         = bvs.lNt;
    bvs.Ni         = bvs.lNi;
    bvs.vf         = bvs.lvf;
    bvs.Pj         = bvs.lPj;
    bvs.phaseExist = bvs.lphaseExist;
    bvs.S          = bvs.lS;
    bvs.vj         = bvs.lvj;
    bvs.xij        = bvs.lxij;
    bvs.rho        = bvs.lrho;
    bvs.xi         = bvs.lxi;
    bvs.mu         = bvs.lmu;

    // derivatives
    bvs.vfP = bvs.lvfP;
    bvs.vfi = bvs.lvfi;

    // Bulk Conn
    rs.conn.vs.upblock      = rs.conn.vs.lupblock;
    rs.conn.vs.dP           = rs.conn.vs.ldP;

    // boundary
    rs.bulk.BOUNDm.ResetToLastTimeStep();
    // Optional Features 
    rs.bulk.optMs.ResetToLastTimeStep();
    rs.conn.optMs.ResetToLastTimeStep();

    // Iters
    NR.ResetIter();
}

void IsoT_IMPEC::UpdateLastTimeStep(Reservoir& rs) const
{

    Bulk&       bk  = rs.bulk;
    BulkVarSet& bvs = bk.vs;

    // Rock
    bvs.lporo   = bvs.poro;
    bvs.lporoP  = bvs.poroP;
    bvs.lrockVp = bvs.rockVp;

    // Fluid
    bvs.lNt         = bvs.Nt;
    bvs.lNi         = bvs.Ni;
    bvs.lvf         = bvs.vf;
    bvs.lP          = bvs.P;
    bvs.lT          = bvs.T;
    bvs.lPj         = bvs.Pj;
    bvs.lPc         = bvs.Pc;
    bvs.lphaseExist = bvs.phaseExist;
    bvs.lS          = bvs.S;
    bvs.lvj         = bvs.vj;
    bvs.lxij        = bvs.xij;
    bvs.lrho        = bvs.rho;
    bvs.lxi         = bvs.xi;
    bvs.lmu         = bvs.mu;
    bvs.lkr         = bvs.kr;

    // derivatives
    bvs.lvfP = bvs.vfP;
    bvs.lvfi = bvs.vfi;

    // boundary
    rs.bulk.BOUNDm.UpdateLastTimeStep();

    BulkConn& conn = rs.conn;

    conn.vs.lupblock    = conn.vs.upblock;
    conn.vs.ldP         = conn.vs.dP;

    rs.allWells.UpdateLastTimeStep();

    rs.bulk.optMs.UpdateLastTimeStep();
    rs.conn.optMs.UpdateLastTimeStep();
}

////////////////////////////////////////////
// IsoT_FIM
////////////////////////////////////////////

void IsoT_FIM::Setup(Reservoir& rs, const OCPControl& ctrl)
{
    // Allocate memory for reservoir
    AllocateReservoir(rs);
}

void IsoT_FIM::InitReservoir(Reservoir& rs)
{
    // Calculate initial bulk pressure and temperature and water saturation
    rs.bulk.Initialize(rs.domain);
    // Initialize rock property
    CalRock(rs.bulk);
    // Initialize fluid properties
    InitFlash(rs.bulk);
    CalKrPc(rs.bulk);
    rs.conn.CalFluxCoeff(rs.bulk);
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
    CalInitRes(rs, dt);
    NR.InitStep(rs.bulk.GetVarSet());
    NR.InitIter();
}

void IsoT_FIM::AssembleMat(LinearSystem&    ls,
                           const Reservoir& rs,
                           const OCP_DBL&   dt) const
{
    // Assemble matrix
    AssembleMatBulks(ls, rs, dt);
    AssembleMatWells(ls, rs, dt);
    // Assemble rhs -- from residual
    ls.CopyRhs(NR.res.resAbs);

    rs.domain.SetNumActWellLocal(rs.GetNumOpenWell());
}

OCP_BOOL IsoT_FIM::SolveLinearSystem(LinearSystem& ls,
                                 Reservoir&    rs,
                                 OCPControl&   ctrl)
{
    GetWallTime timer;
    timer.Start();
    ls.AssembleMatLinearSolver();
    OCPTIME_CONVERT_MAT_FOR_LS_IF += timer.Stop();


    // Solve linear system
    //ls.OutputLinearSystem("proc" + to_string(CURRENT_RANK) + "_A_ddm.out",
    //    "proc" + to_string(CURRENT_RANK) + "_b_ddm.out");

    timer.Start();
    int status = ls.Solve();
    // Record time, iterations
    OCPTIME_LSOLVER += timer.Stop();

    if (status < 0) {
        ls.ClearData();
        ctrl.time.CutDt(-1.0);
        ResetToLastTimeStep(rs, ctrl);     
        return OCP_FALSE;
    }

    NR.UpdateIter(abs(status));

    // ls.OutputSolution("proc" + to_string(CURRENT_RANK) + "_x_ddm.out");

#ifdef DEBUG
     //// Output A, b, x
     //ls.OutputLinearSystem("proc" + to_string(CURRENT_RANK) + "_testA_FIM.out",
     //                      "proc" + to_string(CURRENT_RANK) + "_testb_FIM.out");
     //MPI_Barrier(rs.domain.global_comm);

     //ls.OutputSolution("proc" + to_string(CURRENT_RANK) + "_testx_FIM.out");
     //OCP_ABORT("Stop");
#endif // DEBUG


    //MPI_Barrier(rs.domain.global_comm);
    //OCP_ABORT("Stop");

    // Get solution from linear system to Reservoir
    timer.Start();
    GetSolution(rs, ls.GetSolution(), ctrl.NR);
    OCPTIME_NRSTEP += timer.Stop();
    // rs.PrintSolFIM(ctrl.workDir + "testPNi.out");
    ls.ClearData();

    return OCP_TRUE;
}

OCP_BOOL IsoT_FIM::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    if (!NR.CheckPhysical(rs, { "BulkNi", "BulkP" }, ctrl.time.GetCurrentDt())) {
        ctrl.time.CutDt(NR);
        ResetToLastTimeStep(rs, ctrl);
        return OCP_FALSE;
    }

    // Update fluid property
    CalFlash(rs.bulk);
    CalKrPc(rs.bulk);
    // Update rock property
    CalRock(rs.bulk);
    // Update well property
    rs.allWells.CalFlux(rs.bulk);

    // Update residual
    CalRes(rs, ctrl.time.GetCurrentDt());

    return OCP_TRUE;
}

OCP_BOOL IsoT_FIM::FinishNR(Reservoir& rs, OCPControl& ctrl)
{

    if (ctrl.time.GetCurrentDt() < 5E-2 + TINY) {
        return OCP_TRUE;
    }

    NR.CalMaxChangeNR(rs);
    const OCPNRStateC conflag = ctrl.CheckConverge(NR, { "res", "d" });

    if (conflag == OCPNRStateC::converge) {
        if (!NR.CheckPhysical(rs, { "WellP" }, ctrl.time.GetCurrentDt())) {
            ctrl.time.CutDt(NR);
            ResetToLastTimeStep(rs, ctrl);
            return OCP_FALSE;
        } else {
            return OCP_TRUE;
        }
    } else if (conflag == OCPNRStateC::not_converge) {
        ctrl.time.CutDt();
        ResetToLastTimeStep(rs, ctrl);
        return OCP_FALSE;
    } else {
        return OCP_FALSE;
    }
}

void IsoT_FIM::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.CalIPRT(ctrl.time.GetCurrentDt());
    NR.CalMaxChangeTime(rs);
    ctrl.CalNextTimeStep(NR, {"dP", "dS", "iter"});
    UpdateLastTimeStep(rs);
}

void IsoT_FIM::AllocateReservoir(Reservoir& rs)
{

    Bulk&       bk  = rs.bulk;
    BulkVarSet& bvs = bk.vs;

    const OCP_USI nb = bvs.nb;
    const USI     np = bvs.np;
    const USI     nc = bvs.nc;

    // Rock
    bvs.poro.resize(nb);
    bvs.rockVp.resize(nb);

    bvs.lporo.resize(nb);
    bvs.lrockVp.resize(nb);

    // derivatives
    bvs.poroP.resize(nb);
    bvs.lporoP.resize(nb);

    // Fluid
    bvs.Nt.resize(nb);
    bvs.Ni.resize(nb * nc);
    bvs.vf.resize(nb);
    bvs.initT.resize(nb);
    bvs.T.resize(nb);
    bvs.P.resize(nb);
    bvs.Pb.resize(nb);
    bvs.Pj.resize(nb * np);
    bvs.Pc.resize(nb * np);
    bvs.phaseExist.resize(nb * np, OCP_FALSE);
    bvs.S.resize(nb * np);
    bvs.xij.resize(nb * np * nc);
    bvs.rho.resize(nb * np);
    bvs.xi.resize(nb * np);
    bvs.mu.resize(nb * np);
    bvs.kr.resize(nb * np);

    bvs.lNt.resize(nb);
    bvs.lNi.resize(nb * nc);
    bvs.lvf.resize(nb);
    bvs.lT.resize(nb);
    bvs.lP.resize(nb);
    bvs.lPj.resize(nb * np);
    bvs.lPc.resize(nb * np);
    bvs.lphaseExist.resize(nb * np, OCP_FALSE);
    bvs.lS.resize(nb * np);
    bvs.lxij.resize(nb * np * nc);
    bvs.lrho.resize(nb * np);
    bvs.lxi.resize(nb * np);
    bvs.lmu.resize(nb * np);
    bvs.lkr.resize(nb * np);

    // derivatives
    bvs.vfP.resize(nb);
    bvs.vfi.resize(nb * nc);
    bvs.rhoP.resize(nb * np);
    bvs.rhox.resize(nb * nc * np);
    bvs.xiP.resize(nb * np);
    bvs.xix.resize(nb * nc * np);
    bvs.muP.resize(nb * np);
    bvs.mux.resize(nb * nc * np);
    bvs.dPcdS.resize(nb * np * np);
    bvs.dKrdS.resize(nb * np * np);

    bvs.lvfP.resize(nb);
    bvs.lvfi.resize(nb * nc);
    bvs.lrhoP.resize(nb * np);
    bvs.lrhox.resize(nb * nc * np);
    bvs.lxiP.resize(nb * np);
    bvs.lxix.resize(nb * nc * np);
    bvs.lmuP.resize(nb * np);
    bvs.lmux.resize(nb * nc * np);
    bvs.ldPcdS.resize(nb * np * np);
    bvs.ldKrdS.resize(nb * np * np);

    // FIM-Specified
    bvs.lendSdP = (nc + 1) * (nc + 1) * np;
    bvs.dSec_dPri.resize(nb * bvs.lendSdP);

    bvs.ldSec_dPri.resize(nb * bvs.lendSdP);

    // BulkConn
    BulkConn& conn = rs.conn;

    conn.vs.upblock.resize(conn.numConn* np);
    conn.vs.dP.resize(conn.numConn* np);
    conn.vs.flux_vj.resize(conn.numConn* np);


    // Allocate Residual
    NR.Setup(OCP_FALSE, bvs, rs.allWells.numWell, rs.domain);
}


void IsoT_FIM::InitFlash(Bulk& bk)
{
    BulkVarSet& bvs = bk.vs;

    for (OCP_USI n = 0; n < bvs.nb; n++) {
        auto PVT = bk.PVTm.GetPVT(n);

        PVT->InitFlashFIM(n, bvs);
        for (USI i = 0; i < bvs.nc; i++) {
            bvs.Ni[n * bvs.nc + i] = PVT->GetNi(i);
        }
        PassFlashValue(bk, n);
    }
}

void IsoT_FIM::CalFlash(Bulk& bk)
{
    const BulkVarSet& bvs = bk.vs;

    for (OCP_USI n = 0; n < bvs.nb; n++) {

        bk.PVTm.GetPVT(n)->FlashFIM(n, bvs);
        PassFlashValue(bk, n);
    }
}

void IsoT_FIM::PassFlashValue(Bulk& bk, const OCP_USI& n)
{
    auto&         bvs = bk.vs;
    const auto    PVT = bk.PVTm.GetPVT(n);

    const auto    np     = bvs.np;
    const auto    nc     = bvs.nc;
    const OCP_USI bIdp   = n * np;

    bvs.Nt[n]       = PVT->GetNt();
    bvs.vf[n]       = PVT->GetVf();

    for (USI j = 0; j < np; j++) {
        // Important! Saturation must be passed no matter if the phase exists. This is
        // because it will be used to calculate relative permeability and capillary
        // pressure at each time step. Make sure that all saturations are updated at
        // each step!
        bvs.S[bIdp + j]          = PVT->GetS(j);
        bvs.phaseExist[bIdp + j] = PVT->GetPhaseExist(j);
        if (bvs.phaseExist[bIdp + j]) {
            bvs.rho[bIdp + j] = PVT->GetRho(j);
            bvs.xi[bIdp + j]  = PVT->GetXi(j);
            bvs.mu[bIdp + j]  = PVT->GetMu(j);

            // Derivatives
            bvs.rhoP[bIdp + j] = PVT->GetRhoP(j);
            bvs.xiP[bIdp + j]  = PVT->GetXiP(j);
            bvs.muP[bIdp + j]  = PVT->GetMuP(j);

            for (USI i = 0; i < nc; i++) {
                bvs.xij[bIdp * nc + j * nc + i]  = PVT->GetXij(j, i);
                bvs.rhox[bIdp * nc + j * nc + i] = PVT->GetRhoX(j, i);
                bvs.xix[bIdp * nc + j * nc + i]  = PVT->GetXiX(j, i);
                bvs.mux[bIdp * nc + j * nc + i]  = PVT->GetMuX(j, i);
            }
        }
    }

    bvs.vfP[n] = PVT->GetVfP();
    for (USI i = 0; i < nc; i++) {
        bvs.vfi[n * nc + i] = PVT->GetVfi(i);
    }

    copy(PVT->GetDXsDXp().begin(), PVT->GetDXsDXp().end(), &bvs.dSec_dPri[n * bvs.lendSdP]);
}

void IsoT_FIM::CalKrPc(Bulk& bk) const
{
    BulkVarSet& bvs = bk.vs;
    const USI&  np  = bvs.np;
    for (OCP_USI n = 0; n < bvs.nb; n++) {
        auto SAT = bk.SATm.GetSAT(n);

        const OCP_USI bId = n * np;
        SAT->CalKrPcFIM(n, &bvs.S[bId]);
        copy(SAT->GetKr().begin(), SAT->GetKr().end(), &bvs.kr[bId]);
        copy(SAT->GetPc().begin(), SAT->GetPc().end(), &bvs.Pc[bId]);
        copy(SAT->GetdKrdS().begin(), SAT->GetdKrdS().end(), &bvs.dKrdS[bId * np]);
        copy(SAT->GetdPcdS().begin(), SAT->GetdPcdS().end(), &bvs.dPcdS[bId * np]);
        for (USI j = 0; j < np; j++) bvs.Pj[bId + j] = bvs.P[n] + bvs.Pc[bId + j];
    }
}

void IsoT_FIM::CalRes(Reservoir& rs, const OCP_DBL& dt, const OCP_BOOL& initRes0)
{
    const Bulk&       bk  = rs.bulk;
    const BulkVarSet& bvs = bk.vs;

    const USI nb  = bvs.nbI;
    const USI np  = bvs.np;
    const USI nc  = bvs.nc;
    const USI len = nc + 1;

    OCPNRresidual&   res = NR.res;
    
    res.SetZero();
  
    // Accumalation Term
    for (OCP_USI n = 0; n < nb; n++) {
        const vector<OCP_DBL>& r = bk.ACCm.GetAccumuTerm()->CalResFIM(n, bvs, dt);
        copy(r.begin(), r.end(), &res.resAbs[n * len]);
    }

    // Flux Term
    OCP_USI         bId, eId;
    BulkConn&       conn = rs.conn;
    BulkConnVarSet& bcvs = conn.vs;
    for (OCP_USI c = 0; c < conn.numConn; c++) {

        bId       = conn.iteratorConn[c].BId();
        eId       = conn.iteratorConn[c].EId();
        auto Flux = conn.FLUXm.GetFlux(c);

        Flux->CalFlux(conn.iteratorConn[c], bk);
        copy(Flux->GetConvectUpblock().begin(), Flux->GetConvectUpblock().end(), &bcvs.upblock[c * np]);
        copy(Flux->GetConvectDP().begin(), Flux->GetConvectDP().end(), &bcvs.dP[c * np]);
        copy(Flux->GetConvectVj().begin(), Flux->GetConvectVj().end(), &bcvs.flux_vj[c * np]);
               
        if (eId < nb) {
            for (USI i = 0; i < nc; i++) {               
                res.resAbs[bId * len + 1 + i] += dt * Flux->GetFluxNi()[i];
                res.resAbs[eId * len + 1 + i] -= dt * Flux->GetFluxNi()[i];
            }
        }
        else {
            for (USI i = 0; i < nc; i++) {
                res.resAbs[bId * len + 1 + i] += dt * Flux->GetFluxNi()[i];
            }
        }
    }

    // Well to Bulk, Well
    USI wId = nb * len;
    for (const auto& wl : rs.allWells.wells) {
        wl->CalResFIM(wId, res, bk, dt);
    }

    // Calculate RelRes
    OCP_DBL tmp;
    for (OCP_USI n = 0; n < nb; n++) {

        for (USI i = 0; i < len; i++) {
            tmp = fabs(res.resAbs[n * len + i] / bvs.rockVp[n]);
            if (res.maxRelRes_V < tmp) {
                res.maxRelRes_V = tmp;
                res.maxId_V     = n;
            }
            res.resRelV[n] += tmp * tmp;
        }
        res.resRelV[n] = sqrt(res.resRelV[n]);

        for (USI i = 1; i < len; i++) {
            tmp = fabs(res.resAbs[n * len + i] / bvs.Nt[n]);
            if (res.maxRelRes_N < tmp) {
                res.maxRelRes_N = tmp;
                res.maxId_N     = n;
            }
            res.resRelN[n] += tmp * tmp;
        }
        res.resRelN[n] = sqrt(res.resRelN[n]);
    }

    Dscalar(res.resAbs.size(), -1.0, res.resAbs.data());

    if (initRes0) {
        GetWallTime timer;
        timer.Start();

        MPI_Allreduce(&res.maxRelRes_V, &res.maxRelRes0_V, 1, OCPMPI_DBL, MPI_MIN, rs.domain.global_comm);

        OCPTIME_COMM_COLLECTIVE += timer.Stop();
    }

}

void IsoT_FIM::AssembleMatBulks(LinearSystem&    ls,
                                const Reservoir& rs,
                                const OCP_DBL&   dt) const
{
    const Bulk& bk  = rs.bulk;
    const BulkVarSet& bvs = bk.vs;

    const USI numWell = rs.GetNumOpenWell();

    const BulkConn& conn   = rs.conn;
    const OCP_USI   nbI    = bvs.nbI;
    const USI       np     = bvs.np;
    const USI       nc     = bvs.nc;
    const USI       ncol   = nc + 1;
    const USI       ncol2  = np * nc + np;
    const USI       bsize  = ncol * ncol;
    const USI       bsize2 = ncol * ncol2;

    ls.AddDim(nbI);


    // Accumulation term
    vector<OCP_DBL> bmat(bsize, 0);
    for (OCP_USI n = 0; n < nbI; n++) {
        ls.NewDiag(n, bk.ACCm.GetAccumuTerm()->CaldFdXpFIM(n, bvs, dt));
    }

    // flux term  
    OCP_USI  bId, eId;
    for (OCP_USI c = 0; c < conn.numConn; c++) {

        bId       = conn.iteratorConn[c].BId();
        eId       = conn.iteratorConn[c].EId();
        auto Flux = conn.FLUXm.GetFlux(c);

        Flux->AssembleMatFIM(conn.iteratorConn[c], c, conn.vs, bk);
        
        bmat = Flux->GetdFdXpB();
        DaABpbC(ncol, ncol, ncol2, 1, Flux->GetdFdXsB().data(), &bvs.dSec_dPri[bId * bsize2], 1,
            bmat.data());
        Dscalar(bsize, dt, bmat.data());

        // Assemble
        // Begin - Begin -- add
        ls.AddDiag(bId, bmat);
        // End - Begin -- insert
        if (eId < nbI) {
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
        bmat = Flux->GetdFdXpE();
        DaABpbC(ncol, ncol, ncol2, 1, Flux->GetdFdXsE().data(), &bvs.dSec_dPri[eId * bsize2], 1,
            bmat.data());
        Dscalar(bsize, dt, bmat.data());
        
        if (eId < nbI) {
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
        wl->AssembleMatFIM(ls, rs.bulk, dt);
    }
}


void IsoT_FIM::GetSolution(Reservoir&        rs,
                           vector<OCP_DBL>&  u,
                           const ControlNR& ctrlNR)
{
    const auto& domain = rs.domain;
    auto&       bk     = rs.bulk;
    auto&       bvs    = bk.vs;
    const auto  nb     = bvs.nb;
    const auto  np     = bvs.np;
    const auto  nc     = bvs.nc;
    const auto  row    = np * (nc + 1);
    const auto  col    = nc + 1;

    // Well first
    USI wId = bvs.nbI * col;
    for (auto& wl : rs.allWells.wells) {
        wl->GetSolutionFIM(u, wId);
    }

    GetWallTime timerT;         ///< total timer
    GetWallTime timerC;         ///< calculation timer
    OCP_DBL     time_cal = 0;   ///< calculation time
    timerT.Start();
    

    USI iter = 0;
    for (const auto& r : domain.recv_element_loc) {
        const auto& rv = r.second;
        MPI_Irecv(&u[rv[0] * col], (rv[1] - rv[0]) * col, OCPMPI_DBL, r.first, 0, domain.global_comm, &domain.recv_request[iter]);
        iter++;
    }

    iter = 0;
    vector<vector<OCP_DBL>> send_buffer(domain.send_element_loc.size());
    for (const auto& s : domain.send_element_loc) {
        const auto& sv = s.second;
        auto&       sb = send_buffer[iter];
        sb.reserve(sv.size() * col);
        for (const auto& sv1 : sv) {
            const OCP_DBL* bId = u.data() + sv1 * col;
            sb.insert(sb.end(), bId, bId + col);
        }
        MPI_Isend(sb.data(), sb.size(), OCPMPI_DBL, s.first, 0, domain.global_comm, &domain.send_request[iter]);
        iter++;
    }


    // Bulk
    const OCP_DBL dSmaxlim = ctrlNR.DSmax();
    const OCP_DBL dPmaxlim = ctrlNR.DPmax();

    vector<OCP_DBL> dtmp(row, 0);
    OCP_DBL         chopmin = 1;
    OCP_DBL         choptmp = 0;

    OCP_USI bId = 0;
    OCP_USI eId = bvs.nbI;

    // interior first, ghost second
    for (USI p = 0; p < 2; p++) {
  
        timerC.Start();

		for (OCP_USI n = bId; n < eId; n++) {
			// const vector<OCP_DBL>& scm = satcm[SATNUM[n]];

			chopmin = 1;
			// compute the chop
			fill(dtmp.begin(), dtmp.end(), 0.0);

			OCP_aAxpby(row, col, static_cast<OCP_DBL>(1.0), &bvs.dSec_dPri[n * bvs.lendSdP], u.data() + n * col, static_cast<OCP_DBL>(1.0), dtmp.data());

			for (USI j = 0; j < np; j++) {
				choptmp = 1;
				if (fabs(dtmp[j]) > dSmaxlim) {
					choptmp = dSmaxlim / fabs(dtmp[j]);
				}
				else if (bvs.S[n * np + j] + dtmp[j] < 0.0) {
					choptmp = 0.9 * bvs.S[n * np + j] / fabs(dtmp[j]);
				}
				// if (fabs(S[n_np_j] - scm[j]) > TINY &&
				//     (S[n_np_j] - scm[j]) / (choptmp * dtmp[js]) < 0)
				//     choptmp *= min(1.0, -((S[n_np_j] - scm[j]) / (choptmp * dtmp[js])));
				chopmin = min(chopmin, choptmp);
			}

			// dS
			for (USI j = 0; j < np; j++) {
				bvs.S[n * np + j] += chopmin * dtmp[j];
			}

			// dxij
			USI js = np;
			for (USI j = 0; j < np; j++) {
				for (USI i = 0; i < bvs.nc; i++) {
					bvs.xij[(n * np + j) * nc + i] += chopmin * dtmp[js];
					js++;
				}
			}

			// dP
			//choptmp = dPmaxlim / fabs(u[n * col]);
			//chopmin = min(chopmin, choptmp);
			bvs.P[n] += u[n * col]; // seems better

			// dNi
			for (USI i = 0; i < nc; i++) {
				bvs.Ni[n * nc + i] += chopmin * u[n * col + 1 + i];

				// if (bvs.Ni[n * nc + i] < 0 && bvs.Ni[n * nc + i] > -1E-3) {
				//     bvs.Ni[n * nc + i] = 1E-20;
				// }
			}
		}

        time_cal += timerC.Stop();

        if (p == 0) {
            bId = eId;
            eId = nb;
            MPI_Waitall(iter, domain.recv_request.data(), MPI_STATUS_IGNORE);
        }
        else {
            break;
        }        
    }

    MPI_Waitall(iter, domain.send_request.data(), MPI_STATUS_IGNORE);

    OCPTIME_COMM_P2P += (timerT.Stop() - time_cal);
    OCPTIME_NRSTEPC  += time_cal;
}

void IsoT_FIM::ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl)
{
    Bulk&       bk = rs.bulk;
    BulkVarSet& bvs = bk.vs;

    // Rock
    bvs.poro   = bvs.lporo;
    bvs.poroP  = bvs.lporoP;
    bvs.rockVp = bvs.lrockVp;
    // Fluid
    bvs.Nt         = bvs.lNt;
    bvs.Ni         = bvs.lNi;
    bvs.vf         = bvs.lvf;
    bvs.P          = bvs.lP;
    bvs.Pj         = bvs.lPj;
    bvs.Pc         = bvs.lPc;
    bvs.phaseExist = bvs.lphaseExist;
    bvs.S          = bvs.lS;
    bvs.xij        = bvs.lxij;
    bvs.rho        = bvs.lrho;
    bvs.xi         = bvs.lxi;
    bvs.mu         = bvs.lmu;
    bvs.kr         = bvs.lkr;
    // derivatives
    bvs.vfP     = bvs.lvfP;
    bvs.vfi     = bvs.lvfi;
    bvs.rhoP    = bvs.lrhoP;
    bvs.rhox    = bvs.lrhox;
    bvs.xiP     = bvs.lxiP;
    bvs.xix     = bvs.lxix;
    bvs.muP     = bvs.lmuP;
    bvs.mux     = bvs.lmux;
    bvs.dPcdS   = bvs.ldPcdS;
    bvs.dKrdS   = bvs.ldKrdS;
    // FIM-Specified
    bvs.dSec_dPri    = bvs.ldSec_dPri;

    // Wells
    rs.allWells.ResetToLastTimeStep(bk);

    // boundary
    rs.bulk.BOUNDm.ResetToLastTimeStep();
    // Optional Features
    rs.bulk.optMs.ResetToLastTimeStep();
    rs.conn.optMs.ResetToLastTimeStep();

    // Residual
    CalInitRes(rs, ctrl.time.GetCurrentDt());

    NR.InitStep(rs.bulk.GetVarSet());
    NR.ResetIter();
}

void IsoT_FIM::UpdateLastTimeStep(Reservoir& rs) const
{
    BulkVarSet& bvs = rs.bulk.vs;
    // Rock
    bvs.lporo   = bvs.poro;
    bvs.lporoP  = bvs.poroP;
    bvs.lrockVp = bvs.rockVp;

    // Fluid
    bvs.lNt         = bvs.Nt;
    bvs.lNi         = bvs.Ni;
    bvs.lvf         = bvs.vf;
    bvs.lP          = bvs.P;
    bvs.lT          = bvs.T;
    bvs.lPj         = bvs.Pj;
    bvs.lPc         = bvs.Pc;
    bvs.lphaseExist = bvs.phaseExist;
    bvs.lS          = bvs.S;
    bvs.lxij        = bvs.xij;
    bvs.lrho        = bvs.rho;
    bvs.lxi         = bvs.xi;
    bvs.lmu         = bvs.mu;
    bvs.lkr         = bvs.kr;

    // derivatives
    bvs.lvfP     = bvs.vfP;
    bvs.lvfi     = bvs.vfi;
    bvs.lrhoP    = bvs.rhoP;
    bvs.lrhox    = bvs.rhox;
    bvs.lxiP     = bvs.xiP;
    bvs.lxix     = bvs.xix;
    bvs.lmuP     = bvs.muP;
    bvs.lmux     = bvs.mux;
    bvs.ldPcdS = bvs.dPcdS;
    bvs.ldKrdS  = bvs.dKrdS;
    // boundary
    rs.bulk.BOUNDm.UpdateLastTimeStep();
    // FIM-Specified
    bvs.ldSec_dPri    = bvs.dSec_dPri;

    rs.allWells.UpdateLastTimeStep();

    rs.bulk.optMs.UpdateLastTimeStep();
    rs.conn.optMs.UpdateLastTimeStep();
}

////////////////////////////////////////////
// IsoT_AIMc
////////////////////////////////////////////

void IsoT_AIMc::Setup(Reservoir& rs, const OCPControl& ctrl)
{
    // Allocate Bulk and BulkConn Memory
    AllocateReservoir(rs);
    // Setup neighbor
    SetupNeighbor(rs);
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


void IsoT_AIMc::InitReservoir(Reservoir& rs)
{
    rs.bulk.Initialize(rs.domain);
    CalRock(rs.bulk);

    IsoT_IMPEC::InitFlash(rs.bulk);
    IsoT_IMPEC::CalKrPc(rs.bulk);
    rs.conn.CalFluxCoeff(rs.bulk);
    rs.allWells.InitBHP(rs.bulk);

    UpdateLastTimeStep(rs);
}

void IsoT_AIMc::Prepare(Reservoir& rs, const OCP_DBL& dt)
{
    rs.allWells.PrepareWell(rs.bulk);
    CalInitRes(rs, dt);

    // Set FIM Bulk
    NR.CalCFL(rs, dt, OCP_FALSE);
    rs.allWells.SetupWellBulk(rs.bulk);
    SetFIMBulk(rs);
    //  Calculate FIM Bulk properties
    CalFlashI(rs.bulk);
    CalKrPcI(rs.bulk);

    UpdateLastTimeStep(rs);
    NR.InitStep(rs.bulk.GetVarSet());
    NR.InitIter();
}

void IsoT_AIMc::AssembleMat(LinearSystem&    ls,
                            const Reservoir& rs,
                            const OCP_DBL&   dt) const
{
    AssembleMatBulks(ls, rs, dt);
    IsoT_FIM::AssembleMatWells(ls, rs, dt);
    ls.CopyRhs(NR.res.resAbs);

    rs.domain.SetNumActWellLocal(rs.GetNumOpenWell());
}

OCP_BOOL IsoT_AIMc::SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl)
{
    GetWallTime timer;

    timer.Start();
    ls.AssembleMatLinearSolver();
    OCPTIME_CONVERT_MAT_FOR_LS_IF += timer.Stop();
    timer.Start();
    int status = ls.Solve();
    OCPTIME_LSOLVER += timer.Stop();

    if (status < 0) {
        ls.ClearData();
        ctrl.time.CutDt(-1.0);
        ResetToLastTimeStep(rs, ctrl);
        return OCP_FALSE;
    }

    NR.UpdateIter(abs(status));

#ifdef DEBUG
    //OCP_INT global_rank = rs.domain.global_rank;
    //ls.OutputLinearSystem("proc" + to_string(CURRENT_RANK) + "_testA_AIMc.out",
    //                      "proc" + to_string(CURRENT_RANK) + "_testb_AIMc.out");
    //MPI_Barrier(rs.domain.global_comm);
    //ls.OutputSolution("proc" + to_string(CURRENT_RANK) + "_testx_AIMc.out");
#endif // DEBUG

    timer.Start();
    GetSolution(rs, ls.GetSolution(), ctrl.NR);
    OCPTIME_NRSTEP += timer.Stop();
    ls.ClearData();

    return OCP_TRUE;
}

OCP_BOOL IsoT_AIMc::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    // First check: Ni check and bulk Pressure check
    if (!NR.CheckPhysical(rs, { "BulkNi", "BulkP" }, ctrl.time.GetCurrentDt())) {
        ctrl.time.CutDt(NR);
        ResetToLastTimeStep(rs, ctrl);
        return OCP_FALSE;
    }

    CalFlashI(rs.bulk);
    CalFlashEp(rs.bulk);
    CalKrPcI(rs.bulk);

    CalRock(rs.bulk);

    rs.allWells.CalFlux(rs.bulk);

    CalRes(rs, ctrl.time.GetCurrentDt());
    return OCP_TRUE;
}

OCP_BOOL IsoT_AIMc::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    NR.CalMaxChangeNR(rs);
    const OCPNRStateC conflag = ctrl.CheckConverge(NR, { "res", "d" });

    if (conflag == OCPNRStateC::converge) {
        if (!NR.CheckPhysical(rs, { "WellP" }, ctrl.time.GetCurrentDt())) {
            ctrl.time.CutDt(NR);
            ResetToLastTimeStep(rs, ctrl);
            return OCP_FALSE;
        } else {
            CalFlashEa(rs.bulk);
            CalKrPcE(rs.bulk);
            return OCP_TRUE;
        }

    } else if (conflag == OCPNRStateC::not_converge) {
        ctrl.time.CutDt();
        ResetToLastTimeStep(rs, ctrl);
        return OCP_FALSE;
    } else {
        return OCP_FALSE;
    }
}

/// Finish a time step.
void IsoT_AIMc::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.CalIPRT(ctrl.time.GetCurrentDt());
    NR.CalMaxChangeTime(rs);
    ctrl.CalNextTimeStep(NR, {"dP", "dS", "iter"});
    UpdateLastTimeStep(rs);
}

/// Allocate memory for reservoir
void IsoT_AIMc::AllocateReservoir(Reservoir& rs)
{
    IsoT_FIM::AllocateReservoir(rs);

    Bulk&         bk  = rs.bulk;
    BulkVarSet&   bvs = bk.vs;
    const OCP_USI nb  = bvs.nb;
    const USI     np  = bvs.np;

    bvs.vj.resize(nb * np);
    bvs.lvj.resize(nb * np);

    bk.bulkTypeAIM.Setup(nb);
}

void IsoT_AIMc::SetFIMBulk(Reservoir& rs)
{
    // IMPORTANT: implicity of the same grid in different processes should be consistent
    
    const OCP_INT nlayers = 2;

    // We just consider at most 1 layer neighbor now

    Bulk&           bk   = rs.bulk;
    BulkVarSet&     bvs  = bk.vs;
    const BulkConn& conn = rs.conn;
    const OCP_USI   nb   = bvs.nbI;
    const USI       np   = bvs.np;
    const USI       nc   = bvs.nc;

    // all impec
    bk.bulkTypeAIM.Init(-1);

    OCP_USI  bIdc;
    OCP_BOOL flag;

    for (OCP_USI n = 0; n < nb; n++) {
        bIdc = n * nc;
        flag = OCP_FALSE;
        // CFL
        for (USI j = 0; j < np; j++) {
            if (NR.GetCFL(n,j) > 0.8) {
                flag = OCP_TRUE;
                break;
            }
        }
        // Volume error
        if (!flag) {
            if ((fabs(bvs.vf[n] - bvs.rockVp[n]) / bvs.rockVp[n]) > 1E-3) {
                flag = OCP_TRUE;
            }
        }

        // NR Step
        if (!flag && OCP_FALSE) {
            // dP
            if (fabs(NR.DP(n) / bvs.P[n]) > 1E-3) {
                flag = OCP_TRUE;
            }
            // dNi
            if (!flag) {
                for (USI i = 0; i < bvs.nc; i++) {
                    if (fabs(NR.DN(n,i) / bvs.Ni[bIdc + i]) > 1E-3) {
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


    USI iter = 0;
    vector<vector<OCP_INT>> recv_buffer(domain.recv_element_loc.size());
    for (const auto& r : domain.recv_element_loc) {
        const auto& rv = r.second;
        auto&       rb = recv_buffer[iter];
        rb.resize(rv[1] - rv[0]);
        MPI_Irecv(&rb[0], rb.size(), OCPMPI_INT, r.first, 0, domain.global_comm, &domain.recv_request[iter]);
        iter++;
    }

    iter = 0;
    vector<vector<OCP_INT>> send_buffer(domain.send_element_loc.size());
    for (const auto& s : domain.send_element_loc) {
        const auto& sv = s.second;
        auto&       sb = send_buffer[iter];
        sb.reserve(sv.size());
        for (const auto& sv1 : sv) {
            sb.push_back(bk.bulkTypeAIM.GetBulkType(sv1));
        }
        MPI_Isend(sb.data(), sb.size(), OCPMPI_INT, s.first, 0, domain.global_comm, &domain.send_request[iter]);
        iter++;
    }


    MPI_Waitall(iter, domain.recv_request.data(), MPI_STATUS_IGNORE);

    iter = 0;
    for (const auto& r : domain.recv_element_loc) {
        const auto& rv = r.second;
        auto&       rb = recv_buffer[iter];

        for (OCP_USI n = 0; n < rb.size(); n++) {
            SetKNeighbor(conn.neighbor, n + rv[0], bk.bulkTypeAIM, rb[n]);
        }

        iter++;
    }

    MPI_Waitall(iter, domain.send_request.data(), MPI_STATUS_IGNORE);


    // Check Consistency
    iter = 0;
    for (const auto& r : domain.recv_element_loc) {
        const auto& rv = r.second;
        auto&       rb = recv_buffer[iter];
        rb.resize(rv[1] - rv[0]);
        MPI_Irecv(&rb[0], rb.size(), OCPMPI_INT, r.first, 0, domain.global_comm, &domain.recv_request[iter]);
        iter++;
    }

    iter = 0;
    for (const auto& s : domain.send_element_loc) {
        const auto& sv = s.second;
        auto&       sb = send_buffer[iter];
        sb.clear();
        sb.reserve(sv.size());
        for (const auto& sv1 : sv) {
            sb.push_back(bk.bulkTypeAIM.GetBulkType(sv1));
        }
        MPI_Isend(sb.data(), sb.size(), OCPMPI_INT, s.first, 0, domain.global_comm, &domain.send_request[iter]);
        iter++;
    }

    MPI_Waitall(iter, domain.recv_request.data(), MPI_STATUS_IGNORE);
     
    iter = 0;
    for (const auto& r : domain.recv_element_loc) {
        const auto& rv = r.second;
        auto&       rb = recv_buffer[iter];

        for (OCP_USI n = 0; n < rb.size(); n++) {
            SetKNeighbor(conn.neighbor, n + rv[0], bk.bulkTypeAIM, rb[n]);  // Maybe not a good idea
        }
        iter++;
    }

    MPI_Waitall(iter, domain.send_request.data(), MPI_STATUS_IGNORE);

    if (OCP_TRUE) {
        cout << fixed << setprecision(2) << "Rank " << CURRENT_RANK << "  " << bk.bulkTypeAIM.GetNumFIMBulk() * 1.0 / bvs.nb * 100 << "% " << endl;
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
    const BulkVarSet& bvs = bk.vs;
    const OCP_USI&    nb  = bvs.nb;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bulkTypeAIM.IfIMPECbulk(n)) {
            // Explicit bulk

            bk.PVTm.GetPVT(n)->FlashIMPEC(n, bvs);
            // bvs.PassFlashValueAIMcEp(n);
            PassFlashValueEp(bk, n);
        }
    }
}

void IsoT_AIMc::CalFlashEa(Bulk& bk)
{
    const BulkVarSet& bvs = bk.vs;
    const OCP_USI&    nb  = bvs.nb;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bulkTypeAIM.IfIMPECbulk(n)) {
            // Explicit bulk

            bk.PVTm.GetPVT(n)->FlashIMPEC(n, bvs);
            // bvs.PassFlashValueAIMcEa(n);

            IsoT_IMPEC::PassFlashValue(bk, n);
        }
    }
}

void IsoT_AIMc::CalFlashI(Bulk& bk)
{
    BulkVarSet&    bvs = bk.vs;
    const OCP_USI& nb  = bvs.nb;
    const USI&     np  = bvs.np;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bulkTypeAIM.IfFIMbulk(n)) {
            // Implicit bulk

            bk.PVTm.GetPVT(n)->FlashFIM(n, bvs);
            IsoT_FIM::PassFlashValue(bk, n);
            for (USI j = 0; j < np; j++) {
                bvs.vj[n * np + j] = bvs.vf[n] * bvs.S[n * np + j];
            }
        }
    }
}

void IsoT_AIMc::PassFlashValueEp(Bulk& bk, const OCP_USI& n)
{
    // only var about volume needs, some flash var also
    OCP_FUNCNAME;
    auto&         bvs  = bk.vs;
    const auto    PVT  = bk.PVTm.GetPVT(n);

    const auto    np   = bvs.np;
    const auto    nc   = bvs.nc;
    const OCP_USI bIdp = n * np;

    bvs.Nt[n]  = PVT->GetNt();
    bvs.vf[n]  = PVT->GetVf();
    bvs.vfP[n] = PVT->GetVfP();
    for (USI i = 0; i < nc; i++) {
        bvs.vfi[n * nc + i] = PVT->GetVfi(i);
    }

    for (USI j = 0; j < np; j++) {
        if (PVT->GetPhaseExist(j)) {

            // IMPORTANT -- need for next Flash
            // But xij in nonlinear equations has been modified
            for (USI i = 0; i < nc; i++) {
                bvs.xij[bIdp * nc + j * nc + i] = PVT->GetXij(j, i);
            }
        }
    }
}

void IsoT_AIMc::CalKrPcE(Bulk& bk)
{
    BulkVarSet&    bvs = bk.vs;
    const OCP_USI& nb = bvs.nb;
    const USI&     np = bvs.np;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bulkTypeAIM.IfIMPECbulk(n)) {

            auto SAT = bk.SATm.GetSAT(n);
            // Explicit bulk
            const OCP_USI bId = n * np;
            SAT->CalKrPc(n, &bvs.S[bId]);
            copy(SAT->GetKr().begin(), SAT->GetKr().end(), &bvs.kr[bId]);
            copy(SAT->GetPc().begin(), SAT->GetPc().end(), &bvs.Pc[bId]);
            for (USI j = 0; j < np; j++) bvs.Pj[bId + j] = bvs.P[n] + bvs.Pc[bId + j];
        }
    }
}

void IsoT_AIMc::CalKrPcI(Bulk& bk)
{
    BulkVarSet&    bvs = bk.vs;
    const OCP_USI& nb = bvs.nb;
    const USI&     np = bvs.np;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bulkTypeAIM.IfFIMbulk(n)) {
            auto SAT = bk.SATm.GetSAT(n);
            // Implicit bulk
            const OCP_USI bId = n * np;
            SAT->CalKrPcFIM(n, &bvs.S[bId]);
            copy(SAT->GetKr().begin(), SAT->GetKr().end(), &bvs.kr[bId]);
            copy(SAT->GetPc().begin(), SAT->GetPc().end(), &bvs.Pc[bId]);
            copy(SAT->GetdKrdS().begin(), SAT->GetdKrdS().end(), &bvs.dKrdS[bId * np]);
            copy(SAT->GetdPcdS().begin(), SAT->GetdPcdS().end(), &bvs.dPcdS[bId * np]);
            for (USI j = 0; j < np; j++) bvs.Pj[bId + j] = bvs.P[n] + bvs.Pc[bId + j];
        }
    }
}

void IsoT_AIMc::AssembleMatBulks(LinearSystem&    ls,
                                 const Reservoir& rs,
                                 const OCP_DBL&   dt) const
{
    const USI numWell = rs.GetNumOpenWell();

    const Bulk&       bk      = rs.bulk;
    const BulkVarSet& bvs     = bk.vs;
    const BulkConn&   conn    = rs.conn;
    const OCP_USI     nbI     = bvs.nbI;
    const USI         np      = bvs.np;
    const USI         nc      = bvs.nc;
    const USI         ncol    = nc + 1;
    const USI         ncol2   = np * nc + np;
    const USI         bsize   = ncol * ncol;
    const USI         bsize2  = ncol * ncol2;

    ls.AddDim(nbI);

    vector<OCP_DBL> bmat(bsize, 0);
    // Accumulation term
    for (USI i = 1; i < nc + 1; i++) {
        bmat[i * ncol + i] = 1;
    }
    for (OCP_USI n = 0; n < nbI; n++) {
        bmat[0] = bvs.v[n] * bvs.poroP[n] - bvs.vfP[n];
        for (USI i = 0; i < nc; i++) {
            bmat[i + 1] = -bvs.vfi[n * nc + i];
        }
        ls.NewDiag(n, bmat);
    }

    // flux term
    OCP_BOOL          bIdFIM, eIdFIM;
    OCP_USI           bId, eId;
      
    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId       = conn.iteratorConn[c].BId();
        eId       = conn.iteratorConn[c].EId();
        auto Flux = conn.FLUXm.GetFlux(c);

        if (bk.bulkTypeAIM.IfFIMbulk(bId))  bIdFIM = OCP_TRUE;
        else                                bIdFIM = OCP_FALSE;

        if (bk.bulkTypeAIM.IfFIMbulk(eId))  eIdFIM = OCP_TRUE;
        else                                eIdFIM = OCP_FALSE;

        Flux->AssembleMatAIM(conn.iteratorConn[c], c, conn.vs, bk);

        // Assemble
        bmat = Flux->GetdFdXpB();
        if (bIdFIM) {
            DaABpbC(ncol, ncol, ncol2, 1, Flux->GetdFdXsB().data(), &bvs.dSec_dPri[bId * bsize2],
                    1, bmat.data());
        }
        Dscalar(bsize, dt, bmat.data());
        // Begin - Begin -- add
        ls.AddDiag(bId, bmat);
        // End - Begin -- insert
        if (eId < nbI) {
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
        bmat = Flux->GetdFdXpE();
        if (eIdFIM) {
            DaABpbC(ncol, ncol, ncol2, 1, Flux->GetdFdXsE().data(), &bvs.dSec_dPri[eId * bsize2],
                    1, bmat.data());
        }
        Dscalar(bsize, dt, bmat.data());

        if (eId < nbI) {
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

void IsoT_AIMc::GetSolution(Reservoir&       rs,
                            vector<OCP_DBL>& u,
                            const ControlNR& ctrlNR)
{
    const Domain&   domain = rs.domain;
    Bulk&           bk     = rs.bulk;
    BulkVarSet&     bvs    = bk.vs;
    const OCP_USI   nb     = bvs.nb;
    const USI       np     = bvs.np;
    const USI       nc     = bvs.nc;
    const USI       row    = np * (nc + 1);
    const USI       col    = nc + 1;

    // Well first
    USI wId = bvs.nbI * col;
    for (auto& wl : rs.allWells.wells) {
        wl->GetSolutionFIM(u, wId);
    }

    // Exchange Solution for ghost grid
    USI iter = 0;
    for (const auto& r : domain.recv_element_loc) {
        const auto& rv = r.second;
        MPI_Irecv(&u[rv[0] * col], (rv[1] - rv[0]) * col, OCPMPI_DBL, r.first, 0, domain.global_comm, &domain.recv_request[iter]);
        iter++;
    }

    iter = 0;
    vector<vector<OCP_DBL>> send_buffer(domain.send_element_loc.size());
    for (const auto& s : domain.send_element_loc) {
        const auto& sv = s.second;
        auto&       sb = send_buffer[iter];
        sb.reserve(sv.size() * col);
        for (const auto& sv1 : sv) {
            const OCP_DBL* bId = u.data() + sv1 * col;
            sb.insert(sb.end(), bId, bId + col);
        }
        MPI_Isend(sb.data(), sb.size(), OCPMPI_DBL, s.first, 0, domain.global_comm, &domain.send_request[iter]);
        iter++;
    }


    // Bulk
    const OCP_DBL dSmaxlim = ctrlNR.DSmax();
    // const OCP_DBL dPmaxlim = ctrlNR.dPmax;

    vector<OCP_DBL> dtmp(row, 0);
    OCP_DBL         chopmin = 1;
    OCP_DBL         choptmp = 0;

    OCP_USI bId = 0;
    OCP_USI eId = bvs.nbI;

    for (USI p = 0; p < 2; p++) {

        for (OCP_USI n = bId; n < eId; n++) {
            if (bk.bulkTypeAIM.IfIMPECbulk(n)) {
                // IMPEC Bulk
                // Pressure
                bvs.P[n] += u[n * col]; // seems better
                // Ni
                for (USI i = 0; i < nc; i++) {
                    bvs.Ni[n * nc + i] += u[n * col + 1 + i];

                    // if (bvs.Ni[n * nc + i] < 0 && bvs.Ni[n * nc + i] > -1E-3) {
                    //     bvs.Ni[n * nc + i] = 1E-20;
                    // }
                }
                // Pj
                for (USI j = 0; j < np; j++) {
                    bvs.Pj[n * np + j] = bvs.P[n] + bvs.Pc[n * np + j];
                }
                continue;
            }

            chopmin = 1;
            // compute the chop
            fill(dtmp.begin(), dtmp.end(), 0.0);
            OCP_aAxpby(row, col, static_cast<OCP_DBL>(1.0), &bvs.dSec_dPri[n * bvs.lendSdP], u.data() + n * col, static_cast<OCP_DBL>(1.0), dtmp.data());

            for (USI j = 0; j < np; j++) {
                choptmp = 1;
                if (fabs(dtmp[j]) > dSmaxlim) {
                    choptmp = dSmaxlim / fabs(dtmp[j]);
                }
                else if (bvs.S[n * np + j] + dtmp[j] < 0.0) {
                    choptmp = 0.9 * bvs.S[n * np + j] / fabs(dtmp[j]);
                }

                // if (fabs(S[n * np + j] - scm[j]) > TINY &&
                //     (S[n * np + j] - scm[j]) / (choptmp * dtmp[js]) < 0)
                //     choptmp *= min(1.0, -((S[n * np + j] - scm[j]) / (choptmp * dtmp[js])));

                chopmin = min(chopmin, choptmp);
            }

            // dS
            for (USI j = 0; j < np; j++) {
                bvs.S[n * np + j] += chopmin * dtmp[j];
            }

            // dxij
			USI js = np;
			for (USI j = 0; j < np; j++) {
				for (USI i = 0; i < bvs.nc; i++) {
					bvs.xij[(n * np + j) * nc + i] += chopmin * dtmp[js];
					js++;
				}
			}

            // dP
            bvs.P[n] += u[n * col]; // seems better

            // dNi
            for (USI i = 0; i < nc; i++) {
                bvs.Ni[n * nc + i] += chopmin * u[n * col + 1 + i];

                // if (bvs.Ni[n * nc + i] < 0 && bvs.Ni[n * nc + i] > -1E-3) {
                //     bvs.Ni[n * nc + i] = 1E-20;
                // }
            }
        }
        if (p == 0) {
            bId = eId;
            eId = nb;
            MPI_Waitall(iter, domain.recv_request.data(), MPI_STATUS_IGNORE);
        }
        else {
            break;
        }
    }

    MPI_Waitall(iter, domain.send_request.data(), MPI_STATUS_IGNORE);
}

void IsoT_AIMc::ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.bulk.vs.vj    = rs.bulk.vs.lvj;
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

    rs.bulk.vs.lvj   = rs.bulk.vs.vj;
}


////////////////////////////////////////////
// IsoT_FIMddm
////////////////////////////////////////////


void IsoT_FIMddm::Setup(Reservoir& rs, const OCPControl& ctrl)
{
    // Allocate memory for reservoir
    AllocateReservoir(rs);
    rs.domain.SetNumNprocNproc();
}


void IsoT_FIMddm::Prepare(Reservoir& rs, const OCP_DBL& dt)
{
    // Calculate well property at the beginning of next time step
    rs.allWells.PrepareWell(rs.bulk);
    // Calculate initial residual
    CalRes(rs, dt, OCP_TRUE);
    NR.InitStep(rs.bulk.GetVarSet());
    NR.InitIter();
    rs.domain.SetLSComm(starBulkSet);
    CalRankSet(rs.domain);
}


OCP_BOOL IsoT_FIMddm::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    if (!NR.CheckPhysical(rs, { "BulkNi", "BulkP" }, ctrl.time.GetCurrentDt())) {
        ctrl.time.CutDt(NR);
        ResetToLastTimeStep(rs, ctrl);
        return OCP_FALSE;
    }

    // Update fluid property
    CalFlash(rs.bulk, rankSetInLS, rs.domain);
    CalKrPc(rs.bulk, rankSetInLS, rs.domain);
    // Update rock property
    CalRock(rs.bulk, rankSetInLS, rs.domain);
    // Update well property
    rs.allWells.CalFlux(rs.bulk);

    // Update residual
    CalRes(rs, ctrl.time.GetCurrentDt(), OCP_FALSE);

    return OCP_TRUE;
}


OCP_BOOL IsoT_FIMddm::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    if (preM) {
        // check residual for each local nonlinear equations
        NR.CalMaxChangeNR(rs);
        const OCPNRStateC conflag = ctrl.CheckConverge(NR, { "res", "d" });

        if (conflag == OCPNRStateC::converge) {
            if (!NR.CheckPhysical(rs, { "WellP" }, ctrl.time.GetCurrentDt())) {
                ctrl.time.CutDt(NR);
                ResetToLastTimeStep(rs, ctrl);
                return OCP_FALSE;
            }
            else {
                // exchange solution
                ExchangePBoundary(rs);
                ExchangeNiBoundary(rs);
                UpdatePropertyBoundary(rs, ctrl);
                return OCP_TRUE;
            }
        }
        else if (conflag == OCPNRStateC::not_converge) {
            ctrl.time.CutDt();
            ResetToLastTimeStep(rs, ctrl);
            return OCP_FALSE;
        }
        else {
            return OCP_FALSE;
        }
    }
    else {
        // check residual for global nonlinear equations
        // exchange solution
        ExchangePBoundary(rs);
        ExchangeNiBoundary(rs);
        UpdatePropertyBoundary(rs, ctrl);

        NR.res.maxRelRes0_V = global_res0;
        IsoT_FIM::CalRes(rs, ctrl.time.GetCurrentDt());
        
        // cout << scientific << setprecision(3);
        // cout << CURRENT_RANK << "   " << global_res0 << "   " << NR.res.maxRelRes_V << endl;

        NR.CalMaxChangeNR(rs);
        const OCPNRStateC conflag = ctrl.CheckConverge(NR, { "res", "d" });
        if (conflag == OCPNRStateC::converge) {
            if (!NR.CheckPhysical(rs, { "WellP" }, ctrl.time.GetCurrentDt())) {
                ctrl.time.CutDt(NR);
                ResetToLastTimeStep(rs, ctrl);
                return OCP_FALSE;
            }
            else {
                return OCP_TRUE;
            }
        }
        else if (conflag == OCPNRStateC::not_converge) {
            ctrl.time.CutDt();
            ResetToLastTimeStep(rs, ctrl);
            return OCP_FALSE;
        }
        else {
            return OCP_FALSE;
        }
    }
}


void IsoT_FIMddm::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.CalIPRT(ctrl.time.GetCurrentDt());
    NR.CalMaxChangeTime(rs);
    ctrl.CalNextTimeStep(NR, { "dP", "dS", "iter" });
    SetStarBulkSet(rs.bulk, rs.domain);

    UpdateLastTimeStep(rs);
}


void IsoT_FIMddm::SetStarBulkSet(const Bulk& bulk, const Domain& domain)
{
    starBulkSet.clear();

    const auto& bvs = bulk.vs;
    for (OCP_USI n = 0; n < bvs.nb; n++) {

        // dS
        for (USI j = 0; j < bvs.np; j++) {
            const OCP_USI n_np_j = n * bvs.np + j;
            if (fabs(bvs.S[n_np_j] - bvs.lS[n_np_j]) > dSlim) {
                starBulkSet.push_back(n);
                break;
            }
        }
    }
}


void IsoT_FIMddm::CalRankSet(const Domain& domain)
{
    rankSetInLS = domain.ls_group_global_rank;
    rankSetOutLS.clear();

	for (const auto& r : domain.recv_element_loc) {
		if (!rankSetInLS.count(r.first)) {
            rankSetOutLS.insert(r.first);
		}
	}

    //if (OCP_TRUE) {
    //    cout << "rank" << CURRENT_RANK << "  ";
    //    for (const auto& s : rankSetInLS) {
    //        cout << s << "   ";
    //    }
    //    cout << endl;
    //}
}


void IsoT_FIMddm::CalFlash(Bulk& bk, const set<OCP_INT>& rankSet, const Domain& domain)
{
    const BulkVarSet& bvs = bk.vs;
    OCP_USI bId, eId;

    for (const auto& p : rankSet) {

        if (p == CURRENT_RANK) {
            bId = 0;
            eId = bvs.nbI;
        }
        else {
            bId = domain.recv_element_loc.at(p)[0];
            eId = domain.recv_element_loc.at(p)[1];
        }

        for (OCP_USI n = bId; n < eId; n++) {
            bk.PVTm.GetPVT(n)->FlashFIM(n, bvs);
            PassFlashValue(bk, n);
        }
    }
}


void IsoT_FIMddm::CalKrPc(Bulk& bk, const set<OCP_INT>& rankSet, const Domain& domain)
{
    BulkVarSet& bvs = bk.vs;
    const USI&  np  = bvs.np;
    OCP_USI     bId, eId;

	for (const auto& p : rankSet) {

		if (p == CURRENT_RANK) {
			bId = 0;
			eId = bvs.nbI;
		}
		else {
			bId = domain.recv_element_loc.at(p)[0];
			eId = domain.recv_element_loc.at(p)[1];
		}

		for (OCP_USI n = bId; n < eId; n++) {
			auto SAT = bk.SATm.GetSAT(n);

			const OCP_USI n_np = n * np;
			SAT->CalKrPcFIM(n, &bvs.S[n_np]);
			copy(SAT->GetKr().begin(), SAT->GetKr().end(), &bvs.kr[n_np]);
			copy(SAT->GetPc().begin(), SAT->GetPc().end(), &bvs.Pc[n_np]);
			copy(SAT->GetdKrdS().begin(), SAT->GetdKrdS().end(), &bvs.dKrdS[n_np * np]);
			copy(SAT->GetdPcdS().begin(), SAT->GetdPcdS().end(), &bvs.dPcdS[n_np * np]);
			for (USI j = 0; j < np; j++) bvs.Pj[n_np + j] = bvs.P[n] + bvs.Pc[n_np + j];
		}
	}
}


void IsoT_FIMddm::CalRock(Bulk& bk, const set<OCP_INT>& rankSet, const Domain& domain)
{
    BulkVarSet& bvs = bk.vs;
    OCP_USI     bId, eId;

    for (const auto& p : rankSet) {

        if (p == CURRENT_RANK) {
            bId = 0;
            eId = bvs.nbI;
        }
        else {
            bId = domain.recv_element_loc.at(p)[0];
            eId = domain.recv_element_loc.at(p)[1];
        }

        for (OCP_USI n = bId; n < eId; n++) {
            auto ROCK = bk.ROCKm.GetROCK(n);

            ROCK->CalPoro(bvs.P[n], bvs.T[n], bvs.poroInit[n], BulkContent::rf);
            bvs.poro[n] = ROCK->GetPoro();
            bvs.poroP[n] = ROCK->GetdPorodP();
            bvs.rockVp[n] = bvs.v[n] * bvs.poro[n];
        }
    }
}


/// Calculate residual
void IsoT_FIMddm::CalRes(Reservoir& rs, const OCP_DBL& dt, const OCP_BOOL& initRes0)
{
    if (boundCondition == constP) {
        CalResConstP(rs, dt, initRes0);
    }
    else if (boundCondition == constV) {
        CalResConstV(rs, dt, initRes0);
    }
    else {
        OCP_ABORT("Not Used!");
    }
}


/// Use Dirichlet boundary with fixed pressure at last time setp
void IsoT_FIMddm::CalResConstP(Reservoir& rs, const OCP_DBL& dt, const OCP_BOOL& initRes0)

{
    const Bulk&       bk = rs.bulk;
    const BulkVarSet& bvs = bk.vs;

    const USI nb  = bvs.nbI;
    const USI np  = bvs.np;
    const USI nc  = bvs.nc;
    const USI len = nc + 1;

    OCPNRresidual& res = NR.res;

    res.SetZero();

    // Accumalation Term
    for (OCP_USI n = 0; n < nb; n++) {
        const vector<OCP_DBL>& r = bk.ACCm.GetAccumuTerm()->CalResFIM(n, bvs, dt);
        copy(r.begin(), r.end(), &res.resAbs[n * len]);
    }

    // Flux Term
    OCP_USI         bId, eId;
    BulkConn& conn = rs.conn;
    BulkConnVarSet& bcvs = conn.vs;
    for (OCP_USI c = 0; c < conn.numConn; c++) {

        bId = conn.iteratorConn[c].BId();
        eId = conn.iteratorConn[c].EId();
        auto Flux = conn.FLUXm.GetFlux(c);

        Flux->CalFlux(conn.iteratorConn[c], bk);
        copy(Flux->GetConvectUpblock().begin(), Flux->GetConvectUpblock().end(), &bcvs.upblock[c * np]);
        copy(Flux->GetConvectDP().begin(), Flux->GetConvectDP().end(), &bcvs.dP[c * np]);
        copy(Flux->GetConvectVj().begin(), Flux->GetConvectVj().end(), &bcvs.flux_vj[c * np]);

        if (eId < nb) {
            for (USI i = 0; i < nc; i++) {
                res.resAbs[bId * len + 1 + i] += dt * Flux->GetFluxNi()[i];
                res.resAbs[eId * len + 1 + i] -= dt * Flux->GetFluxNi()[i];
            }
        }
        else {
            for (USI i = 0; i < nc; i++) {
                res.resAbs[bId * len + 1 + i] += dt * Flux->GetFluxNi()[i];
            }
        }
    }

    // Well to Bulk, Well
    USI wId = nb * len;
    for (const auto& wl : rs.allWells.wells) {
        wl->CalResFIM(wId, res, bk, dt);
    }

    // Calculate RelRes
    OCP_DBL tmp;
    for (OCP_USI n = 0; n < nb; n++) {

        for (USI i = 0; i < len; i++) {
            tmp = fabs(res.resAbs[n * len + i] / bvs.rockVp[n]);
            if (res.maxRelRes_V < tmp) {
                res.maxRelRes_V = tmp;
                res.maxId_V = n;
            }
            res.resRelV[n] += tmp * tmp;
        }
        res.resRelV[n] = sqrt(res.resRelV[n]);

        for (USI i = 1; i < len; i++) {
            tmp = fabs(res.resAbs[n * len + i] / bvs.Nt[n]);
            if (res.maxRelRes_N < tmp) {
                res.maxRelRes_N = tmp;
                res.maxId_N = n;
            }
            res.resRelN[n] += tmp * tmp;
        }
        res.resRelN[n] = sqrt(res.resRelN[n]);
    }

    Dscalar(res.resAbs.size(), -1.0, res.resAbs.data());

    if (initRes0) {
        res.maxRelRes0_V = res.maxRelRes_V;

        GetWallTime timer;
        timer.Start();
        MPI_Allreduce(&res.maxRelRes_V, &global_res0, 1, OCPMPI_DBL, MPI_MIN, rs.domain.global_comm);
        OCPTIME_COMM_COLLECTIVE += timer.Stop();
    }
}

/// Use Dirichlet boundary with fixed flow rate at last time setp
void IsoT_FIMddm::CalResConstV(Reservoir& rs, const OCP_DBL& dt, const OCP_BOOL& initRes0)
{

}


/// Assemble linear system for bulks
void IsoT_FIMddm::AssembleMatBulks(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const
{
    if (boundCondition == constP) {
        AssembleMatBulksConstP(ls ,rs, dt);
    }
    else if (boundCondition == constV) {
        AssembleMatBulksConstV(ls, rs, dt);
    }
    else {
        OCP_ABORT("Not Used!");
    }
}


/// Use Dirichlet boundary with fixed pressure at last time setp
void IsoT_FIMddm::AssembleMatBulksConstP(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const
{
    const Bulk&       bk  = rs.bulk;
    const BulkVarSet& bvs = bk.vs;

    const USI numWell = rs.GetNumOpenWell();

    const BulkConn& conn   = rs.conn;
    const OCP_USI   nbI    = bvs.nbI;
    const USI       np     = bvs.np;
    const USI       nc     = bvs.nc;
    const USI       ncol   = nc + 1;
    const USI       ncol2  = np * nc + np;
    const USI       bsize  = ncol * ncol;
    const USI       bsize2 = ncol * ncol2;

    ls.AddDim(nbI);


    // Accumulation term
    vector<OCP_DBL> bmat(bsize, 0);
    for (OCP_USI n = 0; n < nbI; n++) {
        ls.NewDiag(n, bk.ACCm.GetAccumuTerm()->CaldFdXpFIM(n, bvs, dt));
    }

    // flux term
    OCP_USI  bId, eId;
    for (OCP_USI c = 0; c < conn.numConn; c++) {

        bId = conn.iteratorConn[c].BId();
        eId = conn.iteratorConn[c].EId();
        auto Flux = conn.FLUXm.GetFlux(c);

        Flux->AssembleMatFIM(conn.iteratorConn[c], c, conn.vs, bk);

        bmat = Flux->GetdFdXpB();
        DaABpbC(ncol, ncol, ncol2, 1, Flux->GetdFdXsB().data(), &bvs.dSec_dPri[bId * bsize2], 1,
            bmat.data());
        Dscalar(bsize, dt, bmat.data());

        // Assemble
        // Begin - Begin -- add
        ls.AddDiag(bId, bmat);
        // End - Begin -- insert
        if (eId < nbI) {
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
        if (eId < nbI) {
            // process Interior grid
            bmat = Flux->GetdFdXpE();
            DaABpbC(ncol, ncol, ncol2, 1, Flux->GetdFdXsE().data(), &bvs.dSec_dPri[eId * bsize2], 1,
                bmat.data());
            Dscalar(bsize, dt, bmat.data());
            // Begin - End -- insert
            ls.NewOffDiag(bId, eId, bmat);
            // End - End -- add
            Dscalar(bsize, -1, bmat.data());
            ls.AddDiag(eId, bmat);
        }
        else if (IfBulkInLS(eId, rs.domain)) {
            // group Interior grid
            bmat = Flux->GetdFdXpE();
            DaABpbC(ncol, ncol, ncol2, 1, Flux->GetdFdXsE().data(), &bvs.dSec_dPri[eId * bsize2], 1,
                bmat.data());
            Dscalar(bsize, dt, bmat.data());
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



/// Use Dirichlet boundary with fixed flow rate at last time setp
void IsoT_FIMddm::AssembleMatBulksConstV(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const
{

}


/// Only get local solution
void IsoT_FIMddm::GetSolution(Reservoir& rs, vector<OCP_DBL>& u, const ControlNR& ctrlNR)
{
    const auto& domain = rs.domain;
    auto&       bk     = rs.bulk;
    auto&       bvs    = bk.vs;
    const auto  nb     = bvs.nb;
    const auto  np     = bvs.np;
    const auto  nc     = bvs.nc;
    const auto  row    = np * (nc + 1);
    const auto  col    = nc + 1;

    // Well first
    USI wId = bvs.nbI * col;
    for (auto& wl : rs.allWells.wells) {
        wl->GetSolutionFIM(u, wId);
    }

    GetWallTime timerT;         ///< total timer
    GetWallTime timerC;         ///< calculation timer
    OCP_DBL     time_cal = 0;   ///< calculation time
    timerT.Start();

    
    USI iter = 0;
    for (const auto& r : domain.recv_element_loc) {

        if (!domain.IfIRankInLSCommGroup(r.first))  continue; 

        const auto& rv = r.second;
        MPI_Irecv(&u[rv[0] * col], (rv[1] - rv[0]) * col, OCPMPI_DBL, r.first, 0, domain.global_comm, &domain.recv_request[iter]);
        iter++;
    }

    iter = 0;
    vector<vector<OCP_DBL>> send_buffer(domain.send_element_loc.size());
    for (const auto& s : domain.send_element_loc) {

        if (!domain.IfIRankInLSCommGroup(s.first))  continue; 

        const auto& sv = s.second;
        auto&       sb = send_buffer[iter];
        sb.reserve(sv.size() * col);
        for (const auto& sv1 : sv) {
            const OCP_DBL* bId = u.data() + sv1 * col;
            sb.insert(sb.end(), bId, bId + col);
        }
        MPI_Isend(sb.data(), sb.size(), OCPMPI_DBL, s.first, 0, domain.global_comm, &domain.send_request[iter]);
        iter++;
    }

    // Bulk
    const OCP_DBL dSmaxlim = ctrlNR.DSmax();
    const OCP_DBL dPmaxlim = ctrlNR.DPmax();

    vector<OCP_DBL> dtmp(row, 0);
    OCP_DBL         chopmin = 1;
    OCP_DBL         choptmp = 0;

    OCP_USI bId, eId;
   
    iter = 0;
    for (const auto& p : rankSetInLS) {

        if (p == CURRENT_RANK) {
            bId = 0;
            eId = bvs.nbI;
        }
        else {
            bId = domain.recv_element_loc.at(p)[0];
            eId = domain.recv_element_loc.at(p)[1];

            MPI_Wait(&domain.recv_request[iter], MPI_STATUS_IGNORE);
            iter++;
        }

        timerC.Start();

        for (OCP_USI n = bId; n < eId; n++) {
            // const vector<OCP_DBL>& scm = satcm[SATNUM[n]];

            chopmin = 1;
            // compute the chop
            fill(dtmp.begin(), dtmp.end(), 0.0);

            OCP_aAxpby(row, col, static_cast<OCP_DBL>(1.0), &bvs.dSec_dPri[n * bvs.lendSdP], u.data() + n * col, static_cast<OCP_DBL>(1.0), dtmp.data());

            for (USI j = 0; j < np; j++) {
                choptmp = 1;
                if (fabs(dtmp[j]) > dSmaxlim) {
                    choptmp = dSmaxlim / fabs(dtmp[j]);
                }
                else if (bvs.S[n * np + j] + dtmp[j] < 0.0) {
                    choptmp = 0.9 * bvs.S[n * np + j] / fabs(dtmp[j]);
                }
                // if (fabs(S[n_np_j] - scm[j]) > TINY &&
                //     (S[n_np_j] - scm[j]) / (choptmp * dtmp[js]) < 0)
                //     choptmp *= min(1.0, -((S[n_np_j] - scm[j]) / (choptmp * dtmp[js])));
                chopmin = min(chopmin, choptmp);
            }

            // dS
            for (USI j = 0; j < np; j++) {
                bvs.S[n * np + j] += chopmin * dtmp[j];
            }

            // dxij
            USI js = np;
            for (USI j = 0; j < np; j++) {
                for (USI i = 0; i < bvs.nc; i++) {
                    bvs.xij[(n * np + j) * nc + i] += chopmin * dtmp[js];
                    js++;
                }
            }

            // dP
            //choptmp = dPmaxlim / fabs(u[n * col]);
            //chopmin = min(chopmin, choptmp);
            bvs.P[n] += u[n * col]; // seems better

            // dNi
            for (USI i = 0; i < nc; i++) {
                bvs.Ni[n * nc + i] += chopmin * u[n * col + 1 + i];

                // if (bvs.Ni[n * nc + i] < 0 && bvs.Ni[n * nc + i] > -1E-3) {
                //     bvs.Ni[n * nc + i] = 1E-20;
                // }
            }
        }

        time_cal += timerC.Stop();
    }
    

    MPI_Waitall(iter, domain.send_request.data(), MPI_STATUS_IGNORE);

    OCPTIME_COMM_P2P += (timerT.Stop() - time_cal);
    OCPTIME_NRSTEPC  += time_cal;
}


void IsoT_FIMddm::UpdatePropertyBoundary(Reservoir& rs, OCPControl& ctrl)
{
    CalFlash(rs.bulk, rankSetOutLS, rs.domain);
    CalKrPc(rs.bulk, rankSetOutLS, rs.domain);
    CalRock(rs.bulk, rankSetOutLS, rs.domain);
}


void IsoT_FIMddm::ExchangePBoundary(Reservoir& rs) const
{
    // Exchange Ghost P
    const Domain& domain = rs.domain;
    BulkVarSet&   bvs    = rs.bulk.vs;

    USI iter = 0;
    for (const auto& r : domain.recv_element_loc) {

        if (domain.IfIRankInLSCommGroup(r.first))  continue;

        const auto& rv = r.second;
        MPI_Irecv(&bvs.P[rv[0]], rv[1] - rv[0], OCPMPI_DBL, r.first, 0, domain.global_comm, &domain.recv_request[iter]);
        iter++;
    }

    iter = 0;
    vector<vector<OCP_DBL>> send_buffer(domain.send_element_loc.size());
    for (const auto& s : domain.send_element_loc) {

        if (domain.IfIRankInLSCommGroup(s.first)) continue; 

        const auto& sv = s.second;
        auto&       sb = send_buffer[iter];
        sb.reserve(sv.size());
        for (const auto& sv1 : sv) {
            sb.push_back(bvs.P[sv1]);
        }
        MPI_Isend(sb.data(), sb.size(), OCPMPI_DBL, s.first, 0, domain.global_comm, &domain.send_request[iter]);
        iter++;
    }

    MPI_Waitall(iter, domain.send_request.data(), MPI_STATUS_IGNORE);
    MPI_Waitall(iter, domain.recv_request.data(), MPI_STATUS_IGNORE);
}


void IsoT_FIMddm::ExchangeNiBoundary(Reservoir& rs) const
{
    // Exchange Ghost Ni
    const Domain& domain = rs.domain;
    BulkVarSet&   bvs = rs.bulk.vs;
    const USI     nc = bvs.nc;

    USI iter = 0;
    for (const auto& r : domain.recv_element_loc) {

        if (domain.IfIRankInLSCommGroup(r.first))  continue; 

        const auto& rv = r.second;
        MPI_Irecv(&bvs.Ni[rv[0] * nc], (rv[1] - rv[0]) * nc, OCPMPI_DBL, r.first, 0, domain.global_comm, &domain.recv_request[iter]);
        iter++;
    }

    iter = 0;
    vector<vector<OCP_DBL>> send_buffer(domain.send_element_loc.size());
    for (const auto& s : domain.send_element_loc) {

        if (domain.IfIRankInLSCommGroup(s.first))  continue; 

        const auto& sv = s.second;
        auto&       sb = send_buffer[iter];
        sb.reserve(sv.size() * nc);
        for (const auto& sv1 : sv) {
            const OCP_DBL* bId = &bvs.Ni[0] + sv1 * nc;
            sb.insert(sb.end(), bId, bId + nc);
        }
        MPI_Isend(sb.data(), sb.size(), OCPMPI_DBL, s.first, 0, domain.global_comm, &domain.send_request[iter]);
        iter++;
    }

    MPI_Waitall(iter, domain.send_request.data(), MPI_STATUS_IGNORE);
    MPI_Waitall(iter, domain.recv_request.data(), MPI_STATUS_IGNORE);
}


OCP_BOOL IsoT_FIMddm::IfBulkInLS(const USI& bId, const Domain& domain) const
{
    for (const auto& p : rankSetInLS) {
        if (p == CURRENT_RANK)  continue;
        if (bId >= domain.recv_element_loc.at(p)[0] &&
            bId <  domain.recv_element_loc.at(p)[1]) {
            return OCP_TRUE;
        }
    }

    return OCP_FALSE;
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/01/2021      Create file                          */
/*  Chensong Zhang      Jan/08/2022      Update output                        */
/*----------------------------------------------------------------------------*/