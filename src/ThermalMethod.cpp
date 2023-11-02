/*! \file    ThermalMethod.cpp
 *  \brief   ThermalMethod class declaration
 *  \author  Shizhe Li
 *  \date    Nov/10/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "ThermalMethod.hpp"

void T_FIM::Setup(Reservoir& rs, LinearSystem& ls, const OCPControl& ctrl)
{
    AllocateReservoir(rs);
    AllocateLinearSystem(ls, rs, ctrl);
}

void T_FIM::InitReservoir(Reservoir& rs)
{
    rs.bulk.Initialize(rs.domain);

    InitRock(rs.bulk);
    CalRock(rs.bulk);

    InitFlash(rs.bulk);
    CalKrPc(rs.bulk);

    rs.bulk.optMs.heatConduct.CalHeatConduct(rs.bulk.vs);

    rs.allWells.InitBHP(rs.bulk);
    UpdateLastTimeStep(rs);
}

void T_FIM::Prepare(Reservoir& rs, const OCPControl& ctrl)
{
    rs.allWells.PrepareWell(rs.bulk);
    CalRes(rs, ctrl.time.GetCurrentDt());
    NR.InitStep(rs.bulk.GetVarSet());
    NR.InitIter();
}

void T_FIM::AssembleMat(LinearSystem&    ls,
                        const Reservoir& rs,
                        const OCP_DBL&   dt)
{
    AssembleMatBulks(ls, rs, dt);
    AssembleMatWells(ls, rs, dt);
    ls.AssembleRhsCopy(NR.res.resAbs);
}

void T_FIM::SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl)
{
#ifdef DEBUG
    ls.CheckEquation();
#endif // DEBUG

    GetWallTime timer;

    // Assemble external linear solver with internal A and b
    timer.Start();
    ls.CalCommTerm(rs.GetNumOpenWell());
    ls.AssembleMatLinearSolver();
    OCPTIME_ASSEMBLE_MAT_FOR_LS += timer.Stop() / TIME_S2MS;

    // Solve linear system  
    timer.Start();
    int status = ls.Solve();
    if (status < 0) {
        status = ls.GetNumIters();
    }
    // Record time, iterations
    OCPTIME_LSOLVER += timer.Stop() / TIME_S2MS;

    NR.UpdateIter(status);

#ifdef DEBUG
    // Output A, b, x
     //ls.OutputLinearSystem("testA_FIMT.out", "testb_FIMT.out");
     //ls.OutputSolution("testx_FIMT.out");
    // Check if inf or nan occurs in solution
    ls.CheckSolution();
#endif // DEBUG
    
    timer.Start();
    GetSolution(rs, ls.GetSolution(), ctrl.NR);
    OCPTIME_NRSTEP += timer.Stop() / TIME_S2MS;
    // rs.PrintSolFIM(ctrl.workDir + "testPNi.out");
    ls.ClearData();
}

OCP_BOOL T_FIM::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{

    if (!NR.CheckPhysical(rs, { "BulkNi", "BulkP", "BulkT" })) {
        ctrl.time.CutDt(NR);
        ResetToLastTimeStep(rs, ctrl);   
        return OCP_FALSE;
    }

    // Update reservoir properties
    CalRock(rs.bulk);

    CalFlash(rs.bulk);
    CalKrPc(rs.bulk);

    rs.bulk.optMs.heatConduct.CalHeatConduct(rs.bulk.vs);
    rs.bulk.optMs.boundary.heatLoss.CalHeatLoss(rs.bulk.vs, ctrl.time.GetCurrentTime() + ctrl.time.GetCurrentDt(), ctrl.time.GetCurrentDt());

    rs.allWells.CalFlux(rs.bulk);

    CalRes(rs, ctrl.time.GetCurrentDt());

    return OCP_TRUE;
}


OCP_BOOL T_FIM::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    NR.CalMaxChangeNR(rs);
    const OCPNRStateC conflag = ctrl.CheckConverge(NR, { "resT", "dT" });

    if (conflag == OCPNRStateC::converge) {
        if (!NR.CheckPhysical(rs, { "WellP" })) {
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

void T_FIM::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.CalIPRT(ctrl.time.GetCurrentDt());
    NR.CalMaxChangeTime(rs);
    ctrl.CalNextTimeStep(NR, {"dP", "dS", "iter"});
    UpdateLastTimeStep(rs);
}

void T_FIM::AllocateReservoir(Reservoir& rs)
{
    Bulk&         bk = rs.bulk;
    BulkVarSet&   bvs = bk.vs;
    const OCP_USI nb = bvs.nb;
    const USI     np = bvs.np;
    const USI     nc = bvs.nc;

    // Rock
    bvs.poro.resize(nb);
    bvs.rockVp.resize(nb);
    bvs.vr.resize(nb);
    bvs.Hr.resize(nb);

    bvs.lporo.resize(nb);
    bvs.lrockVp.resize(nb);
    bvs.lvr.resize(nb);
    bvs.lHr.resize(nb);

    // derivatives
    bvs.poroP.resize(nb);
    bvs.poroT.resize(nb);
    bvs.vrP.resize(nb);
    bvs.vrT.resize(nb);
    bvs.HrT.resize(nb);

    bvs.lporoP.resize(nb);
    bvs.lporoT.resize(nb);
    bvs.lvrP.resize(nb);
    bvs.lvrT.resize(nb);
    bvs.lHrT.resize(nb);

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
    bvs.phaseExist.resize(nb * np);
    bvs.S.resize(nb * np);
    bvs.xij.resize(nb * np * nc);
    bvs.rho.resize(nb * np);
    bvs.xi.resize(nb * np);
    bvs.mu.resize(nb * np);
    bvs.kr.resize(nb * np);
    bvs.Uf.resize(nb);
    bvs.H.resize(nb * np);

    bvs.lNt.resize(nb);
    bvs.lNi.resize(nb * nc);
    bvs.lvf.resize(nb);
    bvs.lT.resize(nb);
    bvs.lP.resize(nb);
    bvs.lPj.resize(nb * np);
    bvs.lPc.resize(nb * np);
    bvs.lphaseExist.resize(nb * np);
    bvs.lS.resize(nb * np);
    bvs.lxij.resize(nb * np * nc);
    bvs.lrho.resize(nb * np);
    bvs.lxi.resize(nb * np);
    bvs.lmu.resize(nb * np);
    bvs.lkr.resize(nb * np);
    bvs.lUf.resize(nb);
    bvs.lH.resize(nb * np);

    // derivatives
    bvs.vfP.resize(nb);
    bvs.vfT.resize(nb);
    bvs.vfi.resize(nb * nc);
    bvs.rhoP.resize(nb * np);
    bvs.rhoT.resize(nb * np);
    bvs.rhox.resize(nb * nc * np);
    bvs.xiP.resize(nb * np);
    bvs.xiT.resize(nb * np);
    bvs.xix.resize(nb * nc * np);
    bvs.muP.resize(nb * np);
    bvs.muT.resize(nb * np);
    bvs.mux.resize(nb * nc * np);
    bvs.dPcdS.resize(nb * np * np);
    bvs.dKrdS.resize(nb * np * np);
    bvs.UfP.resize(nb);
    bvs.UfT.resize(nb);
    bvs.Ufi.resize(nb * nc);
    bvs.HT.resize(nb * np);
    bvs.Hx.resize(nb * np * nc);

    bvs.lvfP.resize(nb);
    bvs.lvfT.resize(nb);
    bvs.lvfi.resize(nb * nc);
    bvs.lrhoP.resize(nb * np);
    bvs.lrhoT.resize(nb * np);
    bvs.lrhox.resize(nb * nc * np);
    bvs.lxiP.resize(nb * np);
    bvs.lxiT.resize(nb * np);
    bvs.lxix.resize(nb * nc * np);
    bvs.lmuP.resize(nb * np);
    bvs.lmuT.resize(nb * np);
    bvs.lmux.resize(nb * nc * np);
    bvs.ldPcdS.resize(nb * np * np);
    bvs.ldKrdS.resize(nb * np * np);
    bvs.UfP.resize(nb);
    bvs.UfT.resize(nb);
    bvs.Ufi.resize(nb * nc);
    bvs.HT.resize(nb * np);
    bvs.Hx.resize(nb * np * nc);

    // FIM-Specified
    bvs.lendSdP = (nc + 2) * (nc + 1) * np;
    bvs.dSec_dPri.resize(nb * bvs.lendSdP);

    bvs.ldSec_dPri.resize(nb * bvs.lendSdP);

    // BulkConn
    BulkConn&      conn    = rs.conn;
    const OCP_USI& numConn = conn.numConn;

    conn.vs.upblock.resize(numConn* np);
    conn.vs.rho.resize(numConn* np);
    conn.vs.flux_vj.resize(numConn* np);

    // Allocate Residual
    NR.Setup(OCP_TRUE, bvs, rs.allWells.numWell, rs.domain);
}

void T_FIM::AllocateLinearSystem(LinearSystem&     ls,
                                 const Reservoir&  rs,
                                 const OCPControl& ctrl)
{
    ls.SetupDomain(rs.domain);
    ls.AllocateRowMem(rs.GetComNum() + 2);
    ls.AllocateColMem();
    ls.SetupLinearSolver(OCPModel::thermal, ctrl.GetWorkDir(), ctrl.GetLsFile());
}

void T_FIM::InitRock(Bulk& bk) const
{
    BulkVarSet& bvs = bk.vs;
    for (OCP_USI n = 0; n < bvs.nb; n++) {
        if (bvs.cType[n] == BulkContent::r) {
            // non fluid bulk
            bvs.poroInit[n] = 0;
            bvs.poro[n]     = 0;
            bvs.rockVp[n]   = 0;
            bvs.vr[n]       = bvs.v[n];
        }
    }
}

void T_FIM::CalRock(Bulk& bk) const
{
    BulkVarSet&   bvs = bk.vs;
    const OCP_USI nb  = bvs.nb;

    for (OCP_USI n = 0; n < nb; n++) {
        auto ROCK = bk.ROCKm.GetROCK(n);

        ROCK->CalPoro(bvs.P[n], bvs.T[n], bvs.poroInit[n], bvs.cType[n]);
        if (bvs.cType[n] == BulkContent::rf) {
            // with fluid           
            bvs.poro[n]   = ROCK->GetPoro();
            bvs.poroP[n]  = ROCK->GetdPorodP();
            bvs.poroT[n]  = ROCK->GetdPorodT();           
            bvs.vr[n]     = bvs.v[n] * ROCK->Get_Poro();
            bvs.vrP[n]    = bvs.v[n] * ROCK->Get_dPorodP();
            bvs.vrT[n]    = bvs.v[n] * ROCK->Get_dPorodT();
            bvs.rockVp[n] = bvs.v[n] * bvs.poro[n];
        }
        bvs.Hr[n]  = ROCK->GetHr();
        bvs.HrT[n] = ROCK->GetdHrdT();
    }
}

void T_FIM::InitFlash(Bulk& bk)
{
    BulkVarSet&    bvs = bk.vs;
    const OCP_USI& nb = bvs.nb;
    const USI&     nc = bvs.nc;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bvs.cType[n] == BulkContent::rf) {
            auto PVT = bk.PVTm.GetPVT(n);

            PVT->InitFlashFIM(n, bvs);
            for (USI i = 0; i < nc; i++) {
                bvs.Ni[n * nc + i] = PVT->GetNi(i);
            }
            PassFlashValue(bk, n);
        }
    }
}

void T_FIM::CalFlash(Bulk& bk)
{
    const BulkVarSet& bvs = bk.vs;
    const OCP_USI&    nb  = bvs.nb;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bvs.cType[n] == BulkContent::rf) {
            bk.PVTm.GetPVT(n)->FlashFIM(n, bvs);
            PassFlashValue(bk, n);
        }
    }
}

void T_FIM::PassFlashValue(Bulk& bk, const OCP_USI& n)
{
    auto&         bvs  = bk.vs;
    const auto    PVT  = bk.PVTm.GetPVT(n);
    const auto&   np   = bvs.np;
    const auto&   nc   = bvs.nc;
    const OCP_USI bIdp = n * np;

    bvs.Nt[n]       = PVT->GetNt();
    bvs.vf[n]       = PVT->GetVf();
    bvs.Uf[n]       = PVT->GetUf();

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
            bvs.H[bIdp + j]   = PVT->GetH(j);

            // Derivatives
            bvs.rhoP[bIdp + j] = PVT->GetRhoP(j);
            bvs.rhoT[bIdp + j] = PVT->GetRhoT(j);
            bvs.xiP[bIdp + j]  = PVT->GetXiP(j);
            bvs.xiT[bIdp + j]  = PVT->GetXiT(j);
            bvs.muP[bIdp + j]  = PVT->GetMuP(j);
            bvs.muT[bIdp + j]  = PVT->GetMuT(j);
            bvs.HT[bIdp + j]   = PVT->GetHT(j);

            for (USI i = 0; i < nc; i++) {
                bvs.xij[bIdp * nc + j * nc + i]  = PVT->GetXij(j, i);
                bvs.rhox[bIdp * nc + j * nc + i] = PVT->GetRhoX(j, i);
                bvs.xix[bIdp * nc + j * nc + i]  = PVT->GetXiX(j, i);
                bvs.mux[bIdp * nc + j * nc + i]  = PVT->GetMuX(j, i);
                bvs.Hx[bIdp * nc + j * nc + i]   = PVT->GetHx(j, i);
            }
        }
    }
    bvs.vfP[n] = PVT->GetVfP();
    bvs.vfT[n] = PVT->GetVfT();
    bvs.UfP[n] = PVT->GetUfP();
    bvs.UfT[n] = PVT->GetUfT();

    for (USI i = 0; i < nc; i++) {
        bvs.vfi[n * nc + i] = PVT->GetVfi(i);
        bvs.Ufi[n * nc + i] = PVT->GetUfi(i);
    }

    copy(PVT->GetDXsDXp().begin(), PVT->GetDXsDXp().end(), &bvs.dSec_dPri[n * bvs.lendSdP]);

}

void T_FIM::CalKrPc(Bulk& bk) const
{
    BulkVarSet& bvs = bk.vs;
    const USI& np = bvs.np;

    for (OCP_USI n = 0; n < bvs.nb; n++) {
        if (bvs.cType[n] == BulkContent::rf) {
            auto SAT = bk.SATm.GetSAT(n);

            OCP_USI bId = n * np;
            SAT->CalKrPcFIM(n, &bvs.S[bId]);
            copy(SAT->GetKr().begin(), SAT->GetKr().end(), &bvs.kr[bId]);
            copy(SAT->GetPc().begin(), SAT->GetPc().end(), &bvs.Pc[bId]);
            copy(SAT->GetdKrdS().begin(), SAT->GetdKrdS().end(), &bvs.dKrdS[bId * np]);
            copy(SAT->GetdPcdS().begin(), SAT->GetdPcdS().end(), &bvs.dPcdS[bId * np]);
            for (USI j = 0; j < np; j++) bvs.Pj[n * np + j] = bvs.P[n] + bvs.Pc[n * np + j];
        }
    }
}


void T_FIM::ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl)
{
    // Bulk
    Bulk& bk = rs.bulk;
    BulkVarSet& bvs = bk.vs;
    // Rock
    bvs.poro   = bvs.lporo;
    bvs.rockVp = bvs.lrockVp;
    bvs.vr     = bvs.lvr;
    bvs.Hr     = bvs.lHr;
    // derivatives
    bvs.poroP = bvs.lporoP;
    bvs.poroT = bvs.lporoT;
    bvs.vrP   = bvs.lvrP;
    bvs.vrT   = bvs.lvrT;
    bvs.HrT   = bvs.lHrT;

    // Fluid
    bvs.Nt         = bvs.lNt;
    bvs.Ni         = bvs.lNi;
    bvs.vf         = bvs.lvf;
    bvs.T          = bvs.lT;
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
    bvs.Uf         = bvs.lUf;
    bvs.H          = bvs.lH;
    // derivatives
    bvs.vfP       = bvs.lvfP;
    bvs.vfT       = bvs.lvfT;
    bvs.vfi       = bvs.lvfi;
    bvs.rhoP      = bvs.lrhoP;
    bvs.rhoT      = bvs.lrhoT;
    bvs.rhox      = bvs.lrhox;
    bvs.xiP       = bvs.lxiP;
    bvs.xiT       = bvs.lxiT;
    bvs.xix       = bvs.lxix;
    bvs.muP       = bvs.lmuP;
    bvs.muT       = bvs.lmuT;
    bvs.mux       = bvs.lmux;
    bvs.dPcdS     = bvs.ldPcdS;
    bvs.dKrdS     = bvs.ldKrdS;
    bvs.UfP       = bvs.lUfP;
    bvs.UfT       = bvs.lUfT;
    bvs.Ufi       = bvs.lUfi;
    bvs.HT        = bvs.lHT;
    bvs.Hx        = bvs.lHx;
    bvs.dSec_dPri = bvs.ldSec_dPri;


    // Wells
    rs.allWells.ResetToLastTimeStep(bk);

    rs.bulk.optMs.ResetToLastTimeStep();

    // Iters
    CalRes(rs, ctrl.time.GetCurrentDt());

    NR.InitStep(rs.bulk.GetVarSet());
    NR.ResetIter();
}

void T_FIM::UpdateLastTimeStep(Reservoir& rs) const
{
    // Bulk
    Bulk& bk = rs.bulk;
    BulkVarSet& bvs = bk.vs;
    // Rock
    bvs.lporo   = bvs.poro;
    bvs.lrockVp = bvs.rockVp;
    bvs.lvr     = bvs.vr;
    bvs.lHr     = bvs.Hr;
    // derivatives
    bvs.lporoP = bvs.poroP;
    bvs.lporoT = bvs.poroT;
    bvs.lvrP   = bvs.vrP;
    bvs.lvrT   = bvs.vrT;
    bvs.lHrT   = bvs.HrT;

    // Fluid
    bvs.lNt         = bvs.Nt;
    bvs.lNi         = bvs.Ni;
    bvs.lvf         = bvs.vf;
    bvs.lT          = bvs.T;
    bvs.lP          = bvs.P;
    bvs.lPj         = bvs.Pj;
    bvs.lPc         = bvs.Pc;
    bvs.lphaseExist = bvs.phaseExist;
    bvs.lS          = bvs.S;
    bvs.lxij        = bvs.xij;
    bvs.lrho        = bvs.rho;
    bvs.lxi         = bvs.xi;
    bvs.lmu         = bvs.mu;
    bvs.lkr         = bvs.kr;
    bvs.lUf         = bvs.Uf;
    bvs.lH          = bvs.H;
    // derivatives
    bvs.lvfP       = bvs.vfP;
    bvs.lvfT       = bvs.vfT;
    bvs.lvfi       = bvs.vfi;
    bvs.lrhoP      = bvs.rhoP;
    bvs.lrhoT      = bvs.rhoT;
    bvs.lrhox      = bvs.rhox;
    bvs.lxiP       = bvs.xiP;
    bvs.lxiT       = bvs.xiT;
    bvs.lxix       = bvs.xix;
    bvs.lmuP       = bvs.muP;
    bvs.lmuT       = bvs.muT;
    bvs.lmux       = bvs.mux;
    bvs.ldPcdS   = bvs.dPcdS;
    bvs.ldKrdS    = bvs.dKrdS;
    bvs.lUfP       = bvs.UfP;
    bvs.lUfT       = bvs.UfT;
    bvs.lUfi       = bvs.Ufi;
    bvs.lHT        = bvs.HT;
    bvs.lHx        = bvs.Hx;
    bvs.ldSec_dPri = bvs.dSec_dPri;


    rs.allWells.UpdateLastTimeStep();
    rs.bulk.optMs.UpdateLastTimeStep();
}

void T_FIM::CalRes(Reservoir& rs, const OCP_DBL& dt)
{
    const Bulk& bk   = rs.bulk;
    const BulkVarSet& bvs = bk.vs;
    const USI   nb   = bvs.nbI;
    const USI   np   = bvs.np;
    const USI   nc   = bvs.nc;
    const USI   len  = nc + 2;

    OCPNRresidual& res = NR.res;
    
    res.SetZero();

    // Bulk to Bulk
    
    // Accumalation Term
    for (OCP_USI n = 0; n < nb; n++) {
        const vector<OCP_DBL>& r = bk.ACCm.GetAccumuTerm()->CalResFIM(n, bvs, dt);
        copy(r.begin(), r.end(), &res.resAbs[n * len]);
    }

    BulkConn&         conn = rs.conn;
    BulkConnVarSet&   bcvs = conn.vs;
    USI               fluxnum;
    OCP_USI           bId, eId;

    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId       = conn.iteratorConn[c].BId();
        eId       = conn.iteratorConn[c].EId();
        fluxnum   = conn.iteratorConn[c].FluxNum();
        auto Flux = conn.flux[fluxnum];

        Flux->CalFlux(conn.iteratorConn[c], bk);

        // Thermal conductive term always exists
        const auto conH = Flux->GetConductH();
        res.resAbs[bId * len + 1 + nc] += conH * dt;
        if (eId < nb) {
            // Interior grid
            res.resAbs[eId * len + 1 + nc] -= conH * dt;
        }

        if (bvs.cType[bId] == BulkContent::rf && bvs.cType[eId] == BulkContent::rf) {
            // with fluid flow
            copy(Flux->GetUpblock().begin(), Flux->GetUpblock().end(), &bcvs.upblock[c * np]);
            copy(Flux->GetRho().begin(), Flux->GetRho().end(), &bcvs.rho[c * np]);
            copy(Flux->GetFluxVj().begin(), Flux->GetFluxVj().end(), &bcvs.flux_vj[c * np]);

			if (eId < nb) {
				// Interior grid
				for (USI i = 0; i < nc; i++) {
					res.resAbs[bId * len + 1 + i] += dt * Flux->GetFluxNi()[i];
					res.resAbs[eId * len + 1 + i] -= dt * Flux->GetFluxNi()[i];
				}
                for (USI j = 0; j < np; j++) {
                    res.resAbs[bId * len + 1 + nc] += dt * Flux->GetFluxHj()[j];
                    res.resAbs[eId * len + 1 + nc] -= dt * Flux->GetFluxHj()[j];
                }
			}
			else {
				// Ghost grid
				for (USI i = 0; i < nc; i++) {
					res.resAbs[bId * len + 1 + i] += dt * Flux->GetFluxNi()[i];
				}
                for (USI j = 0; j < np; j++) {
                    res.resAbs[bId * len + 1 + nc] += dt * Flux->GetFluxHj()[j];
                }				
			}
        }
    }

    // Well to Bulk
    USI wId = nb * len;
    for (const auto& wl : rs.allWells.wells) {
        wl->CalResFIM(wId, res, bk, dt);
    }

    // Calculate RelRes
    OCP_DBL tmp;
    for (OCP_USI n = 0; n < nb; n++) {

        // Energy equations always exist
        OCP_DBL eT = bvs.vf[n] * bvs.Uf[n] + bvs.vr[n] * bvs.Hr[n];
        tmp        = fabs(res.resAbs[n * len + nc + 1] / eT);
        if (res.maxRelRes_E < tmp) {
            res.maxRelRes_E = tmp;
            res.maxId_E     = n;
        }

        if (bvs.cType[n] == BulkContent::rf) {
            // Fluid Bulk
            for (USI i = 0; i < len - 1; i++) {
                tmp = fabs(res.resAbs[n * len + i] / bvs.rockVp[n]);
                if (res.maxRelRes_V < tmp) {
                    res.maxRelRes_V = tmp;
                    res.maxId_V     = n;
                }
            }

            for (USI i = 1; i < len - 1; i++) {
                tmp = fabs(res.resAbs[n * len + i] / bvs.Nt[n]);
                if (res.maxRelRes_N < tmp) {
                    res.maxRelRes_N = tmp;
                    res.maxId_N     = n;
                }
            }
        }
    }
    Dscalar(res.resAbs.size(), -1.0, res.resAbs.data());
}

void T_FIM::AssembleMatBulks(LinearSystem&    ls,
                             const Reservoir& rs,
                             const OCP_DBL&   dt) const
{
    const USI numWell = rs.GetNumOpenWell();

    const Bulk&     bk     = rs.bulk;
    const BulkVarSet& bvs = bk.vs;
    const BulkConn& conn   = rs.conn;
    const OCP_USI   nbI    = bvs.nbI;
    const USI       np     = bvs.np;
    const USI       nc     = bvs.nc;
    const USI       ncol   = nc + 2;
    const USI       ncol2  = np * nc + np;
    const USI       bsize  = ncol * ncol;
    const USI       bsize2 = ncol * ncol2;

    ls.AddDim(nbI);

   
    // Accumulation term
    for (OCP_USI n = 0; n < nbI; n++) {
        ls.NewDiag(n, bk.ACCm.GetAccumuTerm()->CaldFdXpFIM(n, bvs, dt));
    }

    // flux term
    vector<OCP_DBL> bmat(bsize, 0);
    OCP_USI  bId, eId;
    USI      fluxnum;
    for (OCP_USI c = 0; c < conn.numConn; c++) {

        bId       = conn.iteratorConn[c].BId();
        eId       = conn.iteratorConn[c].EId();
        fluxnum   = conn.iteratorConn[c].FluxNum();
        auto Flux = conn.flux[fluxnum];

        Flux->AssembleMatFIM(conn.iteratorConn[c], c, conn.vs, bk);

        bmat = Flux->GetdFdXpB();
        DaABpbC(ncol, ncol, ncol2, 1, Flux->GetdFdXsB().data(), &bvs.dSec_dPri[bId * bsize2], 1,
            bmat.data());

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
            OCP_ABORT("INF or INF in bmat !");
        }
#endif

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
            // Ghost grid
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

void T_FIM::AssembleMatWells(LinearSystem&    ls,
                             const Reservoir& rs,
                             const OCP_DBL&   dt) const
{
    for (auto& wl : rs.allWells.wells) {
        wl->AssembleMatFIM(ls, rs.bulk, dt);
    }
}


void T_FIM::GetSolution(Reservoir&       rs,
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
    const USI       col    = nc + 2;

    // Well first
    USI wId = bvs.nbI * col;
    for (auto& wl : rs.allWells.wells) {
        wl->GetSolutionFIM(u, wId);
    }

    // Exchange Solution for ghost grid
    for (USI i = 0; i < domain.numRecvProc; i++) {
        const vector<OCP_USI>& rel = domain.recv_element_loc[i];
        MPI_Irecv(&u[rel[1] * col], (rel[2] - rel[1]) * col, OCPMPI_DBL, rel[0], 0, domain.myComm, &domain.recv_request[i]);
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
        MPI_Isend(s.data() + 1, s.size() - 1, OCPMPI_DBL, s[0], 0, domain.myComm, &domain.send_request[i]);
    }

    // Bulk
    const OCP_DBL dSmaxlim = ctrlNR.DSmax();
    // const OCP_DBL dPmaxlim = ctrlNR.dPmax;

    vector<OCP_DBL> dtmp(row, 0);
    OCP_DBL         chopmin = 1;
    OCP_DBL         choptmp = 0;

    OCP_USI bId = 0;
    OCP_USI eId = bvs.nbI;

    // interior first, ghost second
    for (USI p = bId; p < eId; p++) {

        for (OCP_USI n = bId; n < eId; n++) {
            // const vector<OCP_DBL>& scm = satcm[SATNUM[n]];

            if (bvs.cType[n] == BulkContent::rf) {
                // Fluid Bulk

                chopmin = 1;
                // compute the chop
                fill(dtmp.begin(), dtmp.end(), 0.0);
                OCP_aAxpby(np, col, static_cast<OCP_DBL>(1.0), &bvs.dSec_dPri[n * bvs.lendSdP], u.data() + n * col, static_cast<OCP_DBL>(1.0),
                    dtmp.data());

                for (USI j = 0; j < np; j++) {
                    choptmp = 1;
                    if (fabs(dtmp[j]) > dSmaxlim) {
                        choptmp = dSmaxlim / fabs(dtmp[j]);
                    }
                    else if (bvs.S[n * np + j] + dtmp[j] < 0.0) {
                        choptmp = 0.9 * bvs.S[n * np + j] / fabs(dtmp[j]);
                    }
                    chopmin = min(chopmin, choptmp);
                }

                // dS
                for (USI j = 0; j < np; j++) {
                    bvs.S[n * np + j] += chopmin * dtmp[j];
                }

                // dP
                bvs.P[n] += u[n * col];

                // dNi
                for (USI i = 0; i < nc; i++) {
                    bvs.Ni[n * nc + i] += chopmin * u[n * col + 1 + i];
                }
            }

            // dT
            bvs.T[n] += u[n * col + col - 1];
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

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/