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
    rs.bulk.InitPTSw(50);

    InitRock(rs.bulk);
    CalRock(rs.bulk);

    InitFlash(rs.bulk);
    CalKrPc(rs.bulk);

    CalThermalConduct(rs.conn, rs.bulk);

    rs.allWells.InitBHP(rs.bulk);
    UpdateLastTimeStep(rs);
}

void T_FIM::Prepare(Reservoir& rs, const OCPControl& ctrl)
{
    rs.allWells.PrepareWell(rs.bulk);
    CalRes(rs, ctrl.GetCurTime() + ctrl.GetCurDt(), ctrl.GetCurDt(), OCP_TRUE);
}

void T_FIM::AssembleMat(LinearSystem&    ls,
                        const Reservoir& rs,
                        const OCP_DBL&   t,
                        const OCP_DBL&   dt)
{
    AssembleMatBulks(ls, rs, t, dt);
    AssembleMatWells(ls, rs, dt);
    ls.AssembleRhsCopy(res.resAbs);
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
     //ls.OutputLinearSystem("testA_FIMT.out", "testb_FIMT.out");
     //ls.OutputSolution("testx_FIMT.out");
    // Check if inf or nan occurs in solution
    ls.CheckSolution();
#endif // DEBUG
    
    timer.Start();
    GetSolution(rs, ls.GetSolution(), ctrl);
    OCPTIME_NRSTEP += timer.Stop() / 1000;
    // rs.PrintSolFIM(ctrl.workDir + "testPNi.out");
    ls.ClearData();
}

OCP_BOOL T_FIM::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    if (!ctrl.Check(rs, {"BulkNi", "BulkP", "BulkT"})) {
        ResetToLastTimeStep(rs, ctrl);
        if (CURRENT_RANK == MASTER_PROCESS) {
            cout << "Cut time step size and repeat! current dt = " << fixed
                << setprecision(3) << dt << " days\n";
        }      
        return OCP_FALSE;
    }

    // Update reservoir properties
    CalRock(rs.bulk);

    CalFlash(rs.bulk);
    CalKrPc(rs.bulk);

    CalThermalConduct(rs.conn, rs.bulk);
    CalHeatLoss(rs.bulk, ctrl.GetCurTime() + dt, dt);

    rs.allWells.CalTrans(rs.bulk);
    rs.allWells.CalFlux(rs.bulk);

    CalRes(rs, ctrl.GetCurTime() + dt, dt, OCP_FALSE);

    return OCP_TRUE;
}

OCP_BOOL T_FIM::FinishNR(Reservoir& rs, OCPControl& ctrl)
{

    NRdSmax = CalNRdSmax();
    // const OCP_DBL NRdNmax = rs.GetNRdNmax();

    OCP_INT conflag_loc = -1;
    if (((res.maxRelRes_V <= res.maxRelRes0_V * ctrl.ctrlNR.NRtol ||
        res.maxRelRes_V <= ctrl.ctrlNR.NRtol ||
        res.maxRelRes_N <= ctrl.ctrlNR.NRtol) &&
        res.maxWellRelRes_mol <= ctrl.ctrlNR.NRtol) ||
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

void T_FIM::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.CalIPRT(ctrl.GetCurDt());
    rs.CalMaxChange();
    UpdateLastTimeStep(rs);
    ctrl.CalNextTimeStep(rs, {"dP", "dS", "iter"});
    ctrl.UpdateIters();
}

void T_FIM::AllocateReservoir(Reservoir& rs)
{
    Bulk&         bk = rs.bulk;
    const OCP_USI nb = bk.nb;
    const USI     np = bk.np;
    const USI     nc = bk.nc;

    // Rock
    bk.poro.resize(nb);
    bk.rockVp.resize(nb);
    bk.vr.resize(nb);
    bk.Hr.resize(nb);

    bk.lporo.resize(nb);
    bk.lrockVp.resize(nb);
    bk.lvr.resize(nb);
    bk.lHr.resize(nb);

    // derivatives
    bk.poroP.resize(nb);
    bk.poroT.resize(nb);
    bk.vrP.resize(nb);
    bk.vrT.resize(nb);
    bk.HrT.resize(nb);

    bk.lporoP.resize(nb);
    bk.lporoT.resize(nb);
    bk.lvrP.resize(nb);
    bk.lvrT.resize(nb);
    bk.lHrT.resize(nb);

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
    bk.Uf.resize(nb);
    bk.H.resize(nb * np);
    bk.kt.resize(nb);

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
    bk.lUf.resize(nb);
    bk.lH.resize(nb * np);
    bk.lkt.resize(nb);

    // derivatives
    bk.vfP.resize(nb);
    bk.vfT.resize(nb);
    bk.vfi.resize(nb * nc);
    bk.rhoP.resize(nb * np);
    bk.rhoT.resize(nb * np);
    bk.rhox.resize(nb * nc * np);
    bk.xiP.resize(nb * np);
    bk.xiT.resize(nb * np);
    bk.xix.resize(nb * nc * np);
    bk.muP.resize(nb * np);
    bk.muT.resize(nb * np);
    bk.mux.resize(nb * nc * np);
    bk.dPcdS.resize(nb * np * np);
    bk.dKrdS.resize(nb * np * np);
    bk.UfP.resize(nb);
    bk.UfT.resize(nb);
    bk.Ufi.resize(nb * nc);
    bk.HT.resize(nb * np);
    bk.Hx.resize(nb * np * nc);
    bk.ktP.resize(nb);
    bk.ktT.resize(nb);
    bk.ktS.resize(nb * np);

    bk.lvfP.resize(nb);
    bk.lvfT.resize(nb);
    bk.lvfi.resize(nb * nc);
    bk.lrhoP.resize(nb * np);
    bk.lrhoT.resize(nb * np);
    bk.lrhox.resize(nb * nc * np);
    bk.lxiP.resize(nb * np);
    bk.lxiT.resize(nb * np);
    bk.lxix.resize(nb * nc * np);
    bk.lmuP.resize(nb * np);
    bk.lmuT.resize(nb * np);
    bk.lmux.resize(nb * nc * np);
    bk.ldPcdS.resize(nb * np * np);
    bk.ldKrdS.resize(nb * np * np);
    bk.UfP.resize(nb);
    bk.UfT.resize(nb);
    bk.Ufi.resize(nb * nc);
    bk.HT.resize(nb * np);
    bk.Hx.resize(nb * np * nc);
    bk.ktP.resize(nb);
    bk.ktT.resize(nb);
    bk.ktS.resize(nb * np);

    // FIM-Specified
    bk.maxLendSdP = (nc + 2) * (nc + 1) * np;
    bk.dSec_dPri.resize(nb * bk.maxLendSdP);

    bk.ldSec_dPri.resize(nb * bk.maxLendSdP);

    // BulkConn
    BulkConn&      conn    = rs.conn;
    const OCP_USI& numConn = conn.numConn;

    conn.bcval.upblock.resize(numConn* np);
    conn.bcval.rho.resize(numConn* np);
    conn.bcval.velocity.resize(numConn* np);


    // NR
    dSNR.resize(nb* np);
    dNNR.resize(nb* nc);
    dPNR.resize(nb);
    dTNR.resize(nb);

    // Allocate Residual
    res.SetupT(bk.nbI, rs.allWells.numWell, nc);
}

void T_FIM::AllocateLinearSystem(LinearSystem&     ls,
                                 const Reservoir&  rs,
                                 const OCPControl& ctrl)
{
    ls.SetupDomain(rs.domain);
    ls.AllocateRowMem(rs.GetComNum() + 2);
    ls.AllocateColMem();
    ls.SetupLinearSolver(THERMALMODEL, ctrl.GetWorkDir(), ctrl.GetLsFile());
}

void T_FIM::InitRock(Bulk& bk) const
{
    for (OCP_USI n = 0; n < bk.nb; n++) {
        if (bk.bType[n] == 0) {
            // non fluid bulk
            bk.poroInit[n] = 0;
            bk.poro[n]     = 0;
            bk.rockVp[n]   = 0;
            bk.vr[n]       = bk.v[n];
        }
    }
}

void T_FIM::CalRock(Bulk& bk) const
{
    const OCP_USI nb = bk.nb;

    for (OCP_USI n = 0; n < nb; n++) {
        bk.rock[bk.ROCKNUM[n]]->CalPoro(bk.P[n], bk.T[n], bk.poroInit[n], bk.bType[n]);
        if (bk.bType[n] > 0) {
            // with fluid           
            bk.poro[n]   = bk.rock[bk.ROCKNUM[n]]->GetPoro();
            bk.poroP[n]  = bk.rock[bk.ROCKNUM[n]]->GetdPorodP();
            bk.poroT[n]  = bk.rock[bk.ROCKNUM[n]]->GetdPorodT();           
            bk.vr[n]     = bk.v[n] * bk.rock[bk.ROCKNUM[n]]->Get_Poro();
            bk.vrP[n]    = bk.v[n] * bk.rock[bk.ROCKNUM[n]]->Get_dPorodP();
            bk.vrT[n]    = bk.v[n] * bk.rock[bk.ROCKNUM[n]]->Get_dPorodT();
            bk.rockVp[n] = bk.v[n] * bk.poro[n];
        }
        bk.Hr[n]  = bk.rock[bk.ROCKNUM[n]]->GetHr();
        bk.HrT[n] = bk.rock[bk.ROCKNUM[n]]->GetdHrdT();
    }
}

void T_FIM::InitFlash(Bulk& bk)
{
    const OCP_USI nb = bk.nb;
    const OCP_USI np = bk.np;
    const OCP_USI nc = bk.nc;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bType[n] > 0) {
            bk.flashCal[bk.PVTNUM[n]]->InitFlashFIM(bk.P[n], bk.Pb[n], bk.T[n],
                                                    &bk.S[n * np], bk.rockVp[n],
                                                    &bk.Ni[n * nc], n);
            for (USI i = 0; i < nc; i++) {
                bk.Ni[n * nc + i] = bk.flashCal[bk.PVTNUM[n]]->GetNi(i);
            }
            PassFlashValue(bk, n);
        }
    }
}

void T_FIM::CalFlash(Bulk& bk)
{
    const OCP_USI nb = bk.nb;
    const OCP_USI np = bk.np;
    const OCP_USI nc = bk.nc;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bType[n] > 0) {
            bk.flashCal[bk.PVTNUM[n]]->FlashFIM(bk.P[n], bk.T[n], &bk.Ni[n * nc],
                                                &bk.S[n * np], bk.phaseNum[n],
                                                &bk.xij[n * np * nc], n);
            PassFlashValue(bk, n);
        }
    }
}

void T_FIM::PassFlashValue(Bulk& bk, const OCP_USI& n)
{
    const USI     np     = bk.np;
    const USI     nc     = bk.nc;
    const OCP_USI bIdp   = n * np;
    const USI     pvtnum = bk.PVTNUM[n];

    bk.phaseNum[n] = 0;
    bk.Nt[n]       = bk.flashCal[pvtnum]->GetNt();
    bk.vf[n]       = bk.flashCal[pvtnum]->GetVf();
    bk.Uf[n]       = bk.flashCal[pvtnum]->GetUf();

    for (USI j = 0; j < np; j++) {
        // Important! Saturation must be passed no matter if the phase exists. This is
        // because it will be used to calculate relative permeability and capillary
        // pressure at each time step. Make sure that all saturations are updated at
        // each step!
        bk.S[bIdp + j] = bk.flashCal[pvtnum]->GetS(j);
        dSNR[bIdp + j] = bk.S[bIdp + j] - dSNR[bIdp + j];
        bk.phaseExist[bIdp + j] = bk.flashCal[pvtnum]->GetPhaseExist(j);
        if (bk.phaseExist[bIdp + j]) {
            bk.phaseNum[n]++;
            bk.rho[bIdp + j] = bk.flashCal[pvtnum]->GetRho(j);
            bk.xi[bIdp + j]  = bk.flashCal[pvtnum]->GetXi(j);
            bk.mu[bIdp + j]  = bk.flashCal[pvtnum]->GetMu(j);
            bk.H[bIdp + j]   = bk.flashCal[pvtnum]->GetH(j);

            // Derivatives
            bk.rhoP[bIdp + j] = bk.flashCal[pvtnum]->GetRhoP(j);
            bk.rhoT[bIdp + j] = bk.flashCal[pvtnum]->GetRhoT(j);
            bk.xiP[bIdp + j]  = bk.flashCal[pvtnum]->GetXiP(j);
            bk.xiT[bIdp + j]  = bk.flashCal[pvtnum]->GetXiT(j);
            bk.muP[bIdp + j]  = bk.flashCal[pvtnum]->GetMuP(j);
            bk.muT[bIdp + j]  = bk.flashCal[pvtnum]->GetMuT(j);
            bk.HT[bIdp + j]   = bk.flashCal[pvtnum]->GetHT(j);

            for (USI i = 0; i < nc; i++) {
                bk.xij[bIdp * nc + j * nc + i]  = bk.flashCal[pvtnum]->GetXij(j, i);
                bk.rhox[bIdp * nc + j * nc + i] = bk.flashCal[pvtnum]->GetRhoX(j, i);
                bk.xix[bIdp * nc + j * nc + i]  = bk.flashCal[pvtnum]->GetXiX(j, i);
                bk.mux[bIdp * nc + j * nc + i]  = bk.flashCal[pvtnum]->GetMuX(j, i);
                bk.Hx[bIdp * nc + j * nc + i]   = bk.flashCal[pvtnum]->GetHx(j, i);
            }
        }
    }
    bk.vfP[n] = bk.flashCal[pvtnum]->GetVfP();
    bk.vfT[n] = bk.flashCal[pvtnum]->GetVfT();
    bk.UfP[n] = bk.flashCal[pvtnum]->GetUfP();
    bk.UfT[n] = bk.flashCal[pvtnum]->GetUfT();

    for (USI i = 0; i < nc; i++) {
        bk.vfi[n * nc + i] = bk.flashCal[pvtnum]->GetVfi(i);
        bk.Ufi[n * nc + i] = bk.flashCal[pvtnum]->GetUfi(i);
    }

    Dcopy(bk.maxLendSdP, &bk.dSec_dPri[n * bk.maxLendSdP],
          &bk.flashCal[pvtnum]->GetDXsDXp()[0]);
}

void T_FIM::CalKrPc(Bulk& bk) const
{
    const USI& np = bk.np;

    for (OCP_USI n = 0; n < bk.nb; n++) {
        if (bk.bType[n] > 0) {
            OCP_USI bId = n * np;
            bk.flow[bk.SATNUM[n]]->CalKrPcFIM(&bk.S[bId], n);
            copy(bk.flow[bk.SATNUM[n]]->GetKr().begin(), bk.flow[bk.SATNUM[n]]->GetKr().end(), &bk.kr[bId]);
            copy(bk.flow[bk.SATNUM[n]]->GetPc().begin(), bk.flow[bk.SATNUM[n]]->GetPc().end(), &bk.Pc[bId]);
            copy(bk.flow[bk.SATNUM[n]]->GetdKrdS().begin(), bk.flow[bk.SATNUM[n]]->GetdKrdS().end(), &bk.dKrdS[bId * np]);
            copy(bk.flow[bk.SATNUM[n]]->GetdPcdS().begin(), bk.flow[bk.SATNUM[n]]->GetdPcdS().end(), &bk.dPcdS[bId * np]);
            for (USI j = 0; j < np; j++) bk.Pj[n * np + j] = bk.P[n] + bk.Pc[n * np + j];
        }
    }
}

void T_FIM::CalThermalConduct(BulkConn& conn, Bulk& bk) const
{
    const OCP_USI nb = bk.nb;
    const OCP_USI np = bk.np;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bType[n] > 0) {
            // fluid bulk
            OCP_DBL tmp = 0;
            for (USI j = 0; j < np; j++) {
                tmp += bk.S[n * np + j] * bk.thconp[j];
                bk.ktS[n * np + j] = bk.poro[n] * bk.thconp[j];
            }
            bk.kt[n]  = bk.poro[n] * tmp + (1 - bk.poro[n]) * bk.thconr[n];
            bk.ktP[n] = bk.poroP[n] * (tmp - bk.thconr[n]);
            bk.ktT[n] = bk.poroT[n] * (tmp - bk.thconr[n]);
        } else {
            // non fluid bulk
            bk.kt[n] = bk.thconr[n];
        }
    }

    //OCP_USI bId, eId;
    //OCP_DBL areaB, areaE, T1, T2;
    //OCP_DBL tmpB, tmpE;

    //for (OCP_USI c = 0; c < conn.numConn; c++) {
    //    bId = conn.iteratorConn[c].BId();
    //    eId = conn.iteratorConn[c].EId();
    //    if (bk.bType[bId] > 0 && bk.bType[eId] > 0) {
    //        // fluid bulk connections
    //        areaB        = conn.iteratorConn[c].AreaB();
    //        areaE        = conn.iteratorConn[c].AreaE();
    //        T1           = bk.kt[bId] * areaB;
    //        T2           = bk.kt[eId] * areaE;
    //        conn.Adkt[c] = 1 / (1 / T1 + 1 / T2);

    //        tmpB                  = pow(conn.Adkt[c], 2) / pow(T1, 2) * areaB;
    //        tmpE                  = pow(conn.Adkt[c], 2) / pow(T2, 2) * areaE;
    //        conn.AdktP[c * 2 + 0] = tmpB * bk.ktP[bId];
    //        conn.AdktP[c * 2 + 1] = tmpE * bk.ktP[eId];
    //        conn.AdktT[c * 2 + 0] = tmpB * bk.ktT[bId];
    //        conn.AdktT[c * 2 + 1] = tmpE * bk.ktT[eId];
    //        for (USI j = 0; j < np; j++) {
    //            conn.AdktS[c * np * 2 + j]      = tmpB * bk.ktS[bId * np + j];
    //            conn.AdktS[c * np * 2 + np + j] = tmpE * bk.ktS[eId * np + j];
    //        }
    //    }
    //}
}

void T_FIM::CalHeatLoss(Bulk& bk, const OCP_DBL& t, const OCP_DBL& dt) const
{
    bk.hLoss.CalHeatLoss(bk.bLocation, bk.T, bk.lT, bk.initT, t, dt);
}

void T_FIM::ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl)
{
    // Bulk
    Bulk& bk = rs.bulk;

    // Rock
    bk.poro   = bk.lporo;
    bk.rockVp = bk.lrockVp;
    bk.vr     = bk.lvr;
    bk.Hr     = bk.lHr;
    // derivatives
    bk.poroP = bk.lporoP;
    bk.poroT = bk.lporoT;
    bk.vrP   = bk.lvrP;
    bk.vrT   = bk.lvrT;
    bk.HrT   = bk.lHrT;

    // Fluid
    bk.phaseNum   = bk.lphaseNum;
    bk.Nt         = bk.lNt;
    bk.Ni         = bk.lNi;
    bk.vf         = bk.lvf;
    bk.T          = bk.lT;
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
    bk.Uf         = bk.lUf;
    bk.H          = bk.lH;
    bk.kt         = bk.lkt;
    // derivatives
    bk.vfP       = bk.lvfP;
    bk.vfT       = bk.lvfT;
    bk.vfi       = bk.lvfi;
    bk.rhoP      = bk.lrhoP;
    bk.rhoT      = bk.lrhoT;
    bk.rhox      = bk.lrhox;
    bk.xiP       = bk.lxiP;
    bk.xiT       = bk.lxiT;
    bk.xix       = bk.lxix;
    bk.muP       = bk.lmuP;
    bk.muT       = bk.lmuT;
    bk.mux       = bk.lmux;
    bk.dPcdS   = bk.ldPcdS;
    bk.dKrdS    = bk.ldKrdS;
    bk.UfP       = bk.lUfP;
    bk.UfT       = bk.lUfT;
    bk.Ufi       = bk.lUfi;
    bk.HT        = bk.lHT;
    bk.Hx        = bk.lHx;
    bk.ktP       = bk.lktP;
    bk.ktT       = bk.lktT;
    bk.ktS       = bk.lktS;
    bk.dSec_dPri = bk.ldSec_dPri;

    bk.hLoss.ResetToLastTimeStep();

    // Wells
    rs.allWells.ResetBHP();
    rs.allWells.CalTrans(bk);
    rs.allWells.CaldG(bk);
    rs.allWells.CalFlux(bk);

    rs.optFeatures.ResetToLastTimeStep();

    // Iters
    ctrl.ResetIterNRLS();

    CalRes(rs, ctrl.GetCurTime() + ctrl.GetCurDt(), ctrl.GetCurDt(), OCP_TRUE);
}

void T_FIM::UpdateLastTimeStep(Reservoir& rs) const
{
    // Bulk
    Bulk& bk = rs.bulk;

    // Rock
    bk.lporo   = bk.poro;
    bk.lrockVp = bk.rockVp;
    bk.lvr     = bk.vr;
    bk.lHr     = bk.Hr;
    // derivatives
    bk.lporoP = bk.poroP;
    bk.lporoT = bk.poroT;
    bk.lvrP   = bk.vrP;
    bk.lvrT   = bk.vrT;
    bk.lHrT   = bk.HrT;

    // Fluid
    bk.lphaseNum   = bk.phaseNum;
    bk.lNt         = bk.Nt;
    bk.lNi         = bk.Ni;
    bk.lvf         = bk.vf;
    bk.lT          = bk.T;
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
    bk.lUf         = bk.Uf;
    bk.lH          = bk.H;
    bk.lkt         = bk.kt;
    // derivatives
    bk.lvfP       = bk.vfP;
    bk.lvfT       = bk.vfT;
    bk.lvfi       = bk.vfi;
    bk.lrhoP      = bk.rhoP;
    bk.lrhoT      = bk.rhoT;
    bk.lrhox      = bk.rhox;
    bk.lxiP       = bk.xiP;
    bk.lxiT       = bk.xiT;
    bk.lxix       = bk.xix;
    bk.lmuP       = bk.muP;
    bk.lmuT       = bk.muT;
    bk.lmux       = bk.mux;
    bk.ldPcdS   = bk.dPcdS;
    bk.ldKrdS    = bk.dKrdS;
    bk.lUfP       = bk.UfP;
    bk.lUfT       = bk.UfT;
    bk.lUfi       = bk.Ufi;
    bk.lHT        = bk.HT;
    bk.lHx        = bk.Hx;
    bk.lktP       = bk.ktP;
    bk.lktT       = bk.ktT;
    bk.lktS       = bk.ktS;
    bk.ldSec_dPri = bk.dSec_dPri;

    bk.hLoss.UpdateLastTimeStep();

    rs.allWells.UpdateLastTimeStepBHP();
    rs.optFeatures.UpdateLastTimeStep();
}

void T_FIM::CalRes(Reservoir&      rs,
                   const OCP_DBL&  t,
                   const OCP_DBL&  dt,
                   const OCP_BOOL& resetRes0)
{
    const Bulk& bk   = rs.bulk;
    const USI   nb   = bk.nbI;
    const USI   np   = bk.np;
    const USI   nc   = bk.nc;
    const USI   len  = nc + 2;
    
    res.SetZero();

    // Bulk to Bulk
    OCP_USI bId, eId, bIdb;
    // Accumalation Term
    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bType[n] > 0) {
            // Fluid bulk
            bId  = n * len;
            bIdb = n * nc;
            // Volume Conservation
            res.resAbs[bId] = bk.rockVp[n] - bk.vf[n];
            // Mass Conservation
            for (USI i = 0; i < nc; i++) {
                res.resAbs[n * len + 1 + i] = bk.Ni[bIdb + i] - bk.lNi[bIdb + i];
            }
            // Energy Conservation
            res.resAbs[n * len + nc + 1] =
                (bk.vf[n] * bk.Uf[n] + bk.vr[n] * bk.Hr[n]) -
                (bk.lvf[n] * bk.lUf[n] + bk.lvr[n] * bk.lHr[n]);
        } else {
            // Non fluid bulk
            res.resAbs[n * len + nc + 1] = bk.vr[n] * bk.Hr[n] - bk.lvr[n] * bk.lHr[n];
        }

        // Heat Loss
        if (bk.hLoss.IfHeatLoss() && bk.bLocation[n] > 0) {
            const OCP_DBL lambda = bk.bLocation[n] == 1 ? bk.hLoss.obD : bk.hLoss.ubD;
            const OCP_DBL kappa  = bk.bLocation[n] == 1 ? bk.hLoss.obK : bk.hLoss.ubK;
            // dT
            res.resAbs[n * len + nc + 1] +=
                dt * kappa *
                (2 * (bk.T[n] - bk.initT[n]) / sqrt(lambda * t) - bk.hLoss.p[n]) *
                bk.dx[n] * bk.dy[n];
        }
    }

    BulkConn&         conn  = rs.conn;
    BulkConnVal&      bcval = conn.bcval;
    vector<OCPFlux*>& flux  = conn.flux;

    OCP_USI uId_np_j;
    OCP_DBL dT, Akdt;
    USI     cType;

    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId   = conn.iteratorConn[c].BId();
        eId   = conn.iteratorConn[c].EId();
        cType = conn.iteratorConn[c].Type();

        flux[cType]->CalFlux(conn.iteratorConn[c], bk);

        // Thermal conductive term always exists
        Akdt = flux[cType]->GetAdkt();
        dT   = bk.T[bId] - bk.T[eId];
        res.resAbs[bId * len + 1 + nc] += Akdt * dT * dt;
        if (eId < nb) {
            // Interior grid
            res.resAbs[eId * len + 1 + nc] -= Akdt * dT * dt;
        }

        if (cType == 0) {
            // with flow
            copy(flux[cType]->GetUpblock().begin(), flux[cType]->GetUpblock().end(), &bcval.upblock[c * np]);
            copy(flux[cType]->GetRho().begin(), flux[cType]->GetRho().end(), &bcval.rho[c * np]);
            copy(flux[cType]->GetFluxVj().begin(), flux[cType]->GetFluxVj().end(), &bcval.velocity[c * np]);

			if (eId < nb) {
				// Interior grid
				for (USI i = 0; i < nc; i++) {
					res.resAbs[bId * len + 1 + i] += dt * flux[cType]->GetFluxNi()[i];
					res.resAbs[eId * len + 1 + i] -= dt * flux[cType]->GetFluxNi()[i];
				}
                for (USI j = 0; j < np; j++) {
                    uId_np_j = bcval.upblock[c * np + j] * np + j;
                    res.resAbs[bId * len + 1 + nc] += dt * bcval.velocity[c * np + j] * bk.xi[uId_np_j] * bk.H[uId_np_j];
                    res.resAbs[eId * len + 1 + nc] -= dt * bcval.velocity[c * np + j] * bk.xi[uId_np_j] * bk.H[uId_np_j];
                }
			}
			else {
				// Ghost grid
				for (USI i = 0; i < nc; i++) {
					res.resAbs[bId * len + 1 + i] += dt * flux[cType]->GetFluxNi()[i];
				}
                for (USI j = 0; j < np; j++) {
                    uId_np_j = bcval.upblock[c * np + j] * np + j;
                    res.resAbs[bId * len + 1 + nc] += dt * bcval.velocity[c * np + j] * bk.xi[uId_np_j] * bk.H[uId_np_j];
                }				
			}
        }
    }

    // Well to Bulk
    USI wId = nb * len;
    for (const auto& wl : rs.allWells.wells) {
        wl.CalResFIM_T(wId, res, bk, dt);
    }

    // Calculate RelRes
    OCP_DBL tmp;
    for (OCP_USI n = 0; n < nb; n++) {

        // Energy equations always exist
        OCP_DBL eT = bk.vf[n] * bk.Uf[n] + bk.vr[n] * bk.Hr[n];
        tmp        = fabs(res.resAbs[n * len + nc + 1] / eT);
        if (res.maxRelRes_E < tmp) {
            res.maxRelRes_E = tmp;
            res.maxId_E     = n;
        }

        if (bk.bType[n] > 0) {
            // Fluid Bulk
            for (USI i = 0; i < len - 1; i++) {
                tmp = fabs(res.resAbs[n * len + i] / bk.rockVp[n]);
                if (res.maxRelRes_V < tmp) {
                    res.maxRelRes_V = tmp;
                    res.maxId_V     = n;
                }
            }

            for (USI i = 1; i < len - 1; i++) {
                tmp = fabs(res.resAbs[n * len + i] / bk.Nt[n]);
                if (res.maxRelRes_N < tmp) {
                    res.maxRelRes_N = tmp;
                    res.maxId_N     = n;
                }
            }
        }
    }

    Dscalar(res.resAbs.size(), -1.0, res.resAbs.data());
    if (resetRes0) {
        res.SetInitRes();

        GetWallTime timer;
        timer.Start();

        OCP_DBL tmploc = res.maxRelRes0_V;
        MPI_Allreduce(&tmploc, &res.maxRelRes0_V, 1, MPI_DOUBLE, MPI_MIN, rs.domain.myComm);

        OCPTIME_COMM_COLLECTIVE += timer.Stop() / 1000;
    }
}

void T_FIM::AssembleMatBulks(LinearSystem&    ls,
                             const Reservoir& rs,
                             const OCP_DBL&   t,
                             const OCP_DBL&   dt) const
{
    const USI numWell = rs.GetNumOpenWell();

    const Bulk&     bk     = rs.bulk;
    const BulkConn& conn   = rs.conn;
    const OCP_USI   nb     = bk.nbI;
    const USI       np     = bk.np;
    const USI       nc     = bk.nc;
    const USI       ncol   = nc + 2;
    const USI       ncol2  = np * nc + np;
    const USI       bsize  = ncol * ncol;
    const USI       bsize2 = ncol * ncol2;

    ls.AddDim(nb);

    vector<OCP_DBL> bmat(bsize, 0);
    // Accumulation term
    for (USI i = 1; i < nc + 1; i++) {
        // Mass consevation
        bmat[i * ncol + i] = 1;
    }
    vector<OCP_DBL> bmatR(bmat);
    bmatR[0] = 1;
    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bType[n] > 0) {
            // Fluid Bulk
            // Volume consevation
            // dP
            bmat[0] = bk.v[n] * bk.poroP[n] - bk.vfP[n];
            // dNi
            for (USI i = 0; i < nc; i++) {
                bmat[i + 1] = -bk.vfi[n * nc + i];
            }
            // dT
            bmat[nc + 1] = bk.v[n] * bk.poroT[n] - bk.vfT[n];
            // Energy consevation
            // dP
            bmat[ncol * (ncol - 1)] =
                bk.vfP[n] * bk.Uf[n] + bk.vf[n] * bk.UfP[n] + bk.vrP[n] * bk.Hr[n];
            // dNi
            for (USI i = 0; i < nc; i++) {
                bmat[ncol * (ncol - 1) + i + 1] =
                    bk.vfi[n * nc + i] * bk.Uf[n] + bk.vf[n] * bk.Ufi[n * nc + i];
            }
            // dT
            bmat[ncol * ncol - 1] = bk.vfT[n] * bk.Uf[n] + bk.vf[n] * bk.UfT[n] +
                                    bk.vrT[n] * bk.Hr[n] + bk.vr[n] * bk.HrT[n];

            // Heat Loss iterm
            if (bk.hLoss.IfHeatLoss() && bk.bLocation[n] > 0) {
                const OCP_DBL lambda =
                    bk.bLocation[n] == 1 ? bk.hLoss.obD : bk.hLoss.ubD;
                const OCP_DBL kappa =
                    bk.bLocation[n] == 1 ? bk.hLoss.obK : bk.hLoss.ubK;
                // dT
                bmat[ncol * ncol - 1] += dt * kappa *
                                         (2 / sqrt(lambda * t) - bk.hLoss.pT[n]) *
                                         bk.dx[n] * bk.dy[n];
            }

            ls.NewDiag(n, bmat);
        } else {
            // Non Fluid Bulk
            // Energy consevation
            // dT
            bmatR[ncol * ncol - 1] = bk.vrT[n] * bk.Hr[n] + bk.vr[n] * bk.HrT[n];

            // Heat Loss iterm
            if (bk.hLoss.IfHeatLoss() && bk.bLocation[n] > 0) {
                const OCP_DBL lambda =
                    bk.bLocation[n] == 1 ? bk.hLoss.obD : bk.hLoss.ubD;
                const OCP_DBL kappa =
                    bk.bLocation[n] == 1 ? bk.hLoss.obK : bk.hLoss.ubK;
                // dT
                bmatR[ncol * ncol - 1] += dt * kappa *
                                          (2 / sqrt(lambda * t) - bk.hLoss.pT[n]) *
                                          bk.dx[n] * bk.dy[n];
            }

            ls.NewDiag(n, bmatR);
        }
    }

    // flux term
    OCP_USI  bId, eId;
    USI      cType;

    vector<OCPFlux*>& flux = conn.flux;
    for (OCP_USI c = 0; c < conn.numConn; c++) {

        bId   = conn.iteratorConn[c].BId();
        eId   = conn.iteratorConn[c].EId();
        cType = conn.iteratorConn[c].Type();

        flux[cType]->AssembleMatFIM(conn.iteratorConn[c], c, conn.bcval, bk);

        bmat = flux[cType]->GetdFdXpB();
        DaABpbC(ncol, ncol, ncol2, 1, flux[cType]->GetdFdXsB().data(), &bk.dSec_dPri[bId * bsize2], 1,
            bmat.data());

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
            OCP_ABORT("INF or INF in bmat !");
        }
#endif

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
        wl.AssembleMatFIM_T(ls, rs.bulk, dt);
    }
}


void T_FIM::GetSolution(Reservoir&             rs,
                        vector<OCP_DBL>& u,
                        const OCPControl&      ctrl)
{

    const Domain&   domain = rs.domain;
    Bulk&           bk     = rs.bulk;
    const OCP_USI   nb     = bk.nb;
    const USI       np     = bk.np;
    const USI       nc     = bk.nc;
    const USI       row    = np * (nc + 1);
    const USI       col    = nc + 2;

    // Well first
    USI wId = bk.nbI * col;
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

    dSNR    = bk.S;
    NRdPmax = 0;
    NRdNmax = 0;
    NRdTmax = 0;

    OCP_USI bId = 0;
    OCP_USI eId = bk.GetInteriorBulkNum();

    // interior first, ghost second
    for (USI p = bId; p < eId; p++) {

        for (OCP_USI n = bId; n < eId; n++) {
            // const vector<OCP_DBL>& scm = satcm[SATNUM[n]];

            if (bk.bType[n] > 0) {
                // Fluid Bulk

                chopmin = 1;
                // compute the chop
                fill(dtmp.begin(), dtmp.end(), 0.0);
                DaAxpby(np, col, 1, &bk.dSec_dPri[n * bk.maxLendSdP], u.data() + n * col, 1,
                    dtmp.data());

                for (USI j = 0; j < np; j++) {
                    choptmp = 1;
                    if (fabs(dtmp[j]) > dSmaxlim) {
                        choptmp = dSmaxlim / fabs(dtmp[j]);
                    }
                    else if (bk.S[n * np + j] + dtmp[j] < 0.0) {
                        choptmp = 0.9 * bk.S[n * np + j] / fabs(dtmp[j]);
                    }
                    chopmin = min(chopmin, choptmp);
                }

                // dS
                for (USI j = 0; j < np; j++) {
                    bk.S[n * np + j] += chopmin * dtmp[j];
                }

                // dP
                OCP_DBL dP = u[n * col];
                if (fabs(NRdPmax) < fabs(dP)) NRdPmax = dP;
                bk.P[n] += dP;
                dPNR[n] = dP;

                // dNi
                for (USI i = 0; i < nc; i++) {
                    dNNR[n * nc + i] = u[n * col + 1 + i] * chopmin;
                    if (fabs(NRdNmax) < fabs(dNNR[n * nc + i]) / bk.Nt[n])
                        NRdNmax = dNNR[n * nc + i] / bk.Nt[n];

                    bk.Ni[n * nc + i] += dNNR[n * nc + i];
                }
            }

            // dT
            OCP_DBL dT = u[n * col + col - 1];
            if (fabs(NRdTmax) < fabs(dT)) NRdTmax = dT;
            bk.T[n] += dT;
            dTNR[n] = dT;
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