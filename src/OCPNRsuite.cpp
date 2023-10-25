/*! \file    OCPNRsuite.cpp
 *  \brief   data structure used in non-linear iterations
 *  \author  Shizhe Li
 *  \date    Oct/07/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


 // Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "OCPNRsuite.hpp"


void OCPNRsuite::Setup(const OCP_BOOL& ifthermal, const BulkVarSet& bvs, const OCP_USI& nw, const Domain& domain)
{
    myComm  = domain.myComm;
    numproc = domain.numproc;
    myrank  = domain.myrank;

    nb        = bvs.nbI;
    np        = bvs.np;
    nc        = bvs.nc;

    ifUseNR   = OCP_TRUE;
    ifThermal = ifthermal;

    if (ifThermal) {
        res.SetupT(nb, nw, nc);
    }
    else {
        res.SetupIsoT(nb, nw, nc);
    }


    lP.resize(nb);
    lT.resize(nb);
    lN.resize(nb * nc);
    lS.resize(nb * np);
    dP.resize(nb);
    dT.resize(nb);
    dN.resize(nb * nc);
    dS.resize(nb * np);

    cfl.resize(nb * np);
}


void OCPNRsuite::Setup(const BulkVarSet& bvs, const Domain& domain)
{
    ifUseNR = OCP_FALSE;

    myComm  = domain.myComm;
    numproc = domain.numproc;
    myrank  = domain.myrank;

    nb = bvs.nbI;
    np = bvs.np;
    nc = bvs.nc;

    cfl.resize(nb * np);
}


void OCPNRsuite::InitStep(const BulkVarSet& bvs)
{
    GetWallTime timer;
    timer.Start();
    
    OCP_DBL tmploc = res.maxRelRes_V;
    MPI_Allreduce(&tmploc, &res.maxRelRes0_V, 1, OCPMPI_DBL, MPI_MIN, myComm);
    
    OCPTIME_COMM_COLLECTIVE += timer.Stop() / TIME_S2MS;


    lP = bvs.lP;
    lT = bvs.lT;
    lN = bvs.lNi;
    lS = bvs.lS;

    dPBmaxNR.clear();
    dPWmaxNR.clear();
    dTmaxNR.clear();
    dNmaxNR.clear();
    dSmaxNR.clear();
    eVmaxNR.clear();
}


void OCPNRsuite::CalMaxChangeNR(const Reservoir& rs)
{
    OCP_DBL dPmaxTmp = 0;
    OCP_DBL dTmaxTmp = 0;
    OCP_DBL dNmaxTmp = 0;
    OCP_DBL dSmaxTmp = 0;
    OCP_DBL eVmaxTmp = 0;

    // for bulk
    const BulkVarSet& bvs = rs.bulk.GetVarSet();

    for (OCP_USI n = 0; n < nb; n++) {
        dP[n] = bvs.P[n] - lP[n];
        if (fabs(dPmaxTmp) < fabs(dP[n]))  dPmaxTmp = dP[n];
        dT[n] = bvs.T[n] - lT[n];
        if (fabs(dTmaxTmp) < fabs(dT[n]))  dTmaxTmp = dT[n];

        OCP_DBL  Nt   = 0;
        for (USI i = 0; i < nc; i++) {
            Nt             += lN[n * nc + i];
            dN[n * nc + i] = bvs.Ni[n * nc + i] - lN[n * nc + i];
        }
        for (USI i = 0; i < nc; i++) {
            if (fabs(dNmaxTmp) < fabs(dN[n * nc + i] / Nt))  dNmaxTmp = dN[n * nc + i] / Nt;
        }

        for (USI j = 0; j < np; j++) {
            dS[n * np + j] = bvs.S[n * np + j] - lS[n * np + j];
            if (fabs(dSmaxTmp) < fabs(dS[n * np + j]))  dSmaxTmp = dS[n * np + j];
        }

        if (fabs(eVmaxTmp) < fabs((bvs.vf[n] - bvs.rockVp[n]) / bvs.rockVp[n])) {
            eVmaxTmp = (bvs.vf[n] - bvs.rockVp[n]) / bvs.rockVp[n];
        }
    }

    lP = bvs.P;
    lT = bvs.T;
    lN = bvs.Ni;
    lS = bvs.S;

    dPBmaxNR.push_back(dPmaxTmp);
    dTmaxNR.push_back(dTmaxTmp);
    dNmaxNR.push_back(dNmaxTmp);
    dSmaxNR.push_back(dSmaxTmp);
    eVmaxNR.push_back(eVmaxTmp);

    // for well   -- wrong now
    dPmaxTmp = 0;
    const auto& wells = rs.allWells.wells;
    for (const auto& w : wells) {
        const OCP_DBL dPw = w->CalMaxChangeNR();
        if (fabs(dPmaxTmp) < fabs(dPw)) {
            dPmaxTmp = dPw;
        }
    }
    dPWmaxNR.push_back(dPmaxTmp);


    //cout << scientific << setprecision(6) << res.maxRelRes0_V << "   " << res.maxRelRes_V << "   ";
    //cout << "dPB: " << dPBmaxNR.back() << "   dPW: " << dPWmaxNR.back()
    //    << "    dS : " << dSmaxNR.back() << endl;
}


void OCPNRsuite::CalMaxChangeTime(const Reservoir& rs)
{
    OCP_FUNCNAME;

    dPmaxT  = 0;
    dPBmaxT = 0;
    dPWmaxT = 0;
    dTmaxT  = 0;
    dNmaxT  = 0;
    dSmaxT  = 0;
    eVmaxT  = 0;

    // for bulk
    const BulkVarSet& bvs = rs.bulk.GetVarSet();

    OCP_USI id;
    for (OCP_USI n = 0; n < nb; n++) {

        // dP
        if (fabs(dPBmaxT) < fabs(bvs.P[n] - bvs.lP[n])) {
            dPBmaxT = bvs.P[n] - bvs.lP[n];
        }

        // dT
        if (fabs(dTmaxT) < fabs(bvs.T[n] - bvs.lT[n])) {
            dTmaxT = bvs.T[n] - bvs.lT[n];
        }

        // dS
        for (USI j = 0; j < np; j++) {
            id = n * np + j;
            if (fabs(dSmaxT) < fabs(bvs.S[id] - bvs.lS[id])) {
                dSmaxT = bvs.S[id] - bvs.lS[id];
            }
        }

        // dN
        for (USI i = 0; i < nc; i++) {
            id = n * nc + i;
            const OCP_DBL tmp = fabs(max(bvs.Ni[id], bvs.lNi[id]));
            if (tmp > TINY) {
                if (fabs(dNmaxT) < fabs((bvs.Ni[id] - bvs.lNi[id]) / tmp)) {
                    dNmaxT = (bvs.Ni[id] - bvs.lNi[id]) / tmp;
                }
            }
        }

        // Ve
        if (fabs(eVmaxT) < fabs((bvs.vf[n] - bvs.rockVp[n]) / bvs.rockVp[n])) {
            eVmaxT = (bvs.vf[n] - bvs.rockVp[n]) / bvs.rockVp[n];
        }
    }

    // for well
    const auto& wells = rs.allWells.wells;
    for (const auto& w : wells) {
        OCP_DBL dPw = w->CalMaxChangeTime();
        if (fabs(dPWmaxT) < fabs(dPw)) {
            dPWmaxT = dPw;
        }
    }

    if (fabs(dPBmaxT) < fabs(dPWmaxT)) {
        dPmaxT = dPWmaxT;
    }
    else {
        dPmaxT = dPBmaxT;
    }
}


/// Calculate CFL number
void OCPNRsuite::CalCFL(const Reservoir& rs, const OCP_DBL& dt, const OCP_BOOL& ifComm)
{
    const BulkVarSet&     bvs     = rs.bulk.GetVarSet();
    const BulkConnVarSet& cvs     = rs.conn.GetVarSet();
    const vector<Well*>&  wells   = rs.allWells.wells;
    const auto            numConn = rs.conn.GetNumConn();

    fill(cfl.begin(), cfl.end(), 0.0);

    for (OCP_USI c = 0; c < numConn; c++) {
        for (USI j = 0; j < np; j++) {
            const OCP_USI uId = cvs.upblock[c * np + j];
            if (uId < nb) {
                cfl[uId * np + j] += fabs(cvs.flux_vj[c * np + j]) * dt;
            }          
        }
    }

    for (const auto& wl : wells) {
        if (wl->IsOpen() && wl->WellType() == WellType::productor) {
            for (USI p = 0; p < wl->PerfNum(); p++) {
                if (wl->PerfState(p) == WellState::open) {
                    const OCP_USI k = wl->PerfLocation(p);

                    for (USI j = 0; j < np; j++) {
                        cfl[k * np + j] += fabs(wl->PerfProdQj_ft3(p, j)) * dt;
                    }
                }
            }
        }
    }

    maxCFL_loc = 0;
    const OCP_USI len = nb * np;
    for (OCP_USI n = 0; n < len; n++) {
        if (bvs.phaseExist[n] && bvs.vj[n] > TINY) {
            cfl[n] /= bvs.vj[n];
#ifdef DEBUG
            if (!isfinite(cfl[n])) {
                OCP_ABORT("cfl is nan!");
            }
#endif // DEBUG
            if (maxCFL_loc < cfl[n]) maxCFL_loc = cfl[n];
        }
    }
    if (ifComm) {

        GetWallTime timer;
        timer.Start();

        MPI_Allreduce(&maxCFL_loc, &maxCFL, 1, OCPMPI_DBL, MPI_MAX, myComm);

        OCPTIME_COMM_COLLECTIVE += timer.Stop() / TIME_S2MS;
    }
    else {
        maxCFL = maxCFL_loc;
    }
}


ReservoirState OCPNRsuite::CheckCFL(const OCP_DBL& cflLim) const
{
    if (maxCFL > cflLim)
        return ReservoirState::bulk_large_CFL;
	else
        return ReservoirState::bulk_success;
}


OCP_BOOL OCPNRsuite::CheckPhysical(Reservoir& rs, const initializer_list<string>& il) const
{
    OCPNRStateP     workState_loc = OCPNRStateP::continueSol;
    ReservoirState rsState;
    for (auto& s : il) {

        if (s == "BulkP")        rsState = rs.bulk.CheckP();
        else if (s == "BulkT")   rsState = rs.bulk.CheckT();
        else if (s == "BulkNi")  rsState = rs.bulk.CheckNi();
        else if (s == "BulkVe")  rsState = rs.bulk.CheckVe(0.01);
        else if (s == "CFL")     rsState = CheckCFL(1.0);
        else if (s == "WellP")   rsState = rs.allWells.CheckP(rs.bulk);
        else                     OCP_ABORT("Check iterm not recognized!");

        switch (rsState) {
            // Bulk
        case ReservoirState::bulk_success:
            break;

        case ReservoirState::bulk_negative_P:
        case ReservoirState::bulk_negative_T:
        case ReservoirState::bulk_negative_N:
        case ReservoirState::bulk_large_EV:
            workState_loc = OCPNRStateP::resetCut;
            break;

        case ReservoirState::bulk_large_CFL:
            workState_loc = OCPNRStateP::resetCutCFL;
            break;

            // Well
        case ReservoirState::well_success:
            break;

        case ReservoirState::well_negative_P:
            workState_loc = OCPNRStateP::resetCut;
            break;

        case ReservoirState::well_switch_BHPm:
        case ReservoirState::well_crossflow:
            workState_loc = OCPNRStateP::reset;
            break;

        default:
            break;
        }
        if (workState_loc != OCPNRStateP::continueSol)
            break;
    }

    GetWallTime timer;
    timer.Start();

    MPI_Allreduce(&workState_loc, &workState, 1, OCPMPI_ENUM, MPI_MAX, myComm);

    OCPTIME_COMM_COLLECTIVE += timer.Stop() / TIME_S2MS;
    OCPTIME_COMM_1ALLREDUCE += timer.Stop() / TIME_S2MS;

    switch (workState)
    {
    case OCPNRStateP::continueSol:
        return OCP_TRUE;

    case OCPNRStateP::reset:
    case OCPNRStateP::resetCut:
    case OCPNRStateP::resetCutCFL:
        return OCP_FALSE;

    default:
        OCP_ABORT("WRONG work state!");
    }
}


void OCPNRsuite::InitIter() {
    iterNR  = 0;
    iterLS  = 0;
    iterNRw = 0;
    iterLSw = 0;
    iterNRLS.clear();
}


void OCPNRsuite::UpdateIter(const USI& lsIter) 
{
    iterNR++;
    iterLS += lsIter;
    iterNRLS.push_back(lsIter);
}


void OCPNRsuite::ResetIter() 
{
    iterNRw += iterNR;
    iterLSw += iterLS;
    iterNR = 0;
    iterLS = 0;
    iterNRLS.clear();
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/07/2023      Create file                          */
/*----------------------------------------------------------------------------*/