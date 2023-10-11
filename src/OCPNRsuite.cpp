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
    MPI_Allreduce(&tmploc, &res.maxRelRes0_V, 1, MPI_DOUBLE, MPI_MIN, myComm);
    
    OCPTIME_COMM_COLLECTIVE += timer.Stop() / 1000;


    lP = bvs.lP;
    lT = bvs.lT;
    lN = bvs.lNi;
    lS = bvs.lS;

    dPBmaxNR.clear();
    dTmaxNR.clear();
    dNmaxNR.clear();
    dSmaxNR.clear();
    dPWmaxNR.clear();
}


void OCPNRsuite::CalMaxChangeNR(const Reservoir& rs)
{
    OCP_DBL dPmaxTmp = 0;
    OCP_DBL dTmaxTmp = 0;
    OCP_DBL dNmaxTmp = 0;
    OCP_DBL dSmaxTmp = 0;

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
    }

    lP = bvs.P;
    lT = bvs.T;
    lN = bvs.Ni;
    lS = bvs.S;

    dPBmaxNR.push_back(dPmaxTmp);
    dTmaxNR.push_back(dTmaxTmp);
    dNmaxNR.push_back(dNmaxTmp);
    dSmaxNR.push_back(dSmaxTmp);

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
}


void OCPNRsuite::CalMaxChangeTime(const Reservoir& rs)
{
    OCP_FUNCNAME;

    dPmaxT = 0;
    dPBmaxT = 0;
    dPWmaxT = 0;
    dTmaxT = 0;
    dNmaxT = 0;
    dSmaxT = 0;
    eVmaxT = 0;

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

        MPI_Allreduce(&maxCFL_loc, &maxCFL, 1, MPI_DOUBLE, MPI_MAX, myComm);

        OCPTIME_COMM_COLLECTIVE += timer.Stop() / 1000;
    }
    else {
        maxCFL = maxCFL_loc;
    }
}


OCP_BOOL OCPNRsuite::CheckCFL(const OCP_DBL& cflLim) const
{
	if (maxCFL > cflLim)
		return BULK_OUTRANGED_CFL;
	else
		return BULK_SUCCESS;
}


OCP_BOOL OCPNRsuite::CheckPhysical(Reservoir& rs, const initializer_list<string>& il) const
{
    OCP_INT workState_loc = OCP_CONTINUE;
    OCP_INT flag;
    for (auto& s : il) {

        if (s == "BulkP")        flag = rs.bulk.CheckP();
        else if (s == "BulkT")   flag = rs.bulk.CheckT();
        else if (s == "BulkNi")  flag = rs.bulk.CheckNi();
        else if (s == "BulkVe")  flag = rs.bulk.CheckVe(0.01);
        else if (s == "CFL")     flag = CheckCFL(1.0);
        else if (s == "WellP")   flag = rs.allWells.CheckP(rs.bulk);
        else                     OCP_ABORT("Check iterm not recognized!");

        switch (flag) {
            // Bulk
        case BULK_SUCCESS:
            break;

        case BULK_NEGATIVE_PRESSURE:
        case BULK_NEGATIVE_TEMPERATURE:
        case BULK_NEGATIVE_COMPONENTS_MOLES:
        case BULK_OUTRANGED_VOLUME_ERROR:
            workState_loc = OCP_RESET_CUTTIME;
            break;

        case BULK_OUTRANGED_CFL:
            workState_loc = OCP_RESET_CUTTIME_CFL;
            break;

            // Well
        case WELL_SUCCESS:
            break;

        case WELL_NEGATIVE_PRESSURE:
            workState_loc = OCP_RESET_CUTTIME;
            break;

        case WELL_SWITCH_TO_BHPMODE:
        case WELL_CROSSFLOW:
            workState_loc = OCP_RESET;
            break;

        default:
            break;
        }
        if (workState_loc != OCP_CONTINUE)
            break;
    }

    GetWallTime timer;
    timer.Start();

    MPI_Allreduce(&workState_loc, &workState, 1, MPI_INT, MPI_MIN, myComm);

    OCPTIME_COMM_COLLECTIVE += timer.Stop() / 1000;
    OCPTIME_COMM_1ALLREDUCE += timer.Stop() / 1000;

    switch (workState)
    {
    case OCP_CONTINUE:
        return OCP_TRUE;

    case OCP_RESET:
    case OCP_RESET_CUTTIME:
    case OCP_RESET_CUTTIME_CFL:
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