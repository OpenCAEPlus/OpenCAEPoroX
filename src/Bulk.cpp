/*! \file    Bulk.cpp
 *  \brief   Bulk class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include <algorithm>
#include <cmath>
#include <ctime>

// OpenCAEPoroX header files
#include "Bulk.hpp"

/////////////////////////////////////////////////////////////////////
// HeatLoss
/////////////////////////////////////////////////////////////////////

void HeatLoss::InputParam(const HLoss& loss)
{
    ifHLoss = loss.ifHLoss;
    if (ifHLoss) {
        obC = loss.obC;
        obK = loss.obK;
        ubC = loss.ubC;
        ubK = loss.ubK;
    }
}

void HeatLoss::Setup(const OCP_USI& numBulk)
{
    if (ifHLoss) {
        obD     = obK / obC;
        ubD     = ubK / ubC;
        nb = numBulk;
        I.resize(nb);
        p.resize(nb);
        pT.resize(nb);
        lI.resize(nb);
        lp.resize(nb);
        lpT.resize(nb);
    }
}

void HeatLoss::CalHeatLoss(const vector<USI>&     location,
                           const vector<OCP_DBL>& T,
                           const vector<OCP_DBL>& lT,
                           const vector<OCP_DBL>& initT,
                           const OCP_DBL&         t,
                           const OCP_DBL&         dt)
{
    if (ifHLoss) {
        OCP_DBL lambda, d, dT, theta;
        OCP_DBL tmp, q;
        for (OCP_USI n = 0; n < nb; n++) {
            if (location[n] > 0) {
                // overburden or underburden
                lambda = location[n] == 1 ? obD : ubD;

                dT    = T[n] - lT[n];
                theta = T[n] - initT[n];
                d     = sqrt(lambda * t) / 2;
                tmp   = 3 * pow(d, 2) + lambda * dt;
                p[n]  = (theta * (lambda * dt / d) + lI[n] -
                        dT * (pow(d, 3) / (lambda * dt))) /
                       tmp;
                pT[n] = (lambda * dt / d - pow(d, 3) / (lambda * dt)) / tmp;
                q     = (2 * p[n] * d - theta + pow(d, 2) * dT / (lambda * dt)) /
                    (2 * pow(d, 2));
                I[n] = theta * d + p[n] * pow(d, 2) + 2 * q * pow(d, 3);
            }
        }
    }
}

void HeatLoss::ResetToLastTimeStep()
{
    I  = lI;
    p  = lp;
    pT = lpT;
}

void HeatLoss::UpdateLastTimeStep()
{
    lI  = I;
    lp  = p;
    lpT = pT;
}

/////////////////////////////////////////////////////////////////////
// Input Param and Setup
/////////////////////////////////////////////////////////////////////

/// Read parameters from rs_param data structure.
void Bulk::InputParam(const ParamReservoir& rs_param, OptionalFeatures& opts)
{
    OCP_FUNCNAME;


    if (rs_param.PBVD_T.data.size() > 0) EQUIL.PBVD.Setup(rs_param.PBVD_T.data[0]);

    if (rs_param.thermal) {
        // ifThermal model
        InputParamTHERMAL(rs_param, opts);
    }
    else if (rs_param.blackOil) {
        // Isothermal blackoil model
        InputParamBLKOIL(rs_param, opts);
    } 
    else if (rs_param.comps) {
        // Isothermal compositional model
        InputParamCOMPS(rs_param, opts);
    }


    PVTm.Setup(rs_param, vs, opts);
    SATm.Setup(rs_param, vs, PVTm.GetMixtureType(), opts);
    ROCKm.Setup(rs_param, vs, opts);
    hLoss.InputParam(rs_param.hLoss);
}

void Bulk::InputParamBLKOIL(const ParamReservoir& rs_param, OptionalFeatures& opts)
{

    if (rs_param.water && !rs_param.oil && !rs_param.gas) {
        // water
        OCP_ABORT("NOT COMPLETED!");
    } 
    else if (rs_param.water && rs_param.oil && !rs_param.gas) {
        EQUIL.Dref = rs_param.EQUIL[0];
        EQUIL.Pref = rs_param.EQUIL[1];
        EQUIL.DOWC = rs_param.EQUIL[2];
        EQUIL.PcOW = rs_param.EQUIL[3];
    } 
    else if (rs_param.water && rs_param.oil && rs_param.gas) {
        EQUIL.Dref = rs_param.EQUIL[0];
        EQUIL.Pref = rs_param.EQUIL[1];
        EQUIL.DOWC = rs_param.EQUIL[2];
        EQUIL.PcOW = rs_param.EQUIL[3];
        EQUIL.DGOC = rs_param.EQUIL[4];
        EQUIL.PcGO = rs_param.EQUIL[5];
    }


    if (CURRENT_RANK == MASTER_PROCESS)
        cout << endl << "BLACKOIL model is selected" << endl;
}

void Bulk::InputParamCOMPS(const ParamReservoir& rs_param, OptionalFeatures& opts)
{
 
    EQUIL.Dref = rs_param.EQUIL[0];
    EQUIL.Pref = rs_param.EQUIL[1];
    EQUIL.DOWC = rs_param.EQUIL[2];
    EQUIL.PcOW = rs_param.EQUIL[3];
    EQUIL.DGOC = rs_param.EQUIL[4];
    EQUIL.PcGO = rs_param.EQUIL[5];

    // Init Zi
    for (auto& v : rs_param.ZMFVD_T.data) {
        initZi_Tab.push_back(OCPTable(v));
    }

    // Init T
    // Use RTEMP
    rsTemp = rs_param.rsTemp;
    vector<vector<OCP_DBL>> temp;
    temp.resize(2);
    // add depth
    temp[0].push_back(0);
    temp[0].push_back(1E8);
    // add temperature
    temp[1].push_back(rsTemp);
    temp[1].push_back(rsTemp);
    initT_Tab.push_back(OCPTable(temp));


    if (CURRENT_RANK == MASTER_PROCESS)
        cout << endl << "COMPOSITIONAL model is selected" << endl;
}

void Bulk::InputParamTHERMAL(const ParamReservoir& rs_param, OptionalFeatures& opts)
{
    // Init T
    rsTemp = rs_param.rsTemp;
    for (auto& v : rs_param.TEMPVD_T.data) {
        initT_Tab.push_back(OCPTable(v));
    }
    if (initT_Tab.size() == 0) {
        // Use RTEMP
        vector<vector<OCP_DBL>> temp;
        temp.resize(2);
        // add depth
        temp[0].push_back(0);
        temp[0].push_back(1E8);
        // add temperature
        temp[1].push_back(rsTemp);
        temp[1].push_back(rsTemp);
        initT_Tab.push_back(OCPTable(temp));
    }
    // ifThermal conductivity
    if (rs_param.oil) {
        thconp.push_back(rs_param.thcono);
    }
    if (rs_param.gas) {
        thconp.push_back(rs_param.thcong);
    }
    if (rs_param.water) {
        thconp.push_back(rs_param.thconw);
    }


    EQUIL.Dref = rs_param.EQUIL[0];
    EQUIL.Pref = rs_param.EQUIL[1];
    EQUIL.DOWC = rs_param.EQUIL[2];
    EQUIL.PcOW = rs_param.EQUIL[3];
    EQUIL.DGOC = rs_param.EQUIL[4];
    EQUIL.PcGO = rs_param.EQUIL[5];


    if (CURRENT_RANK == MASTER_PROCESS)
        cout << endl << "THERMAL model is selected" << endl;
}


/// Setup bulk information.
void Bulk::SetupIsoT(const Domain& domain)
{
    OCP_FUNCNAME;

    myComm  = domain.myComm;
    numproc = domain.numproc;
    myrank  = domain.myrank;
    
    // Set defaulted information
    if (vs.ntg.empty()) {
        vs.ntg.resize(vs.nb, 1);
    }
}

/// Allocate memory for fluid grid for Thermal model
void Bulk::SetupT(const Domain& domain)
{
    SetupIsoT(domain);
    vs.cType.resize(vs.nb, BulkContent::rf);
    for (OCP_USI n = 0; n < vs.nb; n++) {
        if (vs.poroInit[n] < 1E-6) {
            vs.cType[n] = BulkContent::r;
        }
    }

    // Setup Heat Loss
    hLoss.Setup(vs.nb);
}


/////////////////////////////////////////////////////////////////////
// Initial Properties
/////////////////////////////////////////////////////////////////////

void Bulk::InitPTSw(const USI& tabrow)
{
    OCP_FUNCNAME;

    initT.resize(vs.nb);

    OCP_DBL Dref = EQUIL.Dref;
    OCP_DBL Pref = EQUIL.Pref;
    OCP_DBL DOWC = EQUIL.DOWC;
    OCP_DBL PcOW = EQUIL.PcOW;
    OCP_DBL DOGC = EQUIL.DGOC;
    OCP_DBL PcGO = EQUIL.PcGO;
    OCP_DBL zRange[2];
    OCP_DBL zRangeTmp[2] = { 1E8,0 };  // min , max


    for (OCP_USI n = 0; n < vs.nb; n++) {
        OCP_DBL temp1 = vs.depth[n] - vs.dz[n] / 2;
        OCP_DBL temp2 = vs.depth[n] + vs.dz[n] / 2;
        zRangeTmp[0]  = zRangeTmp[0] < temp1 ? zRangeTmp[0] : temp1;
        zRangeTmp[1]  = zRangeTmp[1] > temp2 ? zRangeTmp[1] : temp2;
    }
    zRangeTmp[1] *= -1;

    MPI_Allreduce(&zRangeTmp, &zRange, 2, MPI_DOUBLE, MPI_MIN, myComm);
    const OCP_DBL Zmin = zRange[0];
    const OCP_DBL Zmax = -zRange[1];
    OCP_DBL tabdz = (Zmax - Zmin) / (tabrow - 1);

    vector<OCP_DBL> Ztmp(tabrow, 0);
    vector<OCP_DBL> Potmp(tabrow, 0);
    vector<OCP_DBL> Pgtmp(tabrow, 0);
    vector<OCP_DBL> Pwtmp(tabrow, 0);

    vector<OCP_DBL> tmpInitZi(vs.nc, 0);

    // cal Tab_Ztmp
    Ztmp[0] = Zmin;
    for (USI i = 1; i < tabrow; i++) {
        Ztmp[i] = Ztmp[i - 1] + tabdz;
    }

    OCP_DBL myTemp = rsTemp;

    // find the RefId
    USI beginId = 0;
    if (Dref <= Ztmp[0]) {
        beginId = 0;
    } else if (Dref >= Ztmp[tabrow - 1]) {
        beginId = tabrow - 1;
    } else {
        beginId =
            distance(Ztmp.begin(), find_if(Ztmp.begin(), Ztmp.end(),
                                           [s = Dref](auto& t) { return t > s; }));
        beginId--;
    }

    // begin calculating oil pressure:
    OCP_DBL Pbb = Pref;
    OCP_DBL gammaOtmp, gammaWtmp, gammaGtmp;
    OCP_DBL Ptmp;
    USI     mynum = 10;
    OCP_DBL mydz  = 0;
    OCP_DBL Poref, Pgref, Pwref;
    OCP_DBL Pbegin = 0;

    const auto initZi_flag = initZi_Tab.size() > 0 ? OCP_TRUE : OCP_FALSE;
    const auto initT_flag  = initT_Tab.size() > 0 ? OCP_TRUE : OCP_FALSE;
    const auto PBVD_flag   = EQUIL.PBVD.IsEmpty() ? OCP_FALSE : OCP_TRUE;

    auto PVT = PVTm.GetPVT(0);

    if (Dref < DOGC) {
        // reference pressure is gas pressure
        Pgref = Pref;
        if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
        if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Dref, 1);


        gammaGtmp = GRAVITY_FACTOR *
                    PVT->RhoPhase(Pgref, Pbb, myTemp, tmpInitZi, PhaseType::gas);
        Pbegin         = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
        Pgtmp[beginId] = Pbegin;

        // find the gas pressure
        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaGtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pgtmp[id], Pbb, myTemp,
                                                               tmpInitZi, PhaseType::gas);
            Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaGtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pgtmp[id], Pbb, myTemp,
                                                               tmpInitZi, PhaseType::gas);
            Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the oil pressure in Dref by Pgref
        Poref       = 0;
        Ptmp        = Pgref;
        mydz        = (DOGC - Dref) / mynum;
        OCP_DBL myz = Dref;

        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

            gammaGtmp = GRAVITY_FACTOR *
                        PVT->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi, PhaseType::gas);
            Ptmp += gammaGtmp * mydz;
            myz += mydz;
        }
        Ptmp -= PcGO;
        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

            gammaOtmp = GRAVITY_FACTOR *
                        PVT->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi, PhaseType::oil);
            Ptmp -= gammaOtmp * mydz;
            myz -= mydz;
        }
        Poref = Ptmp;

        // find the oil pressure in tab
        if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
        if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Dref, 1);

        gammaOtmp = GRAVITY_FACTOR *
                    PVT->RhoPhase(Poref, Pbb, myTemp, tmpInitZi, PhaseType::oil);
        Pbegin         = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * PVT->RhoPhase(Potmp[id], Pbb, myTemp,
                                                               tmpInitZi, PhaseType::oil);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * PVT->RhoPhase(Potmp[id], Pbb, myTemp,
                                                               tmpInitZi, PhaseType::oil);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the water pressure in Dref by Poref
        Pwref = 0;
        Ptmp  = Poref;
        mydz  = (DOWC - Dref) / mynum;
        myz   = Dref;

        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

            gammaOtmp = GRAVITY_FACTOR * PVT->RhoPhase(Poref, Pbb, myTemp,
                                                               tmpInitZi, PhaseType::oil);
            Ptmp += gammaOtmp * mydz;
            myz += mydz;
        }
        Ptmp -= PcOW;
        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);

            gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Ptmp, Pbb, myTemp,
                                                               tmpInitZi, PhaseType::wat);
            Ptmp -= gammaWtmp * mydz;
            myz -= mydz;
        }
        Pwref = Ptmp;

        // find the water pressure in tab
        if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);

        gammaWtmp = GRAVITY_FACTOR *
                    PVT->RhoPhase(Pwref, Pbb, myTemp, tmpInitZi, PhaseType::wat);
        Pbegin         = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pwtmp[id], Pbb, myTemp,
                                                               tmpInitZi, PhaseType::wat);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pwtmp[id], Pbb, myTemp,
                                                               tmpInitZi, PhaseType::wat);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }
    } 
    else if (Dref > DOWC) {
        OCP_DBL myz;
        // reference pressure is water pressure
        if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);

        Pwref     = Pref;
        gammaWtmp = GRAVITY_FACTOR *
                    PVT->RhoPhase(Pwref, Pbb, myTemp, tmpInitZi, PhaseType::wat);
        Pbegin         = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        // find the water pressure
        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pwtmp[id], Pbb, myTemp,
                                                               tmpInitZi, PhaseType::wat);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }
        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pwtmp[id], Pbb, myTemp,
                                                               tmpInitZi, PhaseType::wat);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the oil pressure in Dref by Pwref
        Poref = 0;
        Ptmp  = Pwref;
        mydz  = (DOWC - Dref) / mynum;
        myz   = Dref;

        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);

            gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Ptmp, Pbb, myTemp,
                                                               tmpInitZi, PhaseType::wat);
            Ptmp += gammaWtmp * mydz;
            myz += mydz;
        }
        Ptmp += PcOW;

        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

            gammaOtmp = GRAVITY_FACTOR *
                        PVT->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi, PhaseType::oil);
            Ptmp -= gammaOtmp * mydz;
            myz -= mydz;
        }
        Poref = Ptmp;

        // find the oil pressure in tab
        if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
        if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Dref, 1);

        gammaOtmp = GRAVITY_FACTOR *
                    PVT->RhoPhase(Poref, Pbb, myTemp, tmpInitZi, PhaseType::oil);
        Pbegin         = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * PVT->RhoPhase(Potmp[id], Pbb, myTemp,
                                                               tmpInitZi, PhaseType::oil);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * PVT->RhoPhase(Potmp[id], Pbb, myTemp,
                                                               tmpInitZi, PhaseType::oil);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        if (vs.gIndex >= 0) {
            // find the gas pressure in Dref by Poref
            Pgref = 0;
            Ptmp  = Poref;
            mydz  = (DOGC - Dref) / mynum;
            myz   = Dref;

            for (USI i = 0; i < mynum; i++) {
                if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
                if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
                if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

                gammaOtmp =
                    GRAVITY_FACTOR *
                    PVT->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi, PhaseType::oil);
                Ptmp += gammaOtmp * mydz;
                myz += mydz;
            }
            Ptmp += PcGO;
            for (USI i = 0; i < mynum; i++) {
                if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
                if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
                if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

                gammaGtmp =
                    GRAVITY_FACTOR *
                    PVT->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi, PhaseType::gas);
                Ptmp -= gammaGtmp * mydz;
                myz -= mydz;
            }
            Pgref = Ptmp;

            // find the gas pressure in tab
            if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Dref, 1);

            gammaGtmp      = GRAVITY_FACTOR * PVT->RhoPhase(Pgref, Pbb, myTemp,
                                                                    tmpInitZi, PhaseType::gas);
            Pbegin         = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
            Pgtmp[beginId] = Pbegin;

            for (USI id = beginId; id > 0; id--) {
                if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
                if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
                if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

                gammaGtmp =
                    GRAVITY_FACTOR * PVT->RhoPhase(Pgtmp[id], Pbb, myTemp,
                                                           tmpInitZi, PhaseType::gas);
                Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
            }
            for (USI id = beginId; id < tabrow - 1; id++) {
                if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
                if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
                if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

                gammaGtmp =
                    GRAVITY_FACTOR * PVT->RhoPhase(Pgtmp[id], Pbb, myTemp,
                                                           tmpInitZi, PhaseType::gas);
                Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
            }
        }

    } 
    else {
        OCP_DBL myz;
        // reference pressure is oil pressure
        Poref = Pref;
        if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
        if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Dref, 1);

        gammaOtmp = GRAVITY_FACTOR *
                    PVT->RhoPhase(Poref, Pbb, myTemp, tmpInitZi, PhaseType::oil);
        Pbegin         = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        // find the oil pressure
        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * PVT->RhoPhase(Potmp[id], Pbb, myTemp,
                                                               tmpInitZi, PhaseType::oil);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }
        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * PVT->RhoPhase(Potmp[id], Pbb, myTemp,
                                                               tmpInitZi, PhaseType::oil);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        if (vs.gIndex >= 0) {
            // find the gas pressure in Dref by Poref
            Pgref = 0;
            Ptmp  = Poref;
            mydz  = (DOGC - Dref) / mynum;
            myz   = Dref;

            for (USI i = 0; i < mynum; i++) {
                if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
                if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
                if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

                gammaOtmp =
                    GRAVITY_FACTOR *
                    PVT->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi, PhaseType::oil);
                Ptmp += gammaOtmp * mydz;
                myz += mydz;
            }
            Ptmp += PcGO;
            for (USI i = 0; i < mynum; i++) {
                if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
                if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
                if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

                gammaGtmp =
                    GRAVITY_FACTOR *
                    PVT->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi, PhaseType::gas);
                Ptmp -= gammaGtmp * mydz;
                myz -= mydz;
            }
            Pgref = Ptmp;

            // find the gas pressure in tab
            if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Dref, 1);

            gammaGtmp      = GRAVITY_FACTOR * PVT->RhoPhase(Pgref, Pbb, myTemp,
                                                                    tmpInitZi, PhaseType::gas);
            Pbegin         = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
            Pgtmp[beginId] = Pbegin;

            for (USI id = beginId; id > 0; id--) {
                if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
                if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
                if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

                gammaGtmp =
                    GRAVITY_FACTOR * PVT->RhoPhase(Pgtmp[id], Pbb, myTemp,
                                                           tmpInitZi, PhaseType::gas);
                Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
            }

            for (USI id = beginId; id < tabrow - 1; id++) {
                if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
                if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
                if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

                gammaGtmp =
                    GRAVITY_FACTOR * PVT->RhoPhase(Pgtmp[id], Pbb, myTemp,
                                                           tmpInitZi, PhaseType::gas);
                Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
            }
        }

        // find the water pressure in Dref by Poref
        Pwref = 0;
        Ptmp  = Poref;
        mydz  = (DOWC - Dref) / mynum;
        myz   = Dref;

        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

            gammaOtmp = GRAVITY_FACTOR *
                        PVT->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi, PhaseType::oil);
            Ptmp += gammaOtmp * mydz;
            myz += mydz;
        }
        Ptmp -= PcOW;
        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);

            gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Ptmp, Pbb, myTemp,
                                                               tmpInitZi, PhaseType::wat);
            Ptmp -= gammaWtmp * mydz;
            myz -= mydz;
        }
        Pwref = Ptmp;

        // find the water pressure in tab
        if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);

        gammaWtmp = GRAVITY_FACTOR *
                    PVT->RhoPhase(Pwref, Pbb, myTemp, tmpInitZi, PhaseType::wat);
        Pbegin         = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pwtmp[id], Pbb, myTemp,
                                                               tmpInitZi, PhaseType::wat);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pwtmp[id], Pbb, myTemp,
                                                               tmpInitZi, PhaseType::wat);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }
    }

    OCPTable DepthP(vector<vector<OCP_DBL>>{Ztmp, Potmp, Pgtmp, Pwtmp});

    if (CURRENT_RANK == MASTER_PROCESS)
        DepthP.Display();

    // calculate Pc from DepthP to calculate Sj
    std::vector<OCP_DBL> data(4, 0), cdata(4, 0);

    for (OCP_USI n = 0; n < vs.nb; n++) {
        if (initZi_flag) {
            initZi_Tab[0].Eval_All0(vs.depth[n], tmpInitZi);
            for (USI i = 0; i < vs.nc; i++) {
                vs.Ni[n * vs.nc + i] = tmpInitZi[i];
            }
        }
        if (initT_flag) {
            myTemp   = initT_Tab[0].Eval(0, vs.depth[n], 1);
            initT[n] = myTemp;
            vs.T[n]     = myTemp;
        }

        DepthP.Eval_All(0, vs.depth[n], data, cdata);
        const auto SAT  = SATm.GetSAT(n);
        auto Po   = data[1];
        auto Pg   = data[2];
        auto Pw   = data[3];
        auto Pcgo = Pg - Po;
        auto Pcow = Po - Pw;
        auto Sw   = SAT->CalSwByPcow(Pcow);
        auto Sg   = 0;
        if (vs.gIndex >= 0) {
            Sg = SAT->CalSgByPcgo(Pcgo);
        }
        if (Sw + Sg > 1) {
            // should be modified
            OCP_DBL Pcgw = Pcow + Pcgo;
            Sw           = SAT->CalSwByPcgw(Pcgw);
            Sg           = 1 - Sw;
        }

        if (1 - Sw < TINY) {
            // all water
            Po = Pw + SAT->CalPcowBySw(1.0);
        } else if (1 - Sg < TINY) {
            // all gas
            Po = Pg - SAT->CalPcgoBySg(1.0);
        } else if (1 - Sw - Sg < TINY) {
            // water and gas
            Po = Pg - SAT->CalPcgoBySg(Sg);
        }
        vs.P[n] = Po;

        if (vs.depth[n] < DOGC) {
            Pbb = Po;
        } else if (PBVD_flag) {
            Pbb = EQUIL.PBVD.Eval(0, vs.depth[n], 1);
        }
        vs.Pb[n] = Pbb;

        // cal Sw
        const OCP_DBL swco = SAT->GetSwco();
        if (fabs(SAT->CalPcowBySw(0.0 - TINY)) < TINY &&
            fabs(SAT->CalPcowBySw(1.0 + TINY) < TINY)) {
            vs.S[n * vs.np + vs.np - 1] = swco;
            continue;
        }

        Sw              = 0;
        Sg              = 0;
        USI     ncut    = 10;
        OCP_DBL avePcow = 0;

        for (USI k = 0; k < ncut; k++) {
            OCP_DBL tmpSw = 0;
            OCP_DBL tmpSg = 0;
            OCP_DBL dep   = vs.depth[n] + vs.dz[n] / ncut * (k - (ncut - 1) / 2.0);
            DepthP.Eval_All(0, dep, data, cdata);
            Po   = data[1];
            Pg   = data[2];
            Pw   = data[3];
            Pcow = Po - Pw;
            Pcgo = Pg - Po;
            avePcow += Pcow;
            tmpSw = SAT->CalSwByPcow(Pcow);
            if (vs.gIndex >= 0) {
                tmpSg = SAT->CalSgByPcgo(Pcgo);
            }
            if (tmpSw + tmpSg > 1) {
                // should be modified
                OCP_DBL Pcgw = Pcow + Pcgo;
                tmpSw        = SAT->CalSwByPcgw(Pcgw);
                tmpSg        = 1 - tmpSw;
            }
            Sw += tmpSw;
            // Sg += tmpSg;
        }
        Sw /= ncut;
        // Sg /= ncut;
        avePcow /= ncut;

        SAT->SetupScale(n, Sw, avePcow);
        vs.S[n * vs.np + vs.np - 1] = Sw;
    }
}


/////////////////////////////////////////////////////////////////////
// Basic Fluid Information
/////////////////////////////////////////////////////////////////////

OCP_DBL Bulk::CalFPR(OCP_DBL& vtmp) const
{
    OCP_FUNCNAME;

    vtmp = 0;
    OCP_DBL ptmp = 0;   
    OCP_DBL tmp  = 0;

    if (vs.np == 3) {
        for (OCP_USI n = 0; n < vs.nbI; n++) {
            tmp = vs.rockVp[n] * (1 - vs.S[n * vs.np + 2]);
            ptmp += vs.P[n] * tmp;
            vtmp += tmp;
        }
    } else if (vs.np < 3) {
        for (OCP_USI n = 0; n < vs.nbI; n++) {
            tmp = vs.rockVp[n] * (vs.S[n * vs.np]);
            ptmp += vs.P[n] * tmp;
            vtmp += tmp;
        }
    } else {
        OCP_ABORT("Number of phases is out of range!");
    }

    return ptmp / vtmp;
}

OCP_DBL Bulk::CalFTR(OCP_DBL& vtmp) const
{
    OCP_FUNCNAME;

    vtmp = 0;
    OCP_DBL Ttmp = 0;
    
    for (OCP_USI n = 0; n < vs.nbI; n++) {
        Ttmp += vs.T[n] * vs.v[n];
        vtmp += vs.v[n];
    }

    return Ttmp / vtmp;
}

/////////////////////////////////////////////////////////////////////
// Important Indicator Variable and Check
/////////////////////////////////////////////////////////////////////


/// Return OCP_TRUE if no negative pressure and OCP_FALSE otherwise.
OCP_INT Bulk::CheckP() const
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < vs.nb; n++) {
        if (vs.P[n] < 0) {
            std::ostringstream PStringSci;
            PStringSci << std::scientific << vs.P[n];
            OCP_WARNING("Negative pressure: P[" + std::to_string(n) +
                        "] = " + PStringSci.str());
            cout << "P = " << vs.P[n] << endl;
            return BULK_NEGATIVE_PRESSURE;
        }
    }

    return BULK_SUCCESS;
}

OCP_INT Bulk::CheckT() const
{
    for (OCP_USI n = 0; n < vs.nb; n++) {
        if (vs.T[n] < 0) {
            std::ostringstream PStringSci;
            PStringSci << std::scientific << vs.T[n];
            OCP_WARNING("Negative pressure: T[" + std::to_string(n) +
                        "] = " + PStringSci.str());
            cout << "T = " << vs.T[n] << endl;
            return BULK_NEGATIVE_TEMPERATURE;
        }
    }

    return BULK_SUCCESS;
}

/// Return OCP_TRUE if no negative Ni and OCP_FALSE otherwise.
OCP_INT Bulk::CheckNi()
{
    OCP_FUNCNAME;

    OCP_USI len = vs.nb * vs.nc;
    for (OCP_USI n = 0; n < len; n++) {
        if (vs.Ni[n] < 0) {
            OCP_USI bId = n / vs.nc;
            if (vs.Ni[n] > -1E-3 * vs.Nt[bId] && OCP_FALSE) {
                vs.Ni[n] = 1E-8 * vs.Nt[bId];
            } else {
                USI                cId = n - bId * vs.nc;
                std::ostringstream NiStringSci;
                NiStringSci << std::scientific << vs.Ni[n];
                OCP_WARNING("Negative Ni: Ni[" + std::to_string(cId) + "] in Bulk[" +
                            std::to_string(bId) + "] = " + NiStringSci.str());

                return BULK_NEGATIVE_COMPONENTS_MOLES;
            }
        }
    }
    return BULK_SUCCESS;
}

/// Return OCP_TRUE if all Ve < Vlim and OCP_FALSE otherwise.
OCP_INT Bulk::CheckVe(const OCP_DBL& Vlim) const
{
    OCP_FUNCNAME;

    OCP_DBL dVe = 0.0;
    for (OCP_USI n = 0; n < vs.nb; n++) {
        dVe = fabs(vs.vf[n] - vs.rockVp[n]) / vs.rockVp[n];
        if (dVe > Vlim) {
            cout << "Volume error at Bulk[" << n << "] = " << setprecision(6) << dVe
                 << " is too big!" << endl;
            return BULK_OUTRANGED_VOLUME_ERROR;
        }
    }
    return BULK_SUCCESS;
}

OCP_INT Bulk::CheckCFL(const OCP_DBL& cflLim) const
{
    if (maxCFL > cflLim)
        return BULK_OUTRANGED_CFL;
    else
        return BULK_SUCCESS;
}

void Bulk::CalMaxChange()
{
    OCP_FUNCNAME;

    dPmax       = 0;
    dTmax       = 0;
    dNmax       = 0;
    dSmax       = 0;
    eVmax       = 0;
    OCP_DBL tmp = 0;
    OCP_USI id;

    for (OCP_USI n = 0; n < vs.nb; n++) {

        // dP
        tmp = fabs(vs.P[n] - vs.lP[n]);
        if (dPmax < tmp) {
            dPmax = tmp;
        }

        // dT
        tmp = fabs(vs.T[n] - vs.lT[n]);
        if (dTmax < tmp) {
            dTmax = tmp;
        }

        // dS
        for (USI j = 0; j < vs.np; j++) {
            id  = n * vs.np + j;
            tmp = fabs(vs.S[id] - vs.lS[id]);
            if (dSmax < tmp) {
                dSmax = tmp;
            }
        }

        // dN
        for (USI i = 0; i < vs.nc; i++) {
            id  = n * vs.nc + i;
            tmp = fabs(max(vs.Ni[id], vs.lNi[id]));
            if (tmp > TINY) {
                tmp = fabs(vs.Ni[id] - vs.lNi[id]) / tmp;
                if (dNmax < tmp) {
                    dNmax = tmp;
                }
            }
        }

        // Ve
        tmp = fabs(vs.vf[n] - vs.rockVp[n]) / vs.rockVp[n];
        if (eVmax < tmp) {
            eVmax = tmp;
        }
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/09/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/