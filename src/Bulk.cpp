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
// Input Param and Setup
/////////////////////////////////////////////////////////////////////

/// Read parameters from rs_param data structure.
void Bulk::Setup(const ParamReservoir& rs_param)
{
    OCP_FUNCNAME;

    rsTemp = rs_param.rsTemp;

    PVTm.Setup(rs_param, vs, optMs);
    SATm.Setup(rs_param, vs, optMs);
    ROCKm.Setup(rs_param, vs.nb, optMs);
    BOUNDm.Setup(rs_param, vs);
    INITm.Setup(rs_param, PVTm.GetMixtureType());

    ACCm.Setup(rs_param, vs, &BOUNDm);


    SetupOthers();
}


/// Setup bulk information.
void Bulk::SetupOthers()
{
    // check
    if (vs.ntg.size() != vs.nb) {
        vs.ntg.clear();
        vs.ntg.resize(vs.nb, 1.0);
    }

    if (PVTm.GetMixtureType() == OCPMixtureType::THERMALK_OW) {
        vs.cType.resize(vs.nb, BulkContent::rf);
        for (OCP_USI n = 0; n < vs.nb; n++) {
            if (vs.poroInit[n] < 1E-6) {
                vs.cType[n] = BulkContent::r;
            }
        }
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
ReservoirState Bulk::CheckP() const
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < vs.nb; n++) {
        if (vs.P[n] < 0) {
            std::ostringstream PStringSci;
            PStringSci << std::scientific << vs.P[n];
            OCP_WARNING("Negative pressure: P[" + std::to_string(n) +
                        "] = " + PStringSci.str());
            return ReservoirState::bulk_negative_P;
        }
    }

    return ReservoirState::bulk_success;
}


ReservoirState Bulk::CheckP(const OCP_DBL& dt)
{
    if (dt > 5E-2 + TINY || OCP_TRUE) {
        return CheckP();
    }
    else {

        OCP_USI count = 0;

        for (OCP_USI n = 0; n < vs.nb; n++) {
            if (vs.P[n] < 0) {
                //std::ostringstream PStringSci;
                //PStringSci << std::scientific << vs.P[n];
                //OCP_WARNING("Negative pressure: P[" + std::to_string(n) +
                //    "] = " + PStringSci.str());
                vs.P[n] = vs.lP[n];

                count++;
            }
        }

        cout << "num of negative P = " << count << endl;

        return ReservoirState::bulk_success;
    }
}


ReservoirState Bulk::CheckT() const
{
    for (OCP_USI n = 0; n < vs.nb; n++) {
        if (vs.T[n] < 0) {
            std::ostringstream PStringSci;
            PStringSci << std::scientific << vs.T[n];
            OCP_WARNING("Negative pressure: T[" + std::to_string(n) +
                        "] = " + PStringSci.str());
            return ReservoirState::bulk_negative_T;
        }
    }

    return ReservoirState::bulk_success;
}

/// Return OCP_TRUE if no negative Ni and OCP_FALSE otherwise.
ReservoirState Bulk::CheckNi()
{
    OCP_FUNCNAME;

    OCP_USI len = vs.nb * vs.nc;
    for (OCP_USI n = 0; n < len; n++) {
        if (vs.Ni[n] < 0) {
            OCP_USI bId = n / vs.nc;
            if (((vs.Ni[n] > -1E-5 * vs.Nt[bId]) ||
                 (vs.lNi[n] < 1E-3 * vs.Nt[bId] && OCP_TRUE)) &&
                OCP_TRUE)
            {
                vs.Ni[n] = 1E-40 * vs.Nt[bId];
                // vs.Ni[n] = 5E-1 * vs.lNi[n];
            }
            else
            {
                USI                cId = n - bId * vs.nc;
                std::ostringstream NiStringSci;
                NiStringSci << std::scientific << vs.Ni[n];
                std::ostringstream lNiStringSci;
                lNiStringSci << std::scientific << vs.lNi[n];
                std::ostringstream NtStringSci;
                NtStringSci << std::scientific << vs.Nt[bId];
                OCP_WARNING("Rank " + to_string(CURRENT_RANK) + " : Negative Ni: Ni[" + std::to_string(cId) + "] in Bulk[" +
                    std::to_string(bId) + "] = " + NiStringSci.str()
                    + ", lNi = " + lNiStringSci.str() + ", Nt = " + NtStringSci.str() 
                    + "  " + to_string(SATm.GetSATNUM(bId)));

                return ReservoirState::bulk_negative_N;
            }
        }
    }
    return ReservoirState::bulk_success;
}


ReservoirState Bulk::CheckNi(const OCP_DBL& dt)
{
    if (OCP_TRUE) {
        return CheckNi();
    }
    else {
        OCP_USI count = 0;

        OCP_USI len = vs.nb * vs.nc;
        for (OCP_USI n = 0; n < len; n++) {
            if (vs.Ni[n] < 0) {
                OCP_USI bId = n / vs.nc;
                vs.Ni[n] = 1E-40 * vs.Nt[bId];

                count++;
            }
        }

        cout << "num of negative Ni = " << count << endl;

        return ReservoirState::bulk_success;
    }
}


/// Return OCP_TRUE if all Ve < Vlim and OCP_FALSE otherwise.
ReservoirState Bulk::CheckVe(const OCP_DBL& Vlim) const
{
    OCP_FUNCNAME;

    OCP_DBL dVe = 0.0;
    for (OCP_USI n = 0; n < vs.nb; n++) {
        dVe = fabs(vs.vf[n] - vs.rockVp[n]) / vs.rockVp[n];
        if (dVe > Vlim) {
            cout << "Volume error at Bulk[" << n << "] = " << setprecision(6) << dVe
                 << " is too big!" << endl;
            return ReservoirState::bulk_large_EV;
        }
    }
    return ReservoirState::bulk_success;
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