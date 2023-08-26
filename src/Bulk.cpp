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
void Bulk::InputParam(const ParamReservoir& rs_param, OptionalFeatures& opts)
{
    OCP_FUNCNAME;

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
    SATm.Setup(rs_param, vs.nb, PVTm.GetMixtureType(), opts);
    ROCKm.Setup(rs_param, vs.nb, opts);
    BCm.Setup(rs_param, vs.nb);
    INITm.Setup(rs_param, PVTm.GetMixtureType());
}

void Bulk::InputParamBLKOIL(const ParamReservoir& rs_param, OptionalFeatures& opts)
{
    if (CURRENT_RANK == MASTER_PROCESS)
        cout << endl << "BLACKOIL model is selected" << endl;
}

void Bulk::InputParamCOMPS(const ParamReservoir& rs_param, OptionalFeatures& opts)
{
    rsTemp = rs_param.rsTemp;

    if (CURRENT_RANK == MASTER_PROCESS)
        cout << endl << "COMPOSITIONAL model is selected" << endl;
}

void Bulk::InputParamTHERMAL(const ParamReservoir& rs_param, OptionalFeatures& opts)
{
    // Init T
    rsTemp = rs_param.rsTemp;
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


    if (CURRENT_RANK == MASTER_PROCESS)
        cout << endl << "THERMAL model is selected" << endl;
}


/// Setup bulk information.
void Bulk::SetupIsoT(const Domain& domain)
{
    OCP_FUNCNAME;
    
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