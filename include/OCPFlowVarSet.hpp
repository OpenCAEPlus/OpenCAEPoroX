/*! \file    OCPFlowVarSet.hpp
 *  \brief   OCPFlowVarSet class declaration
 *  \author  Shizhe Li
 *  \date    Oct/04/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPFLOWVARSET_HEADER__
#define __OCPFLOWVARSET_HEADER__

#include "OCPConst.hpp"
#include <vector>

using namespace std;


enum class OCPFlowType
{
    /// single phase
    SP,
    /// oil and gas
    OG,
    /// oil and water
    OW,
    /// gas and water
    GW,
    /// oil, gas and water
    OGW
};

/// flow vars suite
class OCPFlowVarSet
{
    /// oil phase is reference phase
public:
    OCPFlowVarSet() { Init0(); }
    void Init0() {
        Swco = 0;
        krocw = 0;
        krog = 0; krow = 0;
        dKrogdSo = 0; dKrogdSg = 0;
        dKrowdSo = 0; dKrowdSw = 0;
        dPcodSo = 0; dPcodSg = 0; dPcodSw = 0;
        dPcgdSo = 0; dPcgdSg = 0; dPcgdSw = 0;
        dPcwdSo = 0; dPcwdSg = 0; dPcwdSw = 0;
    }
    void Init(const OCPFlowType& fType, const USI& numPhase, const USI& numCom) {
        flowType = fType;
        np       = numPhase;
        nc       = numCom;

        switch (flowType)
        {
        case OCPFlowType::SP:
            o = 0;
            g = -1;
            w = -1;
            break;
        case OCPFlowType::OG:
            o = 0;
            g = 1;
            w = -1;
            break;
        case OCPFlowType::OW:
            o = 0;
            w = 1;
            g = -1;
            break;
        case OCPFlowType::GW:
            g = 0;
            w = 1;
            o = -1;
            break;
        case OCPFlowType::OGW:
            o = 0;
            g = 1;
            w = 2;
            break;
        default:
            OCP_ABORT("Inavailable Flow Type!");
            break;
        }
        oo = o * np + o;
        og = o * np + g;
        ow = o * np + w;
        go = g * np + o;
        gg = g * np + g;
        gw = g * np + w;
        wo = w * np + o;
        wg = w * np + g;
        ww = w * np + w;

        // vars
        Init0();
        S.resize(np);
        kr.resize(np);
        Pc.resize(np);
        dKrdS.resize(np * np);
        dPcdS.resize(np * np);
    }

public:
    /// flow type
    OCPFlowType     flowType;
    /// index of oil, gas, water
    INT             o, g, w;
    /// auxiliary index
    INT             oo, og, ow, go, gg, gw, wo, wg, ww;
    /// num of phase, components
    USI             np, nc;
    /// saturations
    vector<OCP_DBL> S;
    /// relatve permeability
    vector<OCP_DBL> kr;
    /// Capillary pressure, relative to reference phase, Pj - Pr
    vector<OCP_DBL> Pc;
    /// dKr / dS
    vector<OCP_DBL> dKrdS;
    /// dPc / dS
    vector<OCP_DBL> dPcdS;
    /// saturaion of connate water
    OCP_DBL         Swco;
    /// oil relative permeability in the presence of connate water only
    OCP_DBL         krocw;
    /// the corresponding oil relative permeability when oil, gas and connate water are present
    OCP_DBL krog;
    /// the corresponding oil relative permeability when only oil and water are present
    OCP_DBL krow;
    /// the corresponding derivatives of permeability
    OCP_DBL dKrogdSo, dKrogdSg;
    OCP_DBL dKrowdSo, dKrowdSw;

    /// the corresponding derivatives of capillary pressure
    OCP_DBL dPcodSo, dPcodSg, dPcodSw;
    OCP_DBL dPcgdSo, dPcgdSg, dPcgdSw;
    OCP_DBL dPcwdSo, dPcwdSg, dPcwdSw;
};



#endif /* end if __OCPFLOWVARSET_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/04/2023      Create file                          */
/*----------------------------------------------------------------------------*/