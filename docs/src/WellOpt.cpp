/*! \file    WellOpt.cpp
 *  \brief   WellOpt class declaration
 *  \author  Shizhe Li
 *  \date    Nov/22/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// OpenCAEPoroX header files
#include "WellOpt.hpp"

WellOpt::WellOpt(const WellOptParam& optParam)
{
    if (optParam.type == "INJ") {
        type = WellType::injector;
    } else if (optParam.type == "PROD") {
        type = WellType::productor;
    } else {
        OCP_ABORT("Wrong well type!");
    }

    if (type == WellType::injector) {
        injFluidName = optParam.fluidType;
        if (injFluidName == "WAT" || injFluidName == "WATER") {
            injFluidName = "WAT";
        }
    }

    if (optParam.state == "OPEN") {
        state = WellState::open;
    } else if (optParam.state == "CLOSE") {
        state = WellState::close;
    } else {
        OCP_ABORT("Wrong state type!");
    }

    if (optParam.mode == "RATE") {
        mode = WellOptMode::irate;
    } else if (optParam.mode == "ORAT") {
        mode = WellOptMode::orate;
    } else if (optParam.mode == "GRAT") {
        mode = WellOptMode::grate;
    } else if (optParam.mode == "WRAT") {
        mode = WellOptMode::wrate;
    } else if (optParam.mode == "LRAT") {
        mode = WellOptMode::lrate;
    } else if (optParam.mode == "BHP") {
        mode = WellOptMode::bhp;
    } else {
        OCP_ABORT("Wrong well option mode!");
    }

    initMode = mode;
    maxRate  = optParam.maxRate;
    maxBHP   = optParam.maxBHP;
    minBHP   = optParam.minBHP;
    injTemp  = optParam.injTemp;


    if (type == WellType::injector) {
        tarRate = maxRate;
        tarBHP  = maxBHP;
    }
    else if (type == WellType::productor) {
        tarRate = maxRate;
        tarBHP  = minBHP;
    }
}


OCP_BOOL WellOpt::operator!=(const WellOpt& opt) const
{
    if (this->type != opt.type) return OCP_TRUE;
    if (this->state != opt.state) return OCP_TRUE;
    if (this->mode != opt.mode) return OCP_TRUE;
    if (this->initMode != opt.initMode) return OCP_TRUE;
    if (fabs(this->maxRate - opt.maxRate) > TINY) return OCP_TRUE;
    if (fabs(this->maxBHP - opt.maxBHP) > TINY) return OCP_TRUE;
    if (fabs(this->minBHP - opt.minBHP) > TINY) return OCP_TRUE;
    if (fabs(this->tarRate - opt.tarRate) > TINY) return OCP_TRUE;
    if (fabs(this->tarBHP - opt.tarBHP) > TINY) return OCP_TRUE;
    for (USI i = 0; i < injZi.size(); i++) {
        if (fabs(injZi[i] - opt.injZi[i]) > TINY) return OCP_TRUE;
    }
    for (USI i = 0; i < this->prodPhaseWeight.size(); i++) {
        if (fabs(this->prodPhaseWeight[i] - opt.prodPhaseWeight[i]) > TINY)
            return OCP_TRUE;
    }
    if (this->injPhase != opt.injPhase) return OCP_TRUE;
    if (fabs(this->injTemp - opt.injTemp) > TINY) return OCP_TRUE;
    return OCP_FALSE;
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           NOV/22/2022      Create file                          */
/*----------------------------------------------------------------------------*/
