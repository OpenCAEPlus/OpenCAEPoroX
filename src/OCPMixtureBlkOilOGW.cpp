///*! \file    OCPMixtureBlkOilOGW.cpp
// *  \brief   OCPMixtureBlkOilOGW class declaration
// *  \author  Shizhe Li
// *  \date    Jul/19/2023
// *
// *-----------------------------------------------------------------------------------
// *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
// *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
// *-----------------------------------------------------------------------------------
// */
//
//#include "OCPMixtureBlkOilOGW.hpp"
//
///////////////////////////////////////////////////////
//// OCPMixtureBlkOilOGWMethod
///////////////////////////////////////////////////////
//
//
//OCP_DBL OCPMixtureBlkOilOGWMethod01::CalRho(const OCP_DBL& P, const OCP_DBL& Pb, const USI& tarPhase)
//{
//    if (tarPhase == OIL) {
//        return PVCO.CalRhoO(P, Pb);
//    }
//    else if (tarPhase == GAS) {
//        return PVDG.CalRhoG(P);
//    }
//    else if (tarPhase == WATER) {
//        return PVTW.CalRhoW(P);
//    }
//    else {
//        OCP_ABORT("WRONG tarPhase!");
//    }
//}
//
//
//OCP_DBL OCPMixtureBlkOilOGWMethod01::CalXi(const OCP_DBL& P, const USI& tarPhase)
//{
//    if (tarPhase == OIL) {
//        return PVCO.CalXiO(P);
//    }
//    if (tarPhase == GAS) {
//        return PVDG.CalXiG(P);
//    }
//    else if (tarPhase == WATER) {
//        return PVTW.CalXiW(P);
//    }
//    else {
//        OCP_ABORT("Wrong tarPhase!");
//    }
//}
//
//
//OCP_DBL OCPMixtureBlkOilOGWMethod01::CalXi(const OCP_DBL& P, const USI& tarPhase)
//{
//
//}
//
///*----------------------------------------------------------------------------*/
///*  Brief Change History of This File                                         */
///*----------------------------------------------------------------------------*/
///*  Author              Date             Actions                              */
///*----------------------------------------------------------------------------*/
///*  Shizhe Li           Jul/19/2023      Create file                          */
///*----------------------------------------------------------------------------*/