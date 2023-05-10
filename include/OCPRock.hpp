/*! \file    OCPRock.hpp
 *  \brief   OCPRock class declaration
 *  \author  Shizhe Li
 *  \date    Nov/15/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __ROCK_HEADER__
#define __ROCK_HEADER__

#include <math.h>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "OCPTable.hpp"
#include "ParamReservoir.hpp"

class OCPRock
{
public:
    OCPRock() = default;

    // Calculate porosity and d porosity / d P for isothermal model
    virtual void CalPoro(const OCP_DBL& P, const OCP_DBL& poroInit) = 0;

    // Calculate porosity and d porosity / d P for ifThermal model
    virtual void CalPoroT(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& poroInit) = 0;

    // Calculate enthalpy of unit volume of rock
    virtual void CalRockHT(const OCP_DBL& T) = 0;

    OCP_DBL GetPoro() const { return poro; }
    OCP_DBL GetdPorodP() const { return dPorodP; }
    OCP_DBL GetdPorodT() const { return dPorodT; }
    OCP_DBL Get_Poro() const { return _poro; }
    OCP_DBL Get_dPorodP() const { return _dPorodP; }
    OCP_DBL Get_dPorodT() const { return _dPorodT; }
    OCP_DBL GetHr() const { return Hr; }
    OCP_DBL GetdHrdT() const { return dHrdT; }

protected:
    OCP_DBL   poro;      ///< porosity of rock
    OCP_DBL   dPorodP;   ///< d poro / d P
    OCP_DBL   dPorodT;   ///< d poro / d T
    OCP_DBL   _poro;     ///< ratio of the non-porous part of the rock
    OCP_DBL   _dPorodP;  ///< d _poro / d P
    OCP_DBL   _dPorodT;  ///< d _poro / d T

    OCP_DBL   Hr;        ///< enthalpy of unit volume of rock
    OCP_DBL   dHrdT;     ///< d Hr / dT
};

class OCPRockIsoT : public OCPRock
{
public:
    OCPRockIsoT() = default;
    void CalPoroT(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& poroInit) override
    {
        OCP_ABORT("Not Used!");
    };
    void CalRockHT(const OCP_DBL& T) override
    {
        OCP_ABORT("Not Used!");
    };
};

class OCPRock_Linear : public OCPRockIsoT
{
    // poro = poroInit * (1 + phi)
    // poro = poroInit * (1 + cp1 * (P - Pref) + cp2 / 2 * (P - Pref) * (P - Pref))
public:
    OCPRock_Linear() = default;
    OCPRock_Linear(const RockParam& param)
    {
        Pref = param.Pref;
        cp1  = param.cp1;
        cp2  = param.cp2;
    };
    void CalPoro(const OCP_DBL& P, const OCP_DBL& poroInit) override;

protected:
    OCP_DBL Pref;
    OCP_DBL cp1;
    OCP_DBL cp2;
};

class OCPRockT : public OCPRock
{
public:
    OCPRockT() = default;
    void Assign(const RockParam& param)
    {
        Pref      = param.Pref;
        Tref      = param.Tref;
        cp        = param.cp1;
        ct        = param.ct;
        cpt       = param.cpt;
        ConstRock = param.ConstRock;
        hcp1      = param.HCP1;
        hcp2      = param.HCP2;
    }
    void CalPoro(const OCP_DBL& P, const OCP_DBL& poroInit) override
    {
        OCP_ABORT("Not Used!");
    };
    void CalRockHT(const OCP_DBL& T) override;

protected:
    OCP_DBL  Pref;      ///< Reference pressure at initial porosity.
    OCP_DBL  Tref;      ///< Reference temperature at initial porosity.
    OCP_DBL  cp;        ///< Compressibility factor of rock in reservoir.
    OCP_DBL  ct;        ///< Expansion factor of rock in reservoir, ifThermal only
    OCP_DBL  cpt;       ///< cross items, ifThermal only
    OCP_BOOL ConstRock; ///< if true, rock volume remains const, else, bulk volume
                        ///< remains const
    OCP_DBL hcp1;       ///< coefficients of the rock enthalpy formula, Btu/ft^3 - F
    OCP_DBL hcp2;       ///< coefficients of the rock enthalpy formula, Btu/ft^3 - F
};

class OCPRockT_Linear : public OCPRockT
{
    // poro = poroInit * (1 + phi)
    // poro = poroInit * (1 + cp*(P-Pref) - ct*(T-Tref) + cpt*(P-Pref)*(T-Tref))
public:
    OCPRockT_Linear() = default;
    OCPRockT_Linear(const RockParam& param) { Assign(param); };
    void CalPoroT(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& poroInit) override;
};

class OCPRockT_Exp : public OCPRockT
{
    // poro = poroInit * (1 + exp(phi))
    // poro = poroInit * exp(cp*(P-Pref) - ct*(T-Tref) + cpt*(P-Pref)*(T-Tref))
public:
    OCPRockT_Exp() = default;
    OCPRockT_Exp(const RockParam& param) { Assign(param); };
    void CalPoroT(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& poroInit) override;
};

#endif /* end if __ROCK_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/15/2022      Create file                          */
/*----------------------------------------------------------------------------*/
