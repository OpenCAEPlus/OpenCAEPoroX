/*! \file    FlowUnit.hpp
 *  \brief   FlowUnit class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __FLOWUNIT_HEADER__
#define __FLOWUNIT_HEADER__

#include <math.h>

// OpenCAEPoroX header files
#include "ParamReservoir.hpp"
#include "OptionalModules.hpp"
#include "OCPFlow.hpp"

/// designed to deal with matters related to saturation table.
/// relative permeability, capillary pressure will be calculated here.
class FlowUnit
{
public:
    FlowUnit(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts);

    void SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin) const;
    void CalKrPc(const OCP_USI& bId, const OCP_DBL* S) const;
    void CalKrPcFIM(const OCP_USI& bId, const OCP_DBL* S) const;
    OCP_DBL GetSwco() const { return flow->GetSwco(); }
    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const { return flow->CalPcowBySw(Sw); }
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const { return flow->CalSwByPcow(Pcow); }
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const { return flow->CalPcgoBySg(Sg); }
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const { return flow->CalSgByPcgo(Pcgo); }
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw)  const { return flow->CalSwByPcgw(Pcgw); }
    const vector<OCP_DBL>& GetKr() const { return flow->GetVarSet().kr; }
    const vector<OCP_DBL>& GetPc() const { return flow->GetVarSet().Pc; }
    const vector<OCP_DBL>& GetdKrdS() const { return flow->GetVarSet().dKrdS; }
    const vector<OCP_DBL>& GetdPcdS() const { return flow->GetVarSet().dPcdS; }

protected:

    OCPFlow*        flow;
    // For scaling the water-oil capillary pressure curves
    /// Scale Pcow
    ScalePcow*      scalePcow;
    /// Index of method
    USI             spMethodIndex;

    // For miscible
    /// Miscible curve correction
    MiscibleCurve*  misCurve;
    /// Index of method
    USI             mcMethodIndex;
};


#endif /* end if __FLOWUNIT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/05/2022      Format file                          */
/*----------------------------------------------------------------------------*/