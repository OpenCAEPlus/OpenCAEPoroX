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
#include "OCPConst.hpp"
#include "OCPFuncSAT.hpp"
#include "OptionalModules.hpp"
#include "ParamReservoir.hpp"
#include "OCPFlow.hpp"

/// designed to deal with matters related to saturation table.
/// relative permeability, capillary pressure will be calculated here.
class FlowUnit
{
public:
    /// Default constructor.
    FlowUnit() = default;
    void Allocate(const USI& np) {
        kr.resize(np, 0);
        pc.resize(np, 0);
        dKrdS.resize(np * np, 0);
        dPcdS.resize(np * np, 0);
    }
    virtual void SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin) = 0;
    virtual OCP_DBL CalPcowBySw(const OCP_DBL& sw) const = 0;
    virtual OCP_DBL CalSwByPcow(const OCP_DBL& pcow) const = 0;
    virtual OCP_DBL CalPcgoBySg(const OCP_DBL& sg) const = 0;
    virtual OCP_DBL CalSgByPcgo(const OCP_DBL& pcgo) const = 0;
    virtual OCP_DBL CalSwByPcgw(const OCP_DBL& pcgw) const = 0;

    /// Calculate relative permeability and capillary pressure.
    virtual void CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId) = 0;

    /// Calculate derivatives of relative permeability and capillary pressure.
    virtual void CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId) = 0;

    virtual OCP_DBL GetSwco() const = 0;
    const vector<OCP_DBL>& GetKr()const { return kr; }
    const vector<OCP_DBL>& GetPc()const { return pc; }
    const vector<OCP_DBL>& GetdKrdS()const { return dKrdS; }
    const vector<OCP_DBL>& GetdPcdS()const { return dPcdS; }

protected:

    OCP_USI         bulkId; ///< index of work bulk

    vector<OCP_DBL> kr;     ///< relative permeability of phase
    vector<OCP_DBL> pc;     ///< capillary pressure
    vector<OCP_DBL> dKrdS;  ///< dKr / dPc
    vector<OCP_DBL> dPcdS;  ///< dKr / dS
};

///////////////////////////////////////////////
// FlowUnit_SP
///////////////////////////////////////////////

/// Single phase
class FlowUnit_SP : public FlowUnit
{
public:
    FlowUnit_SP(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts){
        Allocate(1);
        kr[0]    = 1;
        pc[0]    = 0;
        dKrdS[0] = 0;
        dPcdS[0] = 0;
    };
    void
    SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin) override{};
    void CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId) override {};
    void CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId) override {};
    OCP_DBL GetSwco() const override { return 0.0; }
    OCP_DBL CalPcowBySw(const OCP_DBL& sw) const override { return 0.0; }
    OCP_DBL CalSwByPcow(const OCP_DBL& pcow) const override { return 1.0; }
    OCP_DBL CalPcgoBySg(const OCP_DBL& sg) const override { return 0.0; }
    OCP_DBL CalSgByPcgo(const OCP_DBL& pcgo) const override { return 0.0; }
    OCP_DBL CalSwByPcgw(const OCP_DBL& pcgw) const override { return 1.0; }
};

///////////////////////////////////////////////
// FlowUnit_OW
///////////////////////////////////////////////

/// There exists oil and water, the reference phase is oil
class FlowUnit_OW : public FlowUnit
{
public:
    FlowUnit_OW(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts) {
        Allocate(2);
        OWF.Setup(rs_param, i);
        // for scalePcow
        scalePcow     = &opts.scalePcow;
        spMethodIndex = scalePcow->Setup(&OWF);
    }
    void SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin) override;

    void CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId) override;
    void CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId) override;

    OCP_DBL GetSwco() const override { return OWF.GetSwco(); }
    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { return OWF.CalPcowBySw(Sw); }
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const override { return OWF.CalSwByPcow(Pcow); }
    // useless
    OCP_DBL CalPcgoBySg(const OCP_DBL& sg) const override { return 0; }
    OCP_DBL CalSgByPcgo(const OCP_DBL& pcgo) const override { return 0; }
    OCP_DBL CalSwByPcgw(const OCP_DBL& pcgw) const override { return 0; }

protected:
    void AssinValueDer();
    void AssinValue();

protected:
    /// phase index
    enum phaseIndex { oIndex, wIndex };
    enum p2p { oo, ow, wo, ww };

    OCPFlowOW  OWF;

    // For scaling the water-oil capillary pressure curves
    ScalePcow* scalePcow;      ///< ptr to ScalePcow modules
    USI        spMethodIndex;  ///< index of scalePcow
};

///////////////////////////////////////////////
// FlowUnit_OG
///////////////////////////////////////////////

/// There exists oil and gas, the reference phase is oil
class FlowUnit_OG : public FlowUnit
{
public:
    FlowUnit_OG(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts) {
        Allocate(2);
        OGF.Setup(rs_param, i);
    }
    void SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin) override{};
    void    CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId) override;
    void    CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId) override;
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const override { return OGF.CalPcgoBySg(Sg); }
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const override { return OGF.CalSgByPcgo(Pcgo); }
    OCP_DBL GetSwco() const override { return 0; }
    // useless
    OCP_DBL CalPcowBySw(const OCP_DBL& sw) const override { return 0; }
    OCP_DBL CalSwByPcow(const OCP_DBL& pcow) const override  { return 0; }
    OCP_DBL CalSwByPcgw(const OCP_DBL& pcgw) const override { return 0; }

protected:
    void AssinValueDer();
    void AssinValue();

protected:
    /// phase index
    enum phaseIndex { oIndex, gIndex };
    enum p2p { oo, og, go, gg };

    OCPFlowOG  OGF;
};


///////////////////////////////////////////////
// FlowUnit_GW
///////////////////////////////////////////////

/// There exists gas and water, the reference phase is gas
class FlowUnit_GW : public FlowUnit
{
public:
    FlowUnit_GW(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts) {
        Allocate(2);
        GWF.Setup(rs_param, i);
    }   
    void    CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId) override;
    void    CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId) override;
    OCP_DBL GetSwco() const override { return 0; }
    // useless
    void SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin) override {};
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const override { return 0; }
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const override { return 0; }
    OCP_DBL CalPcowBySw(const OCP_DBL& sw) const override { return 0; }
    OCP_DBL CalSwByPcow(const OCP_DBL& pcow) const override { return 0; }
    OCP_DBL CalSwByPcgw(const OCP_DBL& pcgw) const override { return 0; }

protected:
    void AssinValueDer();
    void AssinValue();

protected:
    /// phase index
    enum phaseIndex { gIndex, wIndex };
    enum p2p { gg, gw, wg, ww };

    OCPFlowGW  GWF;
};


///////////////////////////////////////////////
// FlowUnit_OGW
///////////////////////////////////////////////

/// There exists oil, gas and water, the reference phase is oil
class FlowUnit_OGW : public FlowUnit
{
public:
    FlowUnit_OGW(const ParamReservoir& rs_param, const USI& i, OptionalModules& opts) {
        Allocate(3);
        OGWF.Setup(rs_param, i);
        // for miscible
        misCurve = &opts.misCur;
        mcMethodIndex = misCurve->Setup(&OGWF, &opts.misFac);
        // for scalePcow
        scalePcow = &opts.scalePcow;
        spMethodIndex = scalePcow->Setup(&OGWF);
    }

    void SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin) override final;
    void CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId) override final;
    void CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId) override final;

    OCP_DBL GetSwco() const override { return OGWF.GetSwco(); }
    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { return OGWF.CalPcowBySw(Sw); }
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const override { return OGWF.CalSwByPcow(Pcow); }
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const override { return OGWF.CalPcgoBySg(Sg); }
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const override { return OGWF.CalSgByPcgo(Pcgo); }
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw)  const override { return OGWF.CalSwByPcgw(Pcgw); }

protected:
    void AssinValueDer();
    void AssinValue();

protected:
    /// phase index
    enum phaseIndex {oIndex, gIndex, wIndex};
    enum p2p { oo, og, ow, go, gg, gw, wo, wg, ww };

    OCPFlowOGW  OGWF;

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