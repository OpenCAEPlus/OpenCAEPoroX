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
#include "OCPTable.hpp"
#include "OCPSATFunc.hpp"
#include "OptionalFeatures.hpp"
#include "ParamReservoir.hpp"
#include "OCP3PhaseFlow.hpp"

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
    virtual void SetupOptionalFeatures(OptionalFeatures& optFeatures) = 0;
    virtual void
    SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin) = 0;
    virtual OCP_DBL GetPcowBySw(const OCP_DBL& sw) const = 0;
    virtual OCP_DBL GetSwByPcow(const OCP_DBL& pcow) const = 0;
    virtual OCP_DBL GetPcgoBySg(const OCP_DBL& sg) const = 0;
    virtual OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo) const = 0;
    virtual OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) const = 0;

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

    vector<OCP_DBL> Scm;    ///< critical saturation when phase becomes mobile / immobile

    OCP_DBL         Swco;   ///< saturaion of connate water

    vector<OCP_DBL> kr;     ///< relative permeability of phase
    vector<OCP_DBL> pc;     ///< capillary pressure
    vector<OCP_DBL> dKrdS;  ///< dKr / dPc
    vector<OCP_DBL> dPcdS;  ///< dKr / dS

    vector<OCP_DBL> data;   ///< container to store the values of interpolation.
    vector<OCP_DBL> cdata;  ///< container to store the slopes of interpolation.
};

///////////////////////////////////////////////
// FlowUnit_W
///////////////////////////////////////////////

class FlowUnit_W : public FlowUnit
{
public:
    FlowUnit_W() = default;
    FlowUnit_W(const ParamReservoir& rs_param, const USI& i){ 
        Allocate(1);
        kr[0]    = 1;
        pc[0]    = 0;
        dKrdS[0] = 0;
        dPcdS[0] = 0;
    };
    void SetupOptionalFeatures(OptionalFeatures& optFeatures) override{};
    void
    SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin) override{};
    void CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId) override;
    void CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId) override;
    OCP_DBL GetSwco() const override { return 0.0; }
    OCP_DBL GetPcowBySw(const OCP_DBL& sw) const override { return 0.0; }
    OCP_DBL GetSwByPcow(const OCP_DBL& pcow) const override { return 1.0; }
    OCP_DBL GetPcgoBySg(const OCP_DBL& sg) const override { return 0.0; }
    OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo) const override { return 0.0; }
    OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) const override { return 1.0; }
};

///////////////////////////////////////////////
// FlowUnit_OW
///////////////////////////////////////////////

class FlowUnit_OW : public FlowUnit
{
public:
    FlowUnit_OW() = default;
    FlowUnit_OW(const ParamReservoir& rs_param, const USI& i);
    void SetupOptionalFeatures(OptionalFeatures& optFeatures) override{};
    void
    SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin) override{};
    void CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId) override;
    void CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId) override;
    OCP_DBL GetPcowBySw(const OCP_DBL& Sw) const override { return SWOF.CalPcow(Sw); }
    OCP_DBL GetSwByPcow(const OCP_DBL& Pcow) const override { return SWOF.CalSw(Pcow); }
    OCP_DBL GetSwco() const override { return SWOF.GetSwco(); }
    // useless
    OCP_DBL GetPcgoBySg(const OCP_DBL& sg) const override { return 0; }
    OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo) const override { return 0; }
    OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) const override { return 0; }

protected:
    OCP_SWOF SWOF; ///< saturation table about water and oil.
};

///////////////////////////////////////////////
// FlowUnit_OG
///////////////////////////////////////////////

class FlowUnit_OG : public FlowUnit
{
public:
    FlowUnit_OG() = default;
    FlowUnit_OG(const ParamReservoir& rs_param, const USI& i);
    void SetupOptionalFeatures(OptionalFeatures& optFeatures) override{};
    void
    SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin) override{};
    void    CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId) override;
    void    CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId) override;
    OCP_DBL GetPcgoBySg(const OCP_DBL& Sg) const override{ return SGOF.CalPcgo(Sg);}
    OCP_DBL GetSgByPcgo(const OCP_DBL& Pcgo) const override{ return SGOF.CalSg(Pcgo); }
    OCP_DBL GetSwco() const override { return 0; }
    // useless
    OCP_DBL GetPcowBySw(const OCP_DBL& sw) const override { return 0; }
    OCP_DBL GetSwByPcow(const OCP_DBL& pcow) const override  { return 0; }
    OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) const override { return 0; }

protected:
    OCP_SGOF SGOF; ///< saturation table about gas and oil.
};

///////////////////////////////////////////////
// FlowUnit_OGW
///////////////////////////////////////////////

class FlowUnit_OGW : public FlowUnit
{
public:
    FlowUnit_OGW(const ParamReservoir& rs_param, const USI& i) {
        Allocate(3);
        PF3.Setup(rs_param, i);
    }

    void SetupOptionalFeatures(OptionalFeatures& optFeatures) override final;
    void SetupScale(const OCP_USI& bId, OCP_DBL& Swinout, const OCP_DBL& Pcowin) override final;
    void CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId) override final;
    void CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId) override final;

    OCP_DBL GetSwco() const override { return PF3.GetSwco(); }
    OCP_DBL GetPcowBySw(const OCP_DBL& Sw) const override { return PF3.CalPcowBySw(Sw); }
    OCP_DBL GetSwByPcow(const OCP_DBL& Pcow) const override { return PF3.CalSwByPcow(Pcow); }
    OCP_DBL GetPcgoBySg(const OCP_DBL& Sg) const override { return PF3.CalPcgoBySg(Sg); }
    OCP_DBL GetSgByPcgo(const OCP_DBL& Pcgo) const override { return PF3.CalSgByPcgo(Pcgo); }
    OCP_DBL GetSwByPcgw(const OCP_DBL& Pcgw)  const override { return PF3.CalSwByPcgw(Pcgw); }

protected:
    void AssinValueDer();
    void AssinValue();

protected:
    /// phase index
    enum phaseIndex {oIndex, gIndex, wIndex};
    enum p2p { oo, og, ow, go, gg, gw, wo, wg, ww };

    OCP3PhaseFlow  PF3;

    // For scaling the water-oil capillary pressure curves
    ScalePcow* scalePcow;      ///< ptr to ScalePcow modules
    USI        spMethodIndex;  ///< index of scalePcow

    // For miscible
    Miscible*  miscible;
    USI        mcMethodIndex;  ///< index of miscible curve
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