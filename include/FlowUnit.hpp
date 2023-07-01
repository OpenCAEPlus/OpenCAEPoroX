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

/// designed to deal with matters related to saturation table.
/// relative permeability, capillary pressure will be calculated here.
class FlowUnit
{
public:
    /// Default constructor.
    FlowUnit()                                                        = default;
    void Allocate(const USI& np) {
        kr.resize(np, 0);
        pc.resize(np, 0);
        dKrdS.resize(np * np, 0);
        dPcdS.resize(np * np, 0);
    }
    virtual void SetupOptionalFeatures(OptionalFeatures& optFeatures) = 0;
    virtual void
    SetupScale(const OCP_USI& bId, OCP_DBL& Swin, const OCP_DBL& Pcowin) = 0;
    /// Pcow = Po - Pw
    virtual OCP_DBL GetPcowBySw(const OCP_DBL& sw)   = 0;
    virtual OCP_DBL GetSwByPcow(const OCP_DBL& pcow) = 0;
    /// Pcgo = Pg - Po
    virtual OCP_DBL GetPcgoBySg(const OCP_DBL& sg)   = 0;
    virtual OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo) = 0;
    /// Pcgw = Pg - Pw
    virtual OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) = 0;
    /// Return the value of Scm
    virtual const vector<OCP_DBL>& GetScm() const = 0;

    /// Calculate relative permeability and capillary pressure.
    virtual void CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId) = 0;

    /// Calculate derivatives of relative permeability and capillary pressure.
    virtual void CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId) = 0;

    OCP_DBL GetSwco() const { return Swco; };
    const vector<OCP_DBL>& GetKr()const { return kr; }
    const vector<OCP_DBL>& GetPc()const { return pc; }
    const vector<OCP_DBL>& GetdKrdS()const { return dKrdS; }
    const vector<OCP_DBL>& GetdPcdS()const { return dPcdS; }

protected:
    OCP_USI         bId;   ///< index of bulk being calculated
    OCP_DBL         Swco;  ///< saturaion of connate water

    vector<OCP_DBL> kr;    ///< relative permeability of phase
    vector<OCP_DBL> pc;    ///< capillary pressure
    vector<OCP_DBL> dKrdS; ///< dKr / dPc
    vector<OCP_DBL> dPcdS; ///< dKr / dS

    vector<OCP_DBL> data;  ///< container to store the values of interpolation.
    vector<OCP_DBL> cdata; ///< container to store the slopes of interpolation.
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
    SetupScale(const OCP_USI& bId, OCP_DBL& Swin, const OCP_DBL& Pcowin) override{};
    void CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId) override;
    void CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId) override;

    OCP_DBL GetPcowBySw(const OCP_DBL& sw) override { return 0; }
    OCP_DBL GetSwByPcow(const OCP_DBL& pcow) override { return 0; }
    OCP_DBL GetPcgoBySg(const OCP_DBL& sg) override { return 0; }
    OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo) override { return 0; }
    OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) override { return 0; }

    const vector<OCP_DBL>& GetScm() const override
    {
        OCP_ABORT("Not Completed in this Condition!");
    }
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
    SetupScale(const OCP_USI& bId, OCP_DBL& Swin, const OCP_DBL& Pcowin) override{};
    void CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId) override;
    void CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId) override;
    OCP_DBL GetPcowBySw(const OCP_DBL& Sw) override { return SWOF.CalPcow(Sw); }
    OCP_DBL GetSwByPcow(const OCP_DBL& Pcow) override { return SWOF.CalSw(Pcow); }

    // useless
    OCP_DBL GetPcgoBySg(const OCP_DBL& sg) override { return 0; }
    OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo) override { return 0; }
    OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) override { return 0; }

    const vector<OCP_DBL>& GetScm() const override
    {
        OCP_ABORT("Not Completed in this Condition!");
    }

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
    SetupScale(const OCP_USI& bId, OCP_DBL& Swin, const OCP_DBL& Pcowin) override{};
    void    CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId) override;
    void    CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId) override;
    OCP_DBL GetPcgoBySg(const OCP_DBL& sg) override { return SGOF.Eval(0, sg, 3); }
    OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo) override { return SGOF.Eval(3, pcgo, 0); }

    // useless
    OCP_DBL GetPcowBySw(const OCP_DBL& sw) override { return 0; }
    OCP_DBL GetSwByPcow(const OCP_DBL& pcow) override { return 0; }
    OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) override { return 0; }

    const vector<OCP_DBL>& GetScm() const override
    {
        OCP_ABORT("Not Completed in this Condition!");
    }

protected:
    OCPTable SGOF; ///< saturation table about gas and oil.
    OCP_DBL  kroMax;
};

///////////////////////////////////////////////
// FlowUnit_ODGW
///////////////////////////////////////////////

class FlowUnit_ODGW : public FlowUnit
{
public:
    OCP_DBL CalKro_Stone2(const OCP_DBL& krow,
                          const OCP_DBL& krog,
                          const OCP_DBL& krw,
                          const OCP_DBL& krg) const;

    OCP_DBL CalKro_Default(const OCP_DBL& Sg,
                           const OCP_DBL& Sw,
                           const OCP_DBL& krog,
                           const OCP_DBL& krow) const;

    const vector<OCP_DBL>& GetScm() const override { return Scm; }

protected:
    /// oil relative permeability in the presence of connate water only, used in stone2
    OCP_DBL         kroMax;
    vector<OCP_DBL> Scm; ///< critical saturation when phase becomes mobile / immobile
};

///////////////////////////////////////////////
// FlowUnit_ODGW01
///////////////////////////////////////////////

class FlowUnit_ODGW01 : public FlowUnit_ODGW
{
public:
    FlowUnit_ODGW01() = default;
    FlowUnit_ODGW01(const ParamReservoir& rs_param, const USI& i);
    void SetupOptionalFeatures(OptionalFeatures& optFeatures) override{};
    void
    SetupScale(const OCP_USI& bId, OCP_DBL& Swin, const OCP_DBL& Pcowin) override{};
    virtual void CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId) override;
    virtual void CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId) override;

    OCP_DBL CalKro_Stone2Der(OCP_DBL  krow,
                             OCP_DBL  krog,
                             OCP_DBL  krw,
                             OCP_DBL  krg,
                             OCP_DBL  dkrwdSw,
                             OCP_DBL  dkrowdSw,
                             OCP_DBL  dkrgdSg,
                             OCP_DBL  dkrogdSg,
                             OCP_DBL& out_dkrodSw,
                             OCP_DBL& out_dkrodSg) const;
    OCP_DBL CalKro_DefaultDer(const OCP_DBL& Sg,
                              const OCP_DBL& Sw,
                              const OCP_DBL& krog,
                              const OCP_DBL& krow,
                              const OCP_DBL& dkrogSg,
                              const OCP_DBL& dkrowSw,
                              OCP_DBL&       dkroSg,
                              OCP_DBL&       dkroSw) const;

    OCP_DBL GetPcowBySw(const OCP_DBL& Sw) override { return SWOF.CalPcow(Sw); }
    OCP_DBL GetSwByPcow(const OCP_DBL& Pcow) override { return SWOF.CalSw(Pcow); }

    OCP_DBL GetPcgoBySg(const OCP_DBL& sg) override { return SGOF.Eval(0, sg, 3); }
    OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo) override { return SGOF.Eval(3, pcgo, 0); }
    OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) override
    {
        return SWPCGW.Eval_Inv(1, pcgw, 0);
    }

    void Generate_SWPCWG();

protected:
    OCPTable SGOF;   ///< saturation table about gas and oil.
    OCP_SWOF SWOF;   ///< saturation table about water and oil.
    OCPTable SWPCGW; ///< auxiliary table: saturation of water vs. capillary
                     ///< pressure between water and gas.
};

///////////////////////////////////////////////
// FlowUnit_ODGW01_Miscible
///////////////////////////////////////////////

class FlowUnit_ODGW01_Miscible : public FlowUnit_ODGW01
{
public:
    FlowUnit_ODGW01_Miscible(const ParamReservoir& rs_param, const USI& i)
        : FlowUnit_ODGW01(rs_param, i)
    {
        // gas is moveable all the time
        Scm[1]  = 0;
        /*maxPcow = SWOF.GetCol(3).front();
        minPcow = SWOF.GetCol(3).back();*/
        maxPcow = SWOF.GetMaxPc();
        minPcow = SWOF.GetMinPc();
    }
    void SetupOptionalFeatures(OptionalFeatures& optFeatures) override
    {
        misTerm   = &optFeatures.miscible;
        scaleTerm = &optFeatures.scalePcow;
        scaleTerm->Setup();
    };
    void SetupScale(const OCP_USI& bId, OCP_DBL& Swin, const OCP_DBL& Pcowin) override;
    void CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId) override;
    void CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId) override;

protected:
    ScalePcow* scaleTerm;
    Miscible*  misTerm;

    OCP_DBL maxPcow;
    OCP_DBL minPcow;

    OCP_DBL Fk;     ///< The relative permeability interpolation parameter
    OCP_DBL Fp;     ///< The capillary pressure interpolation parameter
    OCP_DBL surTen; ///< Surface tension

    /*
    OCP_DBL kroMis{0};  ///< miscible oil relative permeability
    OCP_DBL krgMis{0};  ///< miscible gas relative permeability
    OCP_DBL PcogMis{0}; ///< miscible gas capillary pressure
    OCP_DBL socrMis{0}; ///< oil critical miscible saturations
    OCP_DBL sgcrMis{0}; ///< gas critical miscible saturations
    */
};

///////////////////////////////////////////////
// FlowUnit_ODGW02
///////////////////////////////////////////////

class FlowUnit_ODGW02 : public FlowUnit_ODGW
{
public:
    FlowUnit_ODGW02() = default;
    FlowUnit_ODGW02(const ParamReservoir& rs_param, const USI& i);
    void SetupOptionalFeatures(OptionalFeatures& optFeatures) override{};
    void
    SetupScale(const OCP_USI& bId, OCP_DBL& Swin, const OCP_DBL& Pcowin) override{};
    void    CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId) override;
    void    CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId) override;
    OCP_DBL CalKro_Stone2Der(OCP_DBL  krow,
                             OCP_DBL  krog,
                             OCP_DBL  krw,
                             OCP_DBL  krg,
                             OCP_DBL  dkrwdSw,
                             OCP_DBL  dkrowdSo,
                             OCP_DBL  dkrgdSg,
                             OCP_DBL  dkrogdSo,
                             OCP_DBL& out_dkrodSo) const;
    OCP_DBL CalKro_DefaultDer(const OCP_DBL& Sg,
                              const OCP_DBL& Sw,
                              const OCP_DBL& krog,
                              const OCP_DBL& krow,
                              const OCP_DBL& dkrogSo,
                              const OCP_DBL& dkrowSo,
                              OCP_DBL&       dkroSo) const;

    OCP_DBL GetPcowBySw(const OCP_DBL& sw) override { return SWFN.Eval(0, sw, 2); }
    OCP_DBL GetSwByPcow(const OCP_DBL& pcow) override
    {
        return SWFN.Eval_Inv(2, pcow, 0);
    }
    OCP_DBL GetPcgoBySg(const OCP_DBL& sg) override { return SGFN.Eval(0, sg, 2); }
    OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo) override { return SGFN.Eval(2, pcgo, 0); }
    OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) override
    {
        return SWPCGW.Eval_Inv(1, pcgw, 0);
    }
    void Generate_SWPCWG();

protected:
    OCPTable SWFN;   ///< saturation table about water.
    OCPTable SGFN;   ///< saturation table about gas.
    OCPTable SOF3;   ///< saturation table about oil.
    OCPTable SWPCGW; ///< auxiliary table: saturation of water vs. capillary
                     ///< pressure between water and gas.
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