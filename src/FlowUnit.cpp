/*! \file    FlowUnit.cpp
 *  \brief   FlowUnit class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "FlowUnit.hpp"

///////////////////////////////////////////////
// FlowUnit_W
///////////////////////////////////////////////

void FlowUnit_W::CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId)
{
    // constant value
    //kr_out[0] = 1;
    //pc_out[0] = 0;
}

void FlowUnit_W::CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId)
{
    // constant value
    //kr_out[0] = 1;
    //pc_out[0] = 0;
    //dkrdS[0]  = 0;
    //dPcjdS[0] = 0;
}

///////////////////////////////////////////////
// FlowUnit_OW
///////////////////////////////////////////////

FlowUnit_OW::FlowUnit_OW(const ParamReservoir& rs_param, const USI& i)
{
    Allocate(2);

    SWOF.Setup(rs_param.SWOF_T.data[i]);
    Swco = SWOF.GetSwco();

    data.resize(4, 0);
    cdata.resize(4, 0);
}

void FlowUnit_OW::CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId)
{
    const OCP_DBL Sw = S_in[1];

    // three phase black oil model using stone 2
    OCP_DBL krw, krow, Pcwo;
    SWOF.CalKrwKrowPcow(Sw, krw, krow, Pcwo);

    kr[0] = krow;
    kr[1] = krw;
    //pc[0] = 0;
    pc[1] = Pcwo;
}

void FlowUnit_OW::CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId)
{
    const OCP_DBL Sw = S_in[1];
    OCP_DBL krw, krow, Pcwo, dKrwdSw, dKrowdSw, dPcwodSw;
    SWOF.CalKrwKrowPcowDer(Sw, krw, krow, Pcwo, dKrwdSw, dKrowdSw, dPcwodSw);

    kr[0] = krow;
    kr[1] = krw;
    //pc[0] = 0;
    pc[1] = Pcwo;

    //dKrdS[0] = 0;
    dKrdS[1] = dKrowdSw;
    //dKrdS[2] = 0;
    dKrdS[3] = dKrwdSw;

    //dPcdS[0] = 0;
    //dPcdS[1] = 0;
    //dPcdS[2] = 0;
    dPcdS[3] = dPcwodSw;
}

///////////////////////////////////////////////
// FlowUnit_OG
///////////////////////////////////////////////

FlowUnit_OG::FlowUnit_OG(const ParamReservoir& rs_param, const USI& i)
{
    Allocate(2);

    SGOF.Setup(rs_param.SGOF_T.data[i]);

    data.resize(4, 0);
    cdata.resize(4, 0);
}

void FlowUnit_OG::CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId)
{
    const OCP_DBL Sg = S_in[1];

    OCP_DBL krg, krog, Pcgo;
    SGOF.CalKrgKrogPcgo(Sg, krg, krog, Pcgo);

    kr[0] = krog;
    kr[1] = krg;
    //pc[0] = 0;
    pc[1] = Pcgo;
}

void FlowUnit_OG::CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId)
{
    OCP_ABORT("Not Completed Now!");
}

///////////////////////////////////////////////
// FlowUnit_OGW
///////////////////////////////////////////////


void FlowUnit_OGW::SetupOptionalFeatures(OptionalFeatures& optFeatures)
{
    // for miscible
    miscible       = &optFeatures.miscible;
    // for scalePcow
    scalePcow      = &optFeatures.scalePcow;
    scalePcowIndex = scalePcow->Setup(Swco, maxPcow, minPcow);
    CalPcow        = bind(&FlowUnit_OGW::GetPcowBySw, this, std::placeholders::_1);
}

void FlowUnit_OGW::SetupScale(const OCP_USI& bId, OCP_DBL& Swin, const OCP_DBL& Pcowin)
{
    scalePcow->SetScaleVal(scalePcowIndex, bId, Swin, Pcowin, CalPcow);
}


OCP_DBL FlowUnit_OGW::CalKro_Stone2(const OCP_DBL& krow,
                                     const OCP_DBL& krog,
                                     const OCP_DBL& krw,
                                     const OCP_DBL& krg) const
{
    // krog : oil relative permeability for a system with oil, gas and connate water
    // krow : oil relative permeability for a system with oil and water only

    OCP_DBL kro =
        krocw * ((krow / krocw + krw) * (krog / krocw + krg) - (krw + krg));
    if (kro < 0) kro = 0;

    return kro;
}

OCP_DBL FlowUnit_OGW::CalKro_Default(const OCP_DBL& Sg,
                                      const OCP_DBL& Sw,
                                      const OCP_DBL& krog,
                                      const OCP_DBL& krow) const
{
    const OCP_DBL tmp = Sg + Sw - Swco;
    if (tmp <= TINY) {
        return krocw;
    }
    const OCP_DBL kro = (Sg * krog + (Sw - Swco) * krow) / tmp;
    return kro;
}

///////////////////////////////////////////////
// FlowUnit_OGW01
///////////////////////////////////////////////

FlowUnit_OGW01::FlowUnit_OGW01(const ParamReservoir& rs_param, const USI& i)
{
    Allocate(3);

    SWOF.Setup(rs_param.SWOF_T.data[i]);
    SGOF.Setup(rs_param.SGOF_T.data[i]);

    krocw = SWOF.GetKrocw();
    Swco  = SWOF.GetSwco();

    data.resize(4, 0);
    cdata.resize(4, 0);

    Generate_SWPCWG();

    Scm.resize(3, 0);
    Scm[1] = SGOF.GetSgcr();
    Scm[2] = SWOF.GetSwcr();
}

void FlowUnit_OGW01::CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId)
{
    const OCP_DBL Sg = S_in[1];
    const OCP_DBL Sw = S_in[2];

    // three phase black oil model using stone 2
    OCP_DBL krw, krow, Pcwo;
    SWOF.CalKrwKrowPcow(Sw, krw, krow, Pcwo);

    OCP_DBL krg, krog, Pcgo;
    SGOF.CalKrgKrogPcgo(Sg, krg, krog, Pcgo);

    const OCP_DBL kro = CalKro_Stone2(krow, krog, krw, krg);
    // OCP_DBL kro = CalKro_Default(Sg, Sw, krog, krow);

    kr[0] = kro;
    kr[1] = krg;
    kr[2] = krw;
    //pc[0] = 0;
    pc[1] = Pcgo;
    pc[2] = Pcwo;
}

void FlowUnit_OGW01::CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId)
{
    const OCP_DBL Sg = S_in[1];
    const OCP_DBL Sw = S_in[2];

    // three phase black oil model using stone 2
    OCP_DBL krw, krow, Pcwo, dKrwdSw, dKrowdSw, dPcwodSw;
    SWOF.CalKrwKrowPcowDer(Sw, krw, krow, Pcwo, dKrwdSw, dKrowdSw, dPcwodSw);

    OCP_DBL krg, krog, Pcgo, dKrgdSg, dKrogdSg, dPcgodSg;
    SGOF.CalKrgKrogPcgoDer(Sg, krg, krog, Pcgo, dKrgdSg, dKrogdSg, dPcgodSg);

    OCP_DBL dKrodSg{0}, dKrodSw{0}, kro{0};

    kro = CalKro_Stone2Der(krow, krog, krw, krg, dKrwdSw, dKrowdSw, dKrgdSg, dKrogdSg,
                           dKrodSw, dKrodSg);
    // if (kro < 0) {
    //    cout << S_in[0] << "   " << S_in[1] << "   " << S_in[2] << endl;
    //    kro = 0;
    //}
    // kro = CalKro_DefaultDer(Sg, Sw, krog, krow, dKrogdSg, dKrowdSw, dKrodSg,
    // dKrodSw);

    kr[0] = kro;
    kr[1] = krg;
    kr[2] = krw;
    //pc[0] = 0;
    pc[1] = Pcgo;
    pc[2] = Pcwo;

    //dKrdS[0] = 0;
    dKrdS[1] = dKrodSg;
    dKrdS[2] = dKrodSw;
    //dKrdS[3] = 0;
    dKrdS[4] = dKrgdSg;
    //dKrdS[5] = 0;
    //dKrdS[6] = 0;
    //dKrdS[7] = 0;
    dKrdS[8] = dKrwdSw;

    //dPcdS[0] = 0;
    //dPcdS[1] = 0;
    //dPcdS[2] = 0;
    //dPcdS[3] = 0;
    dPcdS[4] = dPcgodSg;
    //dPcdS[5] = 0;
    //dPcdS[6] = 0;
    //dPcdS[7] = 0;
    dPcdS[8] = dPcwodSw;
}

OCP_DBL FlowUnit_OGW01::CalKro_Stone2Der(const OCP_DBL& krow,
                                         const OCP_DBL& krog,
                                         const OCP_DBL& krw,
                                         const OCP_DBL& krg,
                                         const OCP_DBL& dkrwdSw,
                                         const OCP_DBL& dkrowdSw,
                                         const OCP_DBL& dkrgdSg,
                                         const OCP_DBL& dkrogdSg,
                                         OCP_DBL& out_dkrodSw,
                                         OCP_DBL& out_dkrodSg) const
{
    OCP_DBL kro, dkrodSw, dkrodSg;
    kro = krocw * ((krow / krocw + krw) * (krog / krocw + krg) - (krw + krg));

    dkrodSw =
        krocw * ((dkrowdSw / krocw + dkrwdSw) * (krog / krocw + krg) - (dkrwdSw));
    dkrodSg =
        krocw * ((krow / krocw + krw) * (dkrogdSg / krocw + dkrgdSg) - (dkrgdSg));

    if (kro < 0) {
        kro     = 0;
        dkrodSg = 0;
        dkrodSw = 0;
    }
    out_dkrodSw = dkrodSw;
    out_dkrodSg = dkrodSg;
    return kro;
}

OCP_DBL FlowUnit_OGW01::CalKro_DefaultDer(const OCP_DBL& Sg,
                                          const OCP_DBL& Sw,
                                          const OCP_DBL& krog,
                                          const OCP_DBL& krow,
                                          const OCP_DBL& dkrogSg,
                                          const OCP_DBL& dkrowSw,
                                          OCP_DBL&       dkroSg,
                                          OCP_DBL&       dkroSw) const
{
    const OCP_DBL tmp = Sg + Sw - Swco;
    if (tmp <= TINY) {
        dkroSg = 0;
        dkroSw = 0;
        return krocw;
    }
    const OCP_DBL kro = (Sg * krog + (Sw - Swco) * krow) / tmp;
    dkroSg      = (krog + Sg * dkrogSg - kro) / tmp;
    dkroSw      = (krow + (Sw - Swco) * dkrowSw - kro) / tmp;
    return kro;
}

void FlowUnit_OGW01::Generate_SWPCWG()
{
    if (SGOF.IsEmpty()) OCP_ABORT("SGOF is missing!");
    if (SWOF.IsEmpty()) OCP_ABORT("SWOF is missing!");

    const std::vector<OCP_DBL> Sw(SWOF.GetSw());
    std::vector<OCP_DBL>       Pcow(SWOF.GetPcow());

    for (USI i = 0; i < Sw.size(); i++) {
        OCP_DBL Pcgo = SGOF.CalPcgo(1 - Sw[i]);
        Pcow[i] += Pcgo; // Pcgw
    }

    SWPCGW.Setup(vector<vector<OCP_DBL>>{Sw, Pcow});
}

///////////////////////////////////////////////
// FlowUnit_OGW01_Miscible
///////////////////////////////////////////////

void FlowUnit_OGW01_Miscible::CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId)
{
    if (!miscible->CalFkFp(bId, Fk, Fp)) {
        FlowUnit_OGW01::CalKrPc(S_in, bId);
    } else {
        const OCP_DBL So = S_in[0];
        const OCP_DBL Sg = S_in[1];
        const OCP_DBL Sw = S_in[2];

        OCP_DBL krw, krow, Pcwo;
        SWOF.CalKrwKrowPcow(Sw, krw, krow, Pcwo);

        OCP_DBL krg, krog, Pcgo;
        SGOF.CalKrgKrogPcgo(Sg, krg, krog, Pcgo);
        Pcgo *= Fp;

        OCP_DBL kro = CalKro_Stone2(krow, krog, krw, krg);
        // OCP_DBL kro = CalKro_Default(Sg, Sw, krog, krow);
        const OCP_DBL krgt = SGOF.CalKrg(1 - Sw);
        const OCP_DBL krh  = 0.5 * (krow + krgt);

        // from CMG, see *SIGMA
        kro = Fk * kro + (1 - Fk) * krh * So / (1 - Sw);
        krg = Fk * krg + (1 - Fk) * krh * Sg / (1 - Sw);

        kr[0] = kro;
        kr[1] = krg;
        kr[2] = krw;
        //pc[0] = 0;
        pc[1] = Pcgo;
        pc[2] = Pcwo;
    }

    scalePcow->Scale(scalePcowIndex, bId, pc[2]);
}

void FlowUnit_OGW01_Miscible::CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId)
{
    if (!miscible->CalFkFp(bId, Fk, Fp)) {
        FlowUnit_OGW01::CalKrPcFIM(S_in, bId);
    } else {
        const OCP_DBL So = S_in[0];
        const OCP_DBL Sg = S_in[1];
        const OCP_DBL Sw = S_in[2];

        OCP_DBL krw, krow, Pcwo, dKrwdSw, dKrowdSw, dPcwodSw;
        SWOF.CalKrwKrowPcowDer(Sw, krw, krow, Pcwo, dKrwdSw, dKrowdSw, dPcwodSw);

        OCP_DBL krg, krog, Pcgo, dKrgdSg, dKrogdSg, dPcgodSg;
        SGOF.CalKrgKrogPcgoDer(Sg, krg, krog, Pcgo, dKrgdSg, dKrogdSg, dPcgodSg);
        Pcgo     *= Fp;
        dPcgodSg *= Fp;

        OCP_DBL dKrodSg{0}, dKrodSw{0}, kro{0};
        kro = CalKro_Stone2Der(krow, krog, krw, krg, dKrwdSw, dKrowdSw, dKrgdSg,
                               dKrogdSg, dKrodSw, dKrodSg);

        OCP_DBL       dkrgd1_Sw = 0;
        const OCP_DBL krgt = SGOF.CalKrg(1 - Sw, dkrgd1_Sw);
        const OCP_DBL krh  = 0.5 * (krow + krgt);

        // from CMG, see *SIGMA
        kro = Fk * kro + (1 - Fk) * krh * So / (1 - Sw);
        krg = Fk * krg + (1 - Fk) * krh * Sg / (1 - Sw);

        kr[0] = kro;
        kr[1] = krg;
        kr[2] = krw;
        //pc[0] = 0;
        pc[1] = Pcgo;
        pc[2] = Pcwo;

        const OCP_DBL dkrhdSw = 0.5 * (dKrowdSw - dkrgd1_Sw);
        const OCP_DBL temp    = (1 - Fk) / (1 - Sw) * (krh / (1 - Sw) + dkrhdSw);

        dKrdS[0] = (1 - Fk) * krh / (1 - Sw);
        dKrdS[1] = Fk * dKrodSg;
        dKrdS[2] = Fk * dKrodSw + temp * So;
        //dKrdS[3] = 0;
        dKrdS[4] = Fk * dKrgdSg + (1 - Fk) * krh / (1 - Sw);
        dKrdS[5] = temp * Sg;
        //dKrdS[6] = 0;
        //dKrdS[7] = 0;
        dKrdS[8] = dKrwdSw;

        //dPcdS[0] = 0;
        //dPcdS[1] = 0;
        //dPcdS[2] = 0;
        //dPcdS[3] = 0;
        dPcdS[4] = dPcgodSg;
        //dPcdS[5] = 0;
        //dPcdS[6] = 0;
        //dPcdS[7] = 0;
        dPcdS[8] = dPcwodSw;
    }

    scalePcow->Scale(scalePcowIndex, bId, pc[2], dPcdS[8]);
}

///////////////////////////////////////////////
// FlowUnit_OGW02
///////////////////////////////////////////////

FlowUnit_OGW02::FlowUnit_OGW02(const ParamReservoir& rs_param, const USI& i)
{
    Allocate(3);

    SWFN.Setup(rs_param.SWFN_T.data[i]);
    SGFN.Setup(rs_param.SGFN_T.data[i]);
    SOF3.Setup(rs_param.SOF3_T.data[i]);

    krocw = SOF3.GetCol(1).back();
    Swco   = SWFN.GetCol(0)[0];

    data.resize(3, 0);
    cdata.resize(3, 0);

    Generate_SWPCWG();
}

void FlowUnit_OGW02::CalKrPc(const OCP_DBL* S_in, const OCP_USI& bId)
{
    const OCP_DBL So = S_in[0];
    const OCP_DBL Sg = S_in[1];
    const OCP_DBL Sw = S_in[2];

    SWFN.Eval_All(0, Sw, data, cdata);
    const OCP_DBL krw  = data[1];
    const OCP_DBL Pcwo = -data[2];

    SGFN.Eval_All(0, Sg, data, cdata);
    const OCP_DBL krg  = data[1];
    const OCP_DBL Pcgo = data[2];

    SOF3.Eval_All(0, So, data, cdata);
    const OCP_DBL krow = data[1];
    const OCP_DBL krog = data[2];

    // OCP_DBL kro = CalKro_Stone2(krow, krog, krw, krg);
    const OCP_DBL kro = CalKro_Default(Sg, Sw, krog, krow);

    kr[0] = kro;
    kr[1] = krg;
    kr[2] = krw;
    //pc[0] = 0;
    pc[1] = Pcgo;
    pc[2] = Pcwo;
}

void FlowUnit_OGW02::CalKrPcFIM(const OCP_DBL* S_in, const OCP_USI& bId)
{
    const OCP_DBL So = S_in[0];
    const OCP_DBL Sg = S_in[1];
    const OCP_DBL Sw = S_in[2];

    SWFN.Eval_All(0, Sw, data, cdata);
    const OCP_DBL krw      = data[1];
    const OCP_DBL dKrwdSw  = cdata[1];
    const OCP_DBL Pcwo     = -data[2];
    const OCP_DBL dPcwodSw = -cdata[2];

    SGFN.Eval_All(0, Sg, data, cdata);
    const OCP_DBL krg      = data[1];
    const OCP_DBL dKrgdSg  = cdata[1];
    const OCP_DBL Pcgo     = data[2];
    const OCP_DBL dPcgodSg = cdata[2];

    SOF3.Eval_All(0, So, data, cdata);
    const OCP_DBL krow     = data[1];
    const OCP_DBL dKrowdSo = cdata[1];
    const OCP_DBL krog     = data[2];
    const OCP_DBL dKrogdSo = cdata[2];

    OCP_DBL dKroSo = 0;
    const OCP_DBL kro    = CalKro_Stone2Der(krow, krog, krw, krg, dKrwdSw, dKrowdSo, dKrgdSg,
                                      dKrogdSo, dKroSo);
    // OCP_DBL kro = CalKro_DefaultDer(Sg, Sw, krog, krow, dKrogdSo, dKrowdSo, dKroSo);

    kr[0] = kro;
    kr[1] = krg;
    kr[2] = krw;
    //pc[0] = 0;
    pc[1] = Pcgo;
    pc[2] = Pcwo;

    dKrdS[0] = dKroSo;
    //dKrdS[1] = 0;
    //dKrdS[2] = 0;
    //dKrdS[3] = 0;
    dKrdS[4] = dKrgdSg;
    //dKrdS[5] = 0;
    //dKrdS[6] = 0;
    //dKrdS[7] = 0;
    dKrdS[8] = dKrwdSw;

    //dPcdS[0] = 0;
    //dPcdS[1] = 0;
    //dPcdS[2] = 0;
    //dPcdS[3] = 0;
    dPcdS[4] = dPcgodSg;
    //dPcdS[5] = 0;
    //dPcdS[6] = 0;
    //dPcdS[7] = 0;
    dPcdS[8] = dPcwodSw;
}

OCP_DBL FlowUnit_OGW02::CalKro_Stone2Der(const OCP_DBL& krow,
                                         const OCP_DBL&  krog,
                                         const OCP_DBL&  krw,
                                         const OCP_DBL&  krg,
                                         const OCP_DBL&  dkrwdSw,
                                         const OCP_DBL&  dkrowdSo,
                                         const OCP_DBL&  dkrgdSg,
                                         const OCP_DBL&  dkrogdSo,
                                         OCP_DBL& out_dkrodSo) const
{
    OCP_DBL kro, dKrodSo;

    kro     = krocw * ((krow / krocw + krw) * (krog / krocw + krg) - (krw + krg));
    dKrodSo = dkrowdSo * (krog / krocw + krg) + dkrogdSo * (krow / krocw + krw);

    if (kro < 0) {
        kro     = 0;
        dKrodSo = 0;
    }
    out_dkrodSo = dKrodSo;
    return kro;
}

OCP_DBL FlowUnit_OGW02::CalKro_DefaultDer(const OCP_DBL& Sg,
                                           const OCP_DBL& Sw,
                                           const OCP_DBL& krog,
                                           const OCP_DBL& krow,
                                           const OCP_DBL& dkrogdSo,
                                           const OCP_DBL& dkrowdSo,
                                           OCP_DBL&       out_dkrodSo) const
{
    const OCP_DBL tmp = Sg + Sw - Swco;
    OCP_DBL kro = (Sg * krog + (Sw - Swco) * krow) / tmp;
    out_dkrodSo = (Sg * dkrogdSo + (Sw - Swco) * dkrowdSo) / tmp;

    if (tmp <= TINY) {
        kro         = krocw;
        out_dkrodSo = 0;
    }
    return kro;
}

void FlowUnit_OGW02::Generate_SWPCWG()
{
    if (SWFN.IsEmpty()) OCP_ABORT("SWFN is missing!");
    if (SGFN.IsEmpty()) OCP_ABORT("SGFN is missing!");

    const std::vector<OCP_DBL> Sw(SWFN.GetCol(0));
    std::vector<OCP_DBL> Pcow(SWFN.GetCol(2));
    USI                  n = Sw.size();
    for (USI i = 0; i < n; i++) {
        OCP_DBL Pcgo = SGFN.Eval(0, 1 - Sw[i], 2);
        Pcow[i] += Pcgo; // pcgw
    }

    SWPCGW.Setup(vector<vector<OCP_DBL>>{Sw, Pcow});
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/