```cpp
/*! \file    OCPFlowMethod.hpp
 *  \brief   OCPFlowMethod 类的声明和相关派生类
 *  \author  Shizhe Li
 *  \date    Oct/04/2023
 *
 *  本文件包含 OCPFlowMethod 类及其多个派生类的声明，这些类用于计算相对渗透率和毛细管压力等。
 *  这些类是石油工程模拟软件中流体流动计算的核心组件。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPFLOWMETHOD_HEADER__
#define __OCPFLOWMETHOD_HEADER__

#include "OCPFlowVarSet.hpp"
#include "OCPFuncSAT.hpp"
#include <vector>

using namespace std;

/*!
 * \class OCPFlowMethod
 * \brief 抽象基类，用于定义计算相对渗透率和毛细管压力的接口
 */
class OCPFlowMethod{
public:
    /// 默认构造函数
    OCPFlowMethod() = default;

    /// 计算相对渗透率和毛细管压力
    virtual void CalKrPc(OCPFlowVarSet& vs) = 0;

    /// 计算相对渗透率和毛细管压力及其导数
    virtual void CalKrPcDer(OCPFlowVarSet& vs) = 0;

    /// 获取束缚水的饱和度
    virtual OCP_DBL GetSwco() const = 0;

    /// 获取水和油之间的最大毛细管压力 (Po - Pw)
    virtual OCP_DBL GetMaxPcow() const = 0;

    /// 获取水和油之间的最小毛细管压力 (Po - Pw)
    virtual OCP_DBL GetMinPcow() const = 0;

    /// 通过Sw计算Pcow
    virtual OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const = 0;

    /// 通过Pcow计算Sw
    virtual OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const = 0;

    /// 通过Sg计算Pcgo
    virtual OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const = 0;

    /// 通过Pcgo计算Sg
    virtual OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const = 0;

    /// 通过Pcgw计算Sw
    virtual OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const = 0;

    /// 通过Sg计算Krg和其导数
    virtual OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const = 0;
};

/// 使用 SGOF, SWOF 表格来计算油气水三相流动情形下的相渗曲线
class OCPFlowMethod_OGW01 : public OCPFlowMethod
{
public:
    /// 构造函数
    OCPFlowMethod_OGW01(const vector<vector<OCP_DBL>>& SGOFin,
        const vector<vector<OCP_DBL>>& SWOFin,
        const USI& i, OCPFlowVarSet& vs);
    /// 计算相对渗透率和毛细管压力
    void CalKrPc(OCPFlowVarSet& vs) override;
    /// 计算相对渗透率和毛细管压力及其导数
    void CalKrPcDer(OCPFlowVarSet& vs) override;
    /// 获取束缚水的饱和度
    OCP_DBL GetSwco() const override { return SWOF.GetSwco(); }
    /// 获取水和油之间的最大毛细管压力 (Po - Pw)
    OCP_DBL GetMaxPcow() const override { return SWOF.GetMaxPc(); }
    /// 获取水和油之间的最小毛细管压力 (Po - Pw)
    OCP_DBL GetMinPcow() const override { return SWOF.GetMinPc(); }
    /// 通过Sw计算Pcow
    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { return SWOF.CalPcow(Sw); }
    /// 通过Pcow计算Sw
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const override { return SWOF.CalSw(Pcow); }
    /// 通过Sg计算Pcgo
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const override { return SGOF.CalPcgo(Sg); }
    /// 通过Pcgo计算Sg
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const override { return SGOF.CalSg(Pcgo); }
    /// 通过Pcgw计算Sw
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const override { return SWPCGW.Eval_Inv(1, Pcgw, 0); }
    /// 通过Sg计算Krg和其导数
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const override { return SGOF.CalKrg(Sg, dKrgdSg); }

protected:
    /// 生成水饱和度与气水毛管力的关系表格
    void Generate_SWPCWG();

protected:
    OCP_SGOF                SGOF;   ///< SGOF表格
    OCP_SWOF                SWOF;   ///< SWOF表格
    OCPTable                SWPCGW; ///< SWPCGW表格，水饱和度与气水毛管力变化关系

protected:
    OCP3POilPerCalculation  opC;    ///< 油气水三相流动情形下油相相对渗透率计算方法
};


/////////////////////////////////////////////////////
// OCPFlowMethod_OGW02
/////////////////////////////////////////////////////


/// 使用 SOF3, SGFN, SWFN 表格来计算油气水三相流动情形下的相渗曲线
class OCPFlowMethod_OGW02 : public OCPFlowMethod
{
public:
    /// 构造函数
    OCPFlowMethod_OGW02(const vector<vector<OCP_DBL>>& SOF3in,
        const vector<vector<OCP_DBL>>& SGFNin,
        const vector<vector<OCP_DBL>>& SWFNin,
        const USI& i, OCPFlowVarSet& vs);
    /// 计算相对渗透率和毛细管压力
    void CalKrPc(OCPFlowVarSet& vs) override;
    /// 计算相对渗透率和毛细管压力及其导数
    void CalKrPcDer(OCPFlowVarSet& vs) override;
    /// 获取束缚水的饱和度
    OCP_DBL GetSwco() const override { return SWFN.GetSwco(); }
    /// 获取水和油之间的最大毛细管压力 (Po - Pw)
    OCP_DBL GetMaxPcow() const override { return SWFN.GetMaxPc(); }
    /// 获取水和油之间的最小毛细管压力 (Po - Pw)
    OCP_DBL GetMinPcow() const override { return SWFN.GetMinPc(); }
    /// 通过Sw计算Pcow
    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { return SWFN.CalPcow(Sw); }
    /// 通过Pcow计算Sw
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const override { return SWFN.CalSw(Pcow); }
    /// 通过Sg计算Pcgo
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const override { return SGFN.CalPcgo(Sg); }
    /// 通过Pcgo计算Sg
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const override { return SGFN.CalSg(Pcgo); }
    /// 通过Pcgw计算Sw
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const override { return SWPCGW.Eval_Inv(1, Pcgw, 0); }
    /// 通过Sg计算Krg和其导数
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const override { return SGFN.CalKrg(Sg, dKrgdSg); }

protected:
    void Generate_SWPCWG();

protected:
    OCP_SOF3            SOF3;        ///< SOF3表格
    OCP_SGFN            SGFN;        ///< SGFN表格
    OCP_SWFN            SWFN;        ///< SWFN表格
    OCPTable            SWPCGW;      ///< SWPCGW表格，水饱和度与气水毛管力变化关系

    
protected:
    OCP3POilPerCalculation  opC;     ///< 油气水三相流动情形下油相相对渗透率计算方法
};


/////////////////////////////////////////////////////
// OCPFlowMethod_OW01
/////////////////////////////////////////////////////


/// 使用 SWOF 表格来计算油水两相流动情形下的相渗曲线
class OCPFlowMethod_OW01 : public OCPFlowMethod
{
public:
    /// 构造函数
    OCPFlowMethod_OW01(const vector<vector<OCP_DBL>>& SWOFin, OCPFlowVarSet& vs);
    /// 计算相对渗透率和毛细管压力
    void CalKrPc(OCPFlowVarSet& vs) override;
    /// 计算相对渗透率和毛细管压力及其导数
    void CalKrPcDer(OCPFlowVarSet& vs) override;
    /// 获取束缚水的饱和度
    OCP_DBL GetSwco() const override { return SWOF.GetSwco(); }
    /// 获取水和油之间的最大毛细管压力 (Po - Pw)
    OCP_DBL GetMaxPcow() const override { return SWOF.GetMaxPc(); }
    /// 获取水和油之间的最小毛细管压力 (Po - Pw)
    OCP_DBL GetMinPcow() const override { return SWOF.GetMinPc(); }
    /// 通过Sw计算Pcow
    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { return SWOF.CalPcow(Sw); }
    /// 通过Pcow计算Sw
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const override { return SWOF.CalSw(Pcow); }
    /// 不可使用
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const  override { OCP_ABORT("Inavailable!"); }
    /// 不可使用
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const override { OCP_ABORT("Inavailable!"); }
    /// 不可使用
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const override { OCP_ABORT("Inavailable!"); }
    /// 不可使用
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const override { OCP_ABORT("Inavailable!"); }

protected:
    OCP_SWOF            SWOF;   ///< SWOF 表格
};


/////////////////////////////////////////////////////
// OCPFlowMethod_OG01
/////////////////////////////////////////////////////


/// 使用 SGOF 表格来计算气油两相流动情形下的相渗曲线
class OCPFlowMethod_OG01 : public OCPFlowMethod
{
public:
    /// 构造函数
    OCPFlowMethod_OG01(const vector<vector<OCP_DBL>>& SGOFin, OCPFlowVarSet& vs);
    /// 计算相对渗透率和毛细管压力
    void CalKrPc(OCPFlowVarSet& vs) override;
    /// 计算相对渗透率和毛细管压力及其导数
    void CalKrPcDer(OCPFlowVarSet& vs) override;
    /// 不可使用
    OCP_DBL GetSwco() const override { OCP_ABORT("Inavailable!"); }
    /// 不可使用
    OCP_DBL GetMaxPcow() const override { OCP_ABORT("Inavailable!"); }
    /// 不可使用
    OCP_DBL GetMinPcow() const override { OCP_ABORT("Inavailable!"); }
    /// 不可使用
    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { OCP_ABORT("Inavailable!"); }
    /// 不可使用
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const override { OCP_ABORT("Inavailable!"); }
    /// 通过Sg计算Pcgo
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const override { return SGOF.CalPcgo(Sg); }
    /// 通过Pcgo计算Sg
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const override { return SGOF.CalSg(Pcgo); }
    /// 不可使用
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const override { OCP_ABORT("Inavailable!"); }
    /// 通过Sg计算Krg和其导数
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const override { return SGOF.CalKrg(Sg, dKrgdSg); }

protected:
    OCP_SGOF            SGOF;  ///< SGOF 表格
};


/////////////////////////////////////////////////////
// OCPFlowMethod_GW01
/////////////////////////////////////////////////////


/// 使用 Brooks-Corey 公式来计算气水两相流动情形下的相渗曲线
class OCPFlowMethod_GW01 : public OCPFlowMethod
{
public:
    /// 构造函数
    OCPFlowMethod_GW01(const BrooksCoreyParam& bcp, OCPFlowVarSet& vs) {
        vs.Init(OCPFlowType::GW, 2, 2);
        bc.Setup(bcp);
    }
    /// 计算相对渗透率和毛细管压力
    void CalKrPc(OCPFlowVarSet& vs) override;
    /// 计算相对渗透率和毛细管压力及其导数
    void CalKrPcDer(OCPFlowVarSet& vs) override;
    /// 获取束缚水的饱和度
    OCP_DBL GetSwco() const override { OCP_ABORT("Inavailable!"); }
    /// 获取水和油之间的最大毛细管压力 (Po - Pw)
    OCP_DBL GetMaxPcow() const override { OCP_ABORT("Inavailable!"); }
    /// 获取水和油之间的最小毛细管压力 (Po - Pw)
    OCP_DBL GetMinPcow() const override { OCP_ABORT("Inavailable!"); }
    /// 不可使用
    OCP_DBL CalPcowBySw(const OCP_DBL& Sw) const override { OCP_ABORT("Inavailable!"); }
    /// 不可使用
    OCP_DBL CalSwByPcow(const OCP_DBL& Pcow) const override { OCP_ABORT("Inavailable!"); }
    /// 不可使用
    OCP_DBL CalPcgoBySg(const OCP_DBL& Sg) const override { OCP_ABORT("Inavailable!"); }
    /// 不可使用
    OCP_DBL CalSgByPcgo(const OCP_DBL& Pcgo) const override { OCP_ABORT("Inavailable!"); }
    /// 不可使用
    OCP_DBL CalSwByPcgw(const OCP_DBL& Pcgw) const override { OCP_ABORT("Inavailable!"); }
    /// 不可使用
    OCP_DBL CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const override { OCP_ABORT("Inavailable!"); }

protected:
    BrooksCorey   bc; ///< Brooks-Corey 公式计算相渗性质
};

#endif /* end if __OCPFLOWMETHOD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/04/2023      Create file                          */
/*----------------------------------------------------------------------------*/
```

