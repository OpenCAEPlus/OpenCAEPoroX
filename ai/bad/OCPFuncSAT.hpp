```cpp
/*! \file    OCPFuncSAT.hpp
 *  \brief   饱和度函数集合
 *  \details 本文件包含了处理油气水三相饱和度相关计算的类和方法。
 *  \author  Shizhe Li
 *  \date    Jun/29/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPFUNCSAT_HEADER__
#define __OCPFUNCSAT_HEADER__

// OpenCAEPoroX header files
#include "OCPFuncTable.hpp"
#include "ParamReservoir.hpp"
#include "OCPFlowVarSet.hpp"

using namespace std;

/*!
 * \brief   SWOF表格
 * \details 计算油水两相的相对渗透率和毛细管力
 */
class OCP_SWOF : public OCPFuncTable
{
	/// 0th column: The water saturation (Sw)
	/// 1th column: The corresponding water relative permeability (Krw)
	/// 2th column: The corresponding oil relative permeability when only oil and water are present (Krow)
	/// 3th column: The corresponding water-oil capillary pressure (Pcow = Po - Pw) (bars (METRIC), psi (FIELD))
	///             Values should be level or decreasing down the column.

public:
	/// 默认构造函数
	OCP_SWOF() = default;
	/// 返回束缚水饱和度
	OCP_DBL  GetSwco() const { return table.GetCol(0)[0]; }
	/// 返回临界水饱和度
	OCP_DBL  GetSwcr() const;
	/// 返回最大油水毛管力 Pcow: Po - Pw
	OCP_DBL  GetMaxPc() const { return table.GetCol(3).front(); }
	/// 返回最小油水毛管力: Po - Pw
	OCP_DBL  GetMinPc() const { return table.GetCol(3).back(); }
	/// 返回仅当束缚水存在时的油相相对渗透率
	OCP_DBL  GetKrocw() const { return table.GetCol(2)[0]; }
	/// 返回水相饱和度
	const vector<OCP_DBL>& GetSw() const { return table.GetCol(0); }
	/// 返回油水毛管力
	const vector<OCP_DBL>& GetPcow() const { return table.GetCol(3); }
	/// 通过油水毛管力计算水相饱和度
	OCP_DBL  CalSw(const OCP_DBL& Pcow) const { return table.Eval_Inv(3, Pcow, 0); }
	/// 通过水相饱和度计算油水毛管力
	OCP_DBL  CalPcow(const OCP_DBL& Sw) const { return table.Eval(0, Sw, 3); }
	/// 通过水相饱和度计算Krw, Krow, Pcwo
	void     CalKrwKrowPcwo(const OCP_DBL& Sw, OCP_DBL& krw, OCP_DBL& krow, OCP_DBL& Pcwo) const {
		table.Eval_All(0, Sw, data);
		krw = data[1];
		krow = data[2];
		Pcwo = -data[3];
	}
	/// 通过水相饱和度计算Krw, Krow, Pcwo和相应的导数
	void     CalKrwKrowPcwoDer(const OCP_DBL& Sw,
		OCP_DBL& krw, OCP_DBL& krow, OCP_DBL& Pcwo,
		OCP_DBL& dkrwdSw, OCP_DBL& dkrowdSw, OCP_DBL& dPcwodSw) const {
		table.Eval_All(0, Sw, data, cdata);
		krw = data[1];
		krow = data[2];
		Pcwo = -data[3];
		dkrwdSw = cdata[1];
		dkrowdSw = cdata[2];
		dPcwodSw = -cdata[3];
	}
};


/*!
 * \brief   SGOF表格
 * \details 计算油气两相的相对渗透率和毛细管力
 */
class OCP_SGOF : public OCPFuncTable
{
	/// 0th column: The gas saturation (Sg)
	/// 1th column: The corresponding gas relative permeability (Krg)
	/// 2th column: The corresponding oil relative permeability when oil, gas and connate water are present. (Krog)
	/// 3th column: The correspondin oil-gas capillary pressure (Pcgo = Pg - Po) (bars (METRIC), psi (FIELD))
	///             Values should be level or increasing down the column.
public:
	/// 默认构造函数
	OCP_SGOF() = default;
	/// 返回临界气饱和度
	OCP_DBL  GetSgcr() const;
	/// 通过油气毛管力计算气相饱和度
	OCP_DBL  CalSg(const OCP_DBL& Pcgo) const { return table.Eval(3, Pcgo, 0); }
	/// 通过气相饱和度计算气相相对渗透率
	OCP_DBL  CalKrg(const OCP_DBL& Sg) const { return table.Eval(0, Sg, 1); }
	/// 通过气相饱和度计算气相相对渗透率
	OCP_DBL  CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const { return table.Eval(0, Sg, 1, dKrgdSg); }
	/// 通过气相饱和度计算油气毛管力
	OCP_DBL  CalPcgo(const OCP_DBL& Sg) const { return table.Eval(0, Sg, 3); }

	/// 通过气相饱和度计算Krw, Krow, Pcwo
	void     CalKrgKrogPcgo(const OCP_DBL& Sg, OCP_DBL& krg, OCP_DBL& krog, OCP_DBL& Pcgo) const {
		table.Eval_All(0, Sg, data);
		krg = data[1];
		krog = data[2];
		Pcgo = data[3];
	}
	/// 通过气相饱和度计算Krw, Krow, Pcwo和相应的导数
	void     CalKrgKrogPcgoDer(const OCP_DBL& Sg,
		OCP_DBL& krg, OCP_DBL& krog, OCP_DBL& Pcgo,
		OCP_DBL& dkrgdSg, OCP_DBL& dkrogdSg, OCP_DBL& dPcgodSg) const {
		table.Eval_All(0, Sg, data, cdata);
		krg = data[1];
		krog = data[2];
		Pcgo = data[3];
		dkrgdSg = cdata[1];
		dkrogdSg = cdata[2];
		dPcgodSg = cdata[3];
	}
};


/*!
 * \brief   SOF3表格
 * \details 计算油气水三相的相对渗透率
 */
class OCP_SOF3 : public OCPFuncTable
{
	/// 0th column: The oil saturation (So)
	/// 1th column: The corresponding oil relative permeability for regions where only oil and water are present (krow)
	/// 3th column: The corresponding oil relative permeability for regions where only oil, gas and connate water are
	///             present (krog)
public:
	/// 默认构造函数
	OCP_SOF3() = default;
	/// 返回仅有束缚水存在时的油相相对渗透率
	OCP_DBL  GetKrocw() const { return table.GetCol(1).back(); }

	/// 通过油相饱和度计算Krow, krog
	void     CalKrowKrog(const OCP_DBL& So, OCP_DBL& krow, OCP_DBL& krog) const {
		table.Eval_All(0, So, data);
		krow = data[1];
		krog = data[2];
	}
	/// 通过油相饱和度计算Krow, krog和相应的导数
	void     CalKrowKrogDer(const OCP_DBL& So, OCP_DBL& krow, OCP_DBL& krog,
		OCP_DBL& dKrowdSo, OCP_DBL& dKrogdSo) const {
		table.Eval_All(0, So, data, cdata);
		krow = data[1];
		krog = data[2];
		dKrowdSo = cdata[1];
		dKrogdSo = cdata[2];
	}
};


/*!
 * \brief   SGFN表格
 * \details 计算油气两相的相对渗透率核毛细管力
 */
class OCP_SGFN : public OCPFuncTable
{
	/// 0th column: The gas saturation (Sg)
	/// 1th column: The corresponding gas relative permeability (krg)
	/// 3th column: The corresponding oil-gas capillary pressure (Pcgo = Pg - Po)
public:
	/// 默认构造函数
	OCP_SGFN() = default;
	/// 通过油气毛管力计算气相饱和度
	OCP_DBL  CalSg(const OCP_DBL& Pcgo) const { return table.Eval(2, Pcgo, 0); }
	/// 通过气相饱和度计算油气毛管力
	OCP_DBL  CalPcgo(const OCP_DBL& Sg) const { return table.Eval(0, Sg, 2); }
	/// 通过气相饱和度计算气相相对渗透率和相应的导数
	OCP_DBL  CalKrg(const OCP_DBL& Sg, OCP_DBL& dKrgdSg) const { return table.Eval(0, Sg, 1, dKrgdSg); }

	/// 通过气相饱和度计算Krg, Pcgo
	void     CalKrgPcgo(const OCP_DBL& Sg, OCP_DBL& krg, OCP_DBL& Pcgo) const {
		table.Eval_All(0, Sg, data);
		krg = data[1];
		Pcgo = data[2];
	}
	/// 通过气相饱和度计算Krg, Pcgo 和相应的导数
	void     CalKrgPcgoDer(const OCP_DBL& Sg, OCP_DBL& krg, OCP_DBL& Pcgo,
		OCP_DBL& dKrgdSg, OCP_DBL& dPcgodSg) const {
		table.Eval_All(0, Sg, data, cdata);
		krg = data[1];
		Pcgo = data[2];
		dKrgdSg = cdata[1];
		dPcgodSg = cdata[2];
	}
};


/*!
 * \brief   SWFN表格
 * \details 计算油水两相的相对渗透率核毛细管力
 */
class OCP_SWFN : public OCPFuncTable
{
	/// 0th column: The water saturation (Sw)
	/// 1th column: The corresponding water relative permeability (krw)
	/// 3th column: The corresponding water-oil capillary pressure (Pcow = Po - Pw)
public:
	/// 默认构造函数
	OCP_SWFN() = default;
	/// 返回束缚水饱和度
	OCP_DBL  GetSwco() const { return table.GetCol(0)[0]; }
	/// 返回最大油水毛管力 pcow
	OCP_DBL  GetMaxPc() const { return table.GetCol(2)[0]; }
	/// 返回最小油水毛管力 pcow
	OCP_DBL  GetMinPc() const { return table.GetCol(2).back(); }
	/// 通过油水毛管力计算水相饱和度
	OCP_DBL  CalSw(const OCP_DBL& Pcow) const { return table.Eval_Inv(2, Pcow, 0); }
	/// 通过水相饱和度计算油水毛管力
	OCP_DBL  CalPcow(const OCP_DBL& Sw) const { return table.Eval(0, Sw, 2); }
	/// 返回表格中水相饱和度序列
	const vector<OCP_DBL>& GetSw() const { return table.GetCol(0); }
	/// 返回表格中油水毛管力序列
	const vector<OCP_DBL>& GetPcow() const { return table.GetCol(2); }

	/// 通过水相饱和度计算Krw, Pcwo
	void     CalKrwPcwo(const OCP_DBL& Sw, OCP_DBL& krw, OCP_DBL& Pcwo) const {
		table.Eval_All(0, Sw, data);
		krw = data[1];
		Pcwo = -data[2];
	}
	/// 通过水相饱和度计算Krw, Pcwo和相应的导数
	void     CalKrwPcwoDer(const OCP_DBL& Sw, OCP_DBL& krw, OCP_DBL& Pcwo,
		OCP_DBL& dKrwdSw, OCP_DBL& dPcwodSw) const {
		table.Eval_All(0, Sw, data, cdata);
		krw = data[1];
		Pcwo = -data[2];
		dKrwdSw = cdata[1];
		dPcwodSw = -cdata[2];
	}
};



/*!
 * \brief   Brooks-Corey模型
 * \details 使用Brooks-Corey模型计算润湿相和非润湿相的相对渗透率和毛细管力
 */
class BrooksCorey
{
public:
	/// 组装参数
	void Setup(const BrooksCoreyParam& bcp);
	/// 计算非润湿相的相对渗透率
	OCP_DBL CalKrN(const OCP_DBL& sn) const;
	/// 计算润湿相的相对渗透率
	OCP_DBL CalKrW(const OCP_DBL& sw) const;
	/// 计算非润湿相的相对渗透率和毛细管力
	void CalKrPcN(const OCP_DBL& sn, OCP_DBL& kr, OCP_DBL& pc) const;
	/// 计算润湿相的相对渗透率和毛细管力
	void CalKrPcW(const OCP_DBL& sw, OCP_DBL& kr, OCP_DBL& pc) const;
	/// 计算非润湿相的相对渗透率和毛细管力和相应的导数
	void CalKrPcDerN(const OCP_DBL& sn, OCP_DBL& kr, OCP_DBL& pc, OCP_DBL& dKrdSn, OCP_DBL& dPcdSn) const;
	/// 计算润湿相的相对渗透率和毛细管力和相应的导数
	void CalKrPcDerW(const OCP_DBL& sw, OCP_DBL& kr, OCP_DBL& pc, OCP_DBL& dKrdSw, OCP_DBL& dPcdSw) const;

protected:
	OCP_DBL sw_imm;	   ///< 不可移动的润湿相饱和度
	OCP_DBL sn_imm;	   ///< 不可移动的非润湿相饱和度
	OCP_DBL Pentry;	   ///< 进入压力
	OCP_DBL Pcmax;	   ///< 最大毛细管力
	OCP_DBL Cw_kr;	   ///< 对于润湿相的相对渗透率计算的指数
	OCP_DBL Cn_kr;	   ///< 对于非润湿相的相对渗透率计算的指数
	OCP_DBL C_pc;      ///< 毛管力计算指数
};



/*!
 * \brief   油相相对渗透率计算方法
 * \details 油气水三相存在时油相的相对渗透率计算方法
 */
class OCP3POilPerMethod
{
public:
	/// 计算油相相对渗透率
	virtual void CalOilPer(OCPFlowVarSet& vs) = 0;
	/// 计算油相相对渗透率和相应的导数
	virtual void CalOilPerDer(OCPFlowVarSet& vs) = 0;
};


/*!
 * \brief   STONE2 公式
 * \details 油气水三相存在时油相的相对渗透率计算方法
 */
class OCP3POilPerMethod01 : public OCP3POilPerMethod
{
public:
	OCP3POilPerMethod01() = default;
	void CalOilPer(OCPFlowVarSet& vs) override;
	void CalOilPerDer(OCPFlowVarSet& vs) override;
};


/*!
 * \brief   ECLIPSE 默认公式
 * \details 油气水三相存在时油相的相对渗透率计算方法
 */
class OCP3POILPerMethod02 : public OCP3POilPerMethod
{
public:
	OCP3POILPerMethod02() = default;
	/// Calculate oil relative permeability
	void CalOilPer(OCPFlowVarSet& vs) override;
	/// Calculate oil relative permeability and derivatives
	void CalOilPerDer(OCPFlowVarSet& vs) override;
};


/*!
 * \brief   油相相对渗透率计算模块
 * \details 油气水三相存在时油相的相对渗透率计算方法
 */
class OCP3POilPerCalculation
{
public:
	/// Setup calculations
	void Setup(const USI& i) {
		if (i == 1) pM = new OCP3POilPerMethod01();
		else        pM = new OCP3POILPerMethod02();
	}
	/// Calculate oil relative permeability
	void CalOilPer(OCPFlowVarSet& vs) { pM->CalOilPer(vs); }
	/// Calculate oil relative permeability and derivatives
	void CalOilPerDer(OCPFlowVarSet& vs) { pM->CalOilPerDer(vs); }

protected:
	OCP3POilPerMethod* pM;  ///< 油相渗透率计算方法集
};



#endif // __OCPFUNCSAT_HEADER__



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jun/29/2023      Create file                          */
/*----------------------------------------------------------------------------*/

