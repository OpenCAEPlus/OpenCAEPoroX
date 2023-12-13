/*! \file    IsothermalSolver.hpp
 *  \brief   IsothermalSolver类的声明
 *  \author  Shizhe Li
 *  \date    Oct/21/2021
 *
 *  本文件包含IsothermalSolver类的声明，用于流体解决方法。
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __ISOTHERMALSOLVER_HEADER__
#define __ISOTHERMALSOLVER_HEADER__

// OpenCAEPoroX头文件
#include "IsoThermalMethod.hpp"

/*!
 *  \class   IsothermalSolver
 *  \brief   用于流体解决方法的IsothermalSolver类。
 */
class IsothermalSolver
{
public:
    /*!
     *  \brief  设置流体求解器。
     *  \param  rs   油藏对象引用。
     *  \param  ctrl 控制参数对象常量引用。
     */
    void SetupMethod(Reservoir& rs, const OCPControl& ctrl);

    /*!
     *  \brief  初始化油藏并为某些方法准备变量。
     *  \param  rs   油藏对象引用。
     */
    void InitReservoir(Reservoir& rs);

    /*!
     *  \brief  准备组装矩阵。
     *  \param  rs   油藏对象引用。
     *  \param  ctrl 控制参数对象引用。
     */
    void Prepare(Reservoir& rs, OCPControl& ctrl);

    /*!
     *  \brief  组装矩阵。
     *  \param  rs   油藏对象常量引用。
     *  \param  ctrl 控制参数对象常量引用。
     */
    void AssembleMat(const Reservoir& rs, OCPControl& ctrl);

    /*!
     *  \brief  在单一问题中求解线性系统。
     *  \param  rs   油藏对象引用。
     *  \param  ctrl 控制参数对象引用。
     */
    void SolveLinearSystem(Reservoir& rs, OCPControl& ctrl);

    /*!
     *  \brief  更新流体属性。
     *  \param  rs   油藏对象引用。
     *  \param  ctrl 控制参数对象引用。
     *  \return       更新成功返回true，否则返回false。
     */
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);

    /*!
     *  \brief  完成牛顿-拉夫森迭代。
     *  \param  rs   油藏对象引用。
     *  \param  ctrl 控制参数对象引用。
     *  \return       完成迭代返回true，否则返回false。
     */
    OCP_BOOL FinishNR(Reservoir& rs, OCPControl& ctrl);

    /*!
     *  \brief  完成当前时间步骤。
     *  \param  rs   油藏对象引用。
     *  \param  ctrl 控制参数对象引用。
     */
    void FinishStep(Reservoir& rs, OCPControl& ctrl);

    /*!
     *  \brief  获取牛顿-拉夫森套件。
     *  \return       返回OCPNRsuite对象的常量引用。
     */
    const OCPNRsuite& GetNRsuite() const;

private:
    OCPNLMethod  curMethod{ OCPNLMethod::none }; //!< 当前方法，默认为none。
    OCPNLMethod  mainMethod{ OCPNLMethod::none }; //!< 主方法，默认为none。
    LinearSystem LSolver;                         //!< 线性系统求解器。
    IsoT_IMPEC   impec;                           //!< IMPEC方法对象。
    IsoT_FIM     fim;                             //!< FIM方法对象。
    IsoT_AIMc    aimc;                            //!< AIMc方法对象。
    IsoT_FIMddm  fim_ddm;                         //!< FIMddm方法对象。
};

#endif /* end if __ISOTHERMALSOLVER_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/21/2021      Create file                          */
/*  Chensong Zhang      Jan/16/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/
