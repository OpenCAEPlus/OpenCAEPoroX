/*! \file    ThermalSolver.hpp 
 *  \brief   ThermalSolver类声明 
 *  \author  Shizhe Li 
 *  \date    Nov/10/2022 
 *
 *----------------------------------------------------------------------------------- 
 *  版权所有 (C) 2021--至今 OpenCAEPoroX团队。保留所有权利。 
 *  根据GNU Lesser General Public License 3.0或更高版本的条款发布。
 *----------------------------------------------------------------------------------- */

#ifndef __THERMALSOLVER_HEADER__
#define __THERMALSOLVER_HEADER__

// OpenCAEPoroX头文件
#include "ThermalMethod.hpp"

/// ThermalSolver类用于流体解决方法。
class ThermalSolver{
public:
    /*! 
     * \brief   设置方法
     * \param   rs  油藏
     * \param   ctrl    控制参数
     */
    void SetupMethod(Reservoir& rs, const OCPControl& ctrl);

    /*! 
     * \brief   初始化油藏
     * \param   rs  油藏
     */
    void InitReservoir(Reservoir& rs);

    /*! 
     * \brief   准备
     * \param   rs  油藏
     * \param   ctrl    控制参数
     */
    void Prepare(Reservoir& rs, const OCPControl& ctrl);

    /*! 
     * \brief   组装矩阵
     * \param   rs  油藏
     * \param   ctrl    控制参数
     */
    void AssembleMat(const Reservoir& rs, OCPControl& ctrl);

    /*! 
     * \brief   解决线性系统
     * \param   rs  油藏
     * \param   ctrl    控制参数
     */
    void SolveLinearSystem(Reservoir& rs, OCPControl& ctrl);

    /*! 
     * \brief   更新流体属性
     * \param   rs  油藏
     * \param   ctrl    控制参数
     * \return  更新是否成功的布尔值
     */
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);

    /*! 
     * \brief   完成牛顿-拉夫逊迭代
     * \param   rs  油藏
     * \param   ctrl    控制参数
     * \return  完成是否成功的布尔值
     */
    OCP_BOOL FinishNR(Reservoir& rs, OCPControl& ctrl);

    /*! 
     * \brief   完成当前时间步
     * \param   rs  油藏
     * \param   ctrl    控制参数
     */
    void FinishStep(Reservoir& rs, OCPControl& ctrl);

    /*! 
     * \brief   获取NRsuite
     * \return  NRsuite对象的常量引用
     */
    const OCPNRsuite& GetNRsuite() const;

protected:
    LinearSystem LSolver;       ///< 线性系统求解器
    LinearSystem auxLSolver;    ///< 辅助线性系统求解器
    T_FIM        fim;           ///< 流体性质
};

#endif /* end if __THERMALSOLVER_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/