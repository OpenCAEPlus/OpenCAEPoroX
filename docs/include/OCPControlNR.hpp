/*! \file    OCPControlNR.hpp
 *  \brief   OCPControlNR类的声明，用于控制牛顿迭代过程的参数和行为。
 *  \author  Shizhe Li
 *  \date    Dec/05/2023
 *
 *  OCPControlNR类及其辅助类ControlNRParam用于管理牛顿迭代过程中的控制参数，
 *  包括迭代次数、收敛容差等，以确保数值解的稳定性和准确性。
 */

#ifndef __OCPCONTROLNR_HEADER__
#define __OCPCONTROLNR_HEADER__

// 标准头文件
#include <vector>

// OpenCAEPoroX头文件
#include "OCPConst.hpp"
#include "ParamControl.hpp"
#include "OCPNRsuite.hpp"

using namespace std;

/// \brief 牛顿迭代控制参数类
/// \note 对于求解方法的收敛至关重要
class ControlNRParam
{
    friend class ControlNR; ///< 允许ControlNR类访问私有成员

public:
    /// 默认构造函数
    ControlNRParam() = default;
    /// 从给定的源向量构造控制参数对象
    ControlNRParam(const vector<OCP_DBL>& src);

protected:
    /// \brief 在一个时间步中牛顿迭代的最大次数
    /// \details 对于算法性能至关重要
    USI     maxIter;    ///< 类型：无符号短整型，默认值未定义

    /// \brief 非线性收敛误差的最大值
    /// \details 用于判断迭代是否收敛
    OCP_DBL tol;        ///< 类型：双精度浮点数，默认值未定义

    /// \brief 一个牛顿迭代中压力变化的最大值
    OCP_DBL dPmax;      ///< 类型：双精度浮点数，默认值未定义

    /// \brief 一个牛顿迭代中饱和度变化的最大值
    OCP_DBL dSmax;      ///< 类型：双精度浮点数，默认值未定义

    /// \brief 一个牛顿迭代中压力变化的最小值
    OCP_DBL dPmin;      ///< 类型：双精度浮点数，默认值未定义

    /// \brief 一个牛顿迭代中饱和度变化的最小值
    OCP_DBL dSmin;      ///< 类型：双精度浮点数，默认值未定义

    /// \brief 一个牛顿步中流体和孔隙体积误差的最大值
    OCP_DBL Verrmax;    ///< 类型：双精度浮点数，默认值未定义
};

/// \brief 牛顿迭代控制类
class ControlNR
{
public:
    /// \brief 设置控制参数
    /// \param src 包含控制参数的双精度浮点数向量
    void SetCtrlParam(const vector<OCP_DBL>& src);

    /// \brief 设置通信器
    /// \param domain 包含通信器和进程信息的域对象
    void SetupComm(const Domain& domain);

    /// \brief 设置下一个时间步的参数
    /// \param i 参数集合中的索引
    void SetNextTSTEP(const USI& i);

    /// \brief 获取最大饱和度变化
    /// \return 返回当前参数集中的最大饱和度变化值
    auto DSmax() const;

    /// \brief 获取最大压力变化
    /// \return 返回当前参数集中的最大压力变化值
    auto DPmax() const;

    /// \brief 检查牛顿迭代是否收敛
    /// \param NRs 牛顿迭代套件对象
    /// \param il 初始化列表，包含收敛检查所需的参数名
    /// \return 返回一个OCPNRStateC类型的对象，代表迭代是否收敛
    OCPNRStateC CheckConverge(const OCPNRsuite& NRs, const initializer_list<string>& il) const;

protected:
    MPI_Comm         myComm;    ///< 类型：MPI通信器，默认值未定义
    OCP_INT          numproc;   ///< 类型：整型，代表进程数，默认值未定义
    OCP_INT          myrank;    ///< 类型：整型，代表当前进程的排名，默认值未定义

protected:
    /// \brief 控制参数集合
    /// \details 存储了多个时间步的牛顿迭代控制参数
    vector<ControlNRParam> ps;  ///< 类型：ControlNRParam类的向量，默认值为空

    /// \brief 当前使用的控制参数
    /// \details 指向当前时间步使用的ControlNRParam对象
    const ControlNRParam* wp;   ///< 类型：指向ControlNRParam的指针，默认值为nullptr
};

#endif /* end if __OCPControlNR_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/05/2023      Create file                          */
/*----------------------------------------------------------------------------*/