```cpp
/*! \file    DenseMat.hpp
 *  \brief   本文件包含了小型稠密矩阵操作的函数定义。
 *  \author  Shizhe Li
 *  \date    Oct/24/2021
 *
 *  本文件定义了多个与小型稠密矩阵操作相关的函数，包括BLAS和LAPACK函数的封装，
 *  矩阵的范数计算，矩阵与向量的操作，以及线性系统的求解等。
 */

#ifndef __DENSEMAT_HEADER__
#define __DENSEMAT_HEADER__

// 标准头文件
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <math.h>
#include <numeric>
#include "OCPDataType.hpp"
#include "UtilError.hpp"

using namespace std;

extern "C" {

// BLAS 函数

/*!
 * \brief   向量乘以常数。
 * \param n     向量的元素数量。
 * \param alpha 乘数。
 * \param x     目标向量。
 * \param incx  向量中每个元素之间的存储间隔。
 */
void dscal_(const int* n, const double* alpha, double* x, const int* incx);

/*!
 * \brief   计算两个向量的点积。
 * \param n     向量的元素数量。
 * \param a     第一个向量。
 * \param inca  第一个向量中每个元素之间的存储间隔。
 * \param b     第二个向量。
 * \param incb  第二个向量中每个元素之间的存储间隔。
 * \return      返回点积结果。
 */
double ddot_(const int* n, double* a, const int* inca, double* b, const int* incb);

/*!
 * \brief   复制向量。
 * \param n     向量的元素数量。
 * \param src   源向量。
 * \param incx  源向量中每个元素之间的存储间隔。
 * \param dst   目标向量。
 * \param incy  目标向量中每个元素之间的存储间隔。
 * \return      返回操作状态码。
 */
int dcopy_(const int* n, const double* src, const int* incx, double* dst, const int* incy);

/*!
 * \brief   向量加法。
 * \param n     向量的元素数量。
 * \param alpha 加数向量的乘数。
 * \param x     加数向量。
 * \param incx  加数向量中每个元素之间的存储间隔。
 * \param y     被加向量。
 * \param incy  被加向量中每个元素之间的存储间隔。
 * \return      返回操作状态码。
 */
int daxpy_(const int* n, const double* alpha, const double* x, const int* incx, double* y, const int* incy);

/*!
 * \brief   计算向量的欧几里得范数。
 * \param n     向量的元素数量。
 * \param x     目标向量。
 * \param incx  向量中每个元素之间的存储间隔。
 * \return      返回向量的欧几里得范数。
 */
double dnrm2_(const int* n, double* x, const int* incx);

/*!
 * \brief   计算向量的绝对值之和。
 * \param n     向量的元素数量。
 * \param x     目标向量。
 * \param incx  向量中每个元素之间的存储间隔。
 * \return      返回向量的绝对值之和。
 */
double dasum_(const int* n, double* x, const int* incx);

/*!
 * \brief   找到具有最大绝对值的元素的索引。
 * \param n     向量的元素数量。
 * \param x     目标向量。
 * \param incx  向量中每个元素之间的存储间隔。
 * \return      返回具有最大绝对值的元素的索引。
 */
int idamax_(const int* n, double* x, const int* incx);

/*!
 * \brief   执行矩阵-矩阵运算。
 * \param transa 指定矩阵A的转置操作。
 * \param transb 指定矩阵B的转置操作。
 * \param m      矩阵A的行数。
 * \param n      矩阵B的列数。
 * \param k      矩阵A的列数和矩阵B的行数。
 * \param alpha  标量乘数。
 * \param A      矩阵A。
 * \param lda    A的前导维度。
 * \param B      矩阵B。
 * \param ldb    B的前导维度。
 * \param beta   标量加数。
 * \param C      目标矩阵C。
 * \param ldc    C的前导维度。
 * \return       返回操作状态码。
 */
int dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k, const double* alpha, const double* A, const int* lda, const double* B, const int* ldb, const double* beta, double* C, const int* ldc);

// LAPACK 函数

/*!
 * \brief   为一般矩阵计算线性方程组的解。
 * \param n     方程组的未知数个数。
 * \param nrhs  右侧矩阵的列数。
 * \param A     系数矩阵A。
 * \param lda   A的前导维度。
 * \param ipiv  存储枢轴信息的数组。
 * \param b     右侧矩阵B。
 * \param ldb   B的前导维度。
 * \param info  输出参数，返回操作状态码。
 */
void dgesv_(const int* n, const int* nrhs, double* A, const int* lda, int* ipiv, double* b, const int* ldb, int* info);

/*!
 * \brief   为对称矩阵计算线性方程组的解。
 * \param uplo  指定矩阵A中存储的是上三角部分还是下三角部分。
 * \param n     方程组的未知数个数。
 * \param nrhs  右侧矩阵的列数。
 * \param A     系数矩阵A。
 * \param lda   A的前导维度。
 * \param ipiv  存储枢轴信息的数组。
 * \param b     右侧矩阵B。
 * \param ldb   B的前导维度。
 * \param work  工作数组。
 * \param lwork 工作数组的长度。
 * \param info  输出参数，返回操作状态码。
 */
void dsysv_(const char* uplo, const int* n, const int* nrhs, double* A, const int* lda, int* ipiv, double* b, const int* ldb, double* work, const int* lwork, int* info);

/*!
 * \brief   为对称矩阵计算特征值和（可选的）左/右特征向量。
 * \param jobz  指定是否计算特征向量。
 * \param uplo  指定矩阵A中存储的是上三角部分还是下三角部分。
 * \param n     矩阵的阶数。
 * \param A     系数矩阵A。
 * \param lda   A的前导维度。
 * \param w     存储特征值的数组。
 * \param work  工作数组。
 * \param lwork 工作数组的长度。
 * \param iwork 整数型工作数组。
 * \param liwork 整数型工作数组的长度。
 * \param info  输出参数，返回操作状态码。
 */
void ssyevd_(const char* jobz, const char* uplo, const int* n, float* A, const int* lda, float* w, float* work, const int* lwork, int* iwork, const int* liwork, int* info);

} // extern "C"

// 向量范数计算

/*!
 * \brief   计算向量的L1范数。
 * \param N 向量的元素数量。
 * \param x 目标向量。
 * \return  返回向量的L1范数。
 */
OCP_DBL Dnorm1(const INT& N, OCP_DBL* x);

/*!
 * \brief   计算向量的L1范数。
 * \param n 向量的元素数量。
 * \param x 目标向量。
 * \return  返回向量的L1范数。
 */
template <typename T1, typename T2>
T2 OCP_norm1(const T1& n, const T2* x) {
    T2 tmp = 0;
    for (T1 i = 0; i < n; i++) {
        tmp += fabs(x[i]);
    }
    return tmp;
}

/*!
 * \brief   计算向量的L2范数。
 * \param N 向量的元素数量。
 * \param x 目标向量。
 * \return  返回向量的L2范数。
 */
OCP_DBL Dnorm2(const INT& N, OCP_DBL* x);

/*!
 * \brief   计算向量的L2范数。
 * \param n 向量的元素数量。
 * \param x 目标向量。
 * \return  返回向量的L2范数。
 */
template <typename T1, typename T2>
T2 OCP_norm2(const T1& n, const T2* x) {
    T2 tmp = 0;
    for (T1 i = 0; i < n; i++) {
        tmp += x[i] * x[i];
    }
    return sqrt(tmp);
}

// 向量和矩阵的操作

/*!
 * \brief   向量乘以常数。
 * \param n     向量的元素数量。
 * \param alpha 乘数。
 * \param x     目标向量。
 */
void Dscalar(const INT& n, const OCP_DBL& alpha, OCP_DBL* x);

/*!
 * \brief   向量乘以常数。
 * \param n 向量的元素数量。
 * \param a 乘数。
 * \param x 目标向量。
 */
template <typename T1, typename T2>
void OCP_scale(const T1& n, const T2& a, T2* x) {
    for (T1 i = 0; i < n; i++) {
        x[i] *= a;
    }
}

/*!
 * \brief   向量加法。
 * \param n     向量的元素数量。
 * \param alpha 加数向量的乘数。
 * \param x     加数向量。
 * \param y     被加向量。
 */
void Daxpy(const INT& n, const OCP_DBL& alpha, const OCP_DBL* x, OCP_DBL* y);

/*!
 * \brief   向量加法。
 * \param n 向量的元素数量。
 * \param a 加数向量的乘数。
 * \param x 加数向量。
 * \param y 被加向量。
 */
template <typename T1, typename T2>
void OCP_axpy(const T1& n, const T2& a, const T2* x, T2* y) {
    for (T1 i = 0; i < n; i++) {
        y[i] += a * x[i];
    }
}

// 矩阵和矩阵的操作

/*!
 * \brief   矩阵乘法并加到另一矩阵。
 * \param m     结果矩阵的行数。
 * \param n     结果矩阵的列数。
 * \param k     A矩阵的列数和B矩阵的行数。
 * \param alpha A矩阵与B矩阵相乘的结果的乘数。
 * \param A     矩阵A。
 * \param B     矩阵B。
 * \param beta  结果矩阵的乘数。
 * \param C     结果矩阵C。
 */
void DaABpbC(const INT& m, const INT& n, const INT& k, const OCP_DBL& alpha, const OCP_DBL* A, const OCP_DBL* B, const OCP_DBL& beta, OCP_DBL* C);

/*!
 * \brief   矩阵乘法并加到另一矩阵。
 * \param m 结果矩阵的行数。
 * \param n 结果矩阵的列数。
 * \param k A矩阵的列数和B矩阵的行数。
 * \param A 矩阵A。
 * \param B 矩阵B。
 * \param C 结果矩阵C。
 */
template <typename T1, typename T2>
void OCP_ABpC(const T1& m, const T1& n, const T1& k, const T2* A, const T2* B, T2* C) {
    // C = AB + C
    // A: m*k  B:k*n  C:m*n
    // all matrix are row majored matrices
    for (T1 i = 0; i < m; i++) {
        for (T1 j = 0; j < n; j++) {
            for (T1 l = 0; l < k; l++) {
                C[i * n + j] += A[i * k + l] * B[l * n + j];
            }
        }
    }
}

// 其他辅助函数

/*!
 * \brief   打印向量。
 * \param N 向量的元素数量。
 * \param x 目标向量。
 */
template <typename T>
void PrintDX(const int& N, const T* x) {
    for (int i = 0; i < N; i++) {
        cout << i << "   " << setprecision(16) << x[i] << endl;
    }
    cout << endl;
}

/*!
 * \brief   检查向量中是否存在NaN值。
 * \param N 向量的元素数量。
 * \param x 目标向量。
 * \return  如果存在NaN值，返回false，否则返回true。
 */
template <typename T>
bool CheckNan(const int& N, const T* x) {
    for (int i = 0; i < N; i++) {
        if (!isfinite(x[i])) {
            return false;
        }
    }
    return true;
}

/*!
 * \brief   交换两个向量的值。
 * \param a 向量a。
 * \param b 向量b。
 * \param n 向量的元素数量。
 * \param w 临时工作向量。
 */
template <typename T>
inline void OCPSwap(T a, T b, const int& n, T w) {
    for (int i = 0; i < n; i++) {
        w[i] = a[i];
        a[i] = b[i];
        b[i] = w[i];
    }
}

// 自定义函数

/*!
 * \brief   自定义的矩阵乘法并加到另一矩阵。
 * \param m     结果矩阵的行数。
 * \param n     结果矩阵的列数。
 * \param k     A矩阵的列数和B矩阵的行数。
 * \param A     矩阵A。
 * \param B     矩阵B。
 * \param C     结果矩阵C。
 * \param flag  操作标志。
 * \param N     标志的数量。
 */
void myDABpCp(const int& m, const int& n, const int& k, const double* A, const double* B, double* C, const int* flag, const int N);

/*!
 * \brief   自定义的矩阵乘法并加到另一矩阵的变体1。
 * \param m     结果矩阵的行数。
 * \param n     结果矩阵的列数。
 * \param k     A矩阵的列数和B矩阵的行数。
 * \param A     矩阵A。
 * \param B     矩阵B.
 * \param C     结果矩阵C。
 * \param flag  操作标志。
 * \param N     标志的数量。
 */

void myDABpCp1(const int& m,
    const int& n,
    const int& k,
    const double* A,
    const double* B,
    double* C,
    const int* flag,
    const int     N);

/*!
 * \brief   自定义的矩阵乘法并加到另一矩阵的变体2。
 * \param m     结果矩阵的行数。
 * \param n     结果矩阵的列数。
 * \param k     A矩阵的列数和B矩阵的行数。
 * \param A     矩阵A。
 * \param B     矩阵B.
 * \param C     结果矩阵C。
 * \param flag  操作标志。
 * \param N     标志的数量。
 */
void myDABpCp2(const int& m,
    const int& n,
    const int& k,
    const double* A,
    const double* B,
    double* C,
    const int* flag,
    const int     N);



#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/24/2021      Create file                          */
/*  Chensong Zhang      Jan/16/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/