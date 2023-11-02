/*! \file    DenseMat.hpp
 *  \brief   Operations about small dense mat
 *  \author  Shizhe Li
 *  \date    Oct/24/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __DENSEMAT_HEADER__
#define __DENSEMAT_HEADER__

// Standard header files
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

////// BLAS functions

/// Scales a vector by a constant.
void dscal_(const int* n, const double* alpha, double* x, const int* incx);

/// Forms the dot product of two vectors.
double ddot_(const int* n, double* a, const int* inca, double* b, const int* incb);

/// Copies a vector, src, to a vector, dst.
int dcopy_(
    const int* n, const double* src, const int* incx, double* dst, const int* incy);

/// Constant times a vector plus a vector.
int daxpy_(const int*    n,
           const double* alpha,
           const double* x,
           const int*    incx,
           double*       y,
           const int*    incy);

/// Computes the Euclidean norm of a vector.
double dnrm2_(const int* n, double* x, const int* incx);

/// Computes the sum of the absolute values of a vector.
double dasum_(const int* n, double* x, const int* incx);

/// Finds the index of element having max absolute value.
int idamax_(const int* n, double* x, const int* incx);

/// Performs matrix-matrix operations C : = alpha * op(A) * op(B) + beta * C.
int dgemm_(const char*   transa,
           const char*   transb,
           const int*    m,
           const int*    n,
           const int*    k,
           const double* alpha,
           const double* A,
           const int*    lda,
           const double* B,
           const int*    ldb,
           const double* beta,
           double*       C,
           const int*    ldc);

////// LAPACK functions

/// Computes the solution to system of linear equations A * X = B for general matrices.
void dgesv_(const int* n,
    const int* nrhs,
    double* A,
    const int* lda,
    int* ipiv,
    double* b,
    const int* ldb,
    int* info);

/// Computes the solution to system of linear equations A * X = B for symm matrices.
void dsysv_(const char* uplo,
    const int* n,
    const int* nrhs,
    double* A,
    const int* lda,
    int* ipiv,
    double* b,
    const int* ldb,
    double* work,
    const int* lwork,
    int* info);

/// Computes the eigenvalues and, optionally, the leftand /or right eigenvectors for SY
/// matrices
void ssyevd_(const char* jobz,
    const char* uplo,
    const int* n,
    float* A,
    const int* lda,
    float* w,
    float* work,
    const int* lwork,
    int* iwork,
    const int* liwork,
    int* info);
}


/// Computes L1-norm of a vector.
OCP_DBL Dnorm1(const INT& N, OCP_DBL* x);

/// Computes L1-norm of a vector.
template <typename T1, typename T2>
T2 OCP_norm1(const T1& n, const T2* x)
{
    T2 tmp = 0;
    for (T1 i = 0; i < n; i++) {
        tmp += fabs(x[i]);
    }
    return tmp;
}


/// Computes L2-norm of a vector.
OCP_DBL Dnorm2(const INT& N, OCP_DBL* x);

/// Computes L2-norm of a vector.
template <typename T1, typename T2>
T2 OCP_norm2(const T1& n, const T2* x)
{
    T2 tmp = 0;
    for (T1 i = 0; i < n; i++) {
        tmp += x[i] * x[i];
    }
    return sqrt(tmp);
}


/// Scales a vector by a constant.
void Dscalar(const INT& n, const OCP_DBL& alpha, OCP_DBL* x);

/// Computes x = ax
template <typename T1, typename T2>
void OCP_scale(const T1& n, const T2& a, T2* x)
{
    for (T1 i = 0; i < n; i++) {
        x[i] *= a;
    }
}


/// Constant times a vector plus a vector.
void Daxpy(const INT& n, const OCP_DBL& alpha, const OCP_DBL* x, OCP_DBL* y);

/// Computes y = ax + y
template <typename T1, typename T2>
void OCP_axpy(const T1& n, const T2& a, const T2* x, T2* y)
{
    for (T1 i = 0; i < n; i++) {
        y[i] += a * x[i];
    }
}


/// Computes C' = alpha B'A' + beta C', all matrices are column-major.
void DaABpbC(const INT&    m,
             const INT&    n,
             const INT&    k,
             const OCP_DBL& alpha,
             const OCP_DBL* A,
             const OCP_DBL* B,
             const OCP_DBL& beta,
             OCP_DBL*       C);

/// Computes C = AB + C
template <typename T1, typename T2>
void OCP_ABpC(const T1& m, const T1& n, const T1& k, const T2* A, const T2* B, T2* C)
{
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


/// Computes y = a A x + b y
template <typename T1, typename T2>
void OCP_aAxpby(const T1& m, const T1& n, const T2& a, const T2* A, const T2* x, const T2& b, T2* y)
{
	for (T1 i = 0; i < m; i++) {
		y[i] = b * y[i];
		for (T1 j = 0; j < n; j++) {
			y[i] += a * A[i * n + j] * x[j];
		}
	}
}

/// Calls dgesv to solve the linear system for general matrices.
void LUSolve(const INT& nrhs, const INT& N, OCP_DBL* A, OCP_DBL* b, INT* pivot);

/// Calls dsysy to solve the linear system for symm matrices.
INT SYSSolve(const INT& nrhs, const OCP_CHAR* uplo, const INT& N, OCP_DBL* A, OCP_DBL* b, INT* pivot, OCP_DBL* work, const INT& lwork);

/// Calculate the minimal eigenvalue for symmetric matrix with mkl lapack
void CalEigenSY(const INT& N, OCP_SIN* A, OCP_SIN* w, OCP_SIN* work, const INT& lwork);


void myDABpCp(const int& m,
    const int& n,
    const int& k,
    const double* A,
    const double* B,
    double* C,
    const int* flag,
    const int     N);


void myDABpCp1(const int& m,
    const int& n,
    const int& k,
    const double* A,
    const double* B,
    double* C,
    const int* flag,
    const int     N);


void myDABpCp2(const int& m,
    const int& n,
    const int& k,
    const double* A,
    const double* B,
    double* C,
    const int* flag,
    const int     N);


/// Prints a vector.
template <typename T>
void PrintDX(const int& N, const T* x)
{
    for (int i = 0; i < N; i++) {
        cout << i << "   " << setprecision(16) << x[i] << endl;
    }
    cout << endl;
}


/// check NaN
template <typename T>
bool CheckNan(const int& N, const T* x)
{
    for (int i = 0; i < N; i++) {
        if (!isfinite(x[i])) {
            return false;
        }
    }
    return true;
}


/// swap value instead of pointer
template <typename T>
inline void OCPSwap(T a, T b, const int& n, T w)
{
    for (int i = 0; i < n; i++) {
        w[i] = a[i];
        a[i] = b[i];
        b[i] = w[i];
    }
}

#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/24/2021      Create file                          */
/*  Chensong Zhang      Jan/16/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/