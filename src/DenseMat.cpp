/*! \file    DenseMat.cpp
 *  \brief   Dense matrix-vector operations
 *  \author  Shizhe Li
 *  \date    Oct/21/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "DenseMat.hpp"


// WARNING: absolute sum!
OCP_DBL Dnorm1(const INT& N, OCP_DBL* x)
{

    return OCP_norm1(N, x);
    //const INT incx = 1;
    //return dasum_(&N, x, &incx);
}

OCP_DBL Dnorm2(const INT& N, OCP_DBL* x)
{
    return OCP_norm2(N, x);
    //const INT incx = 1;
    //return dnrm2_(&N, x, &incx);
}

void Dscalar(const INT& n, const OCP_DBL& alpha, OCP_DBL* x)
{

    OCP_scale(n, alpha, x);

    //// x = a x
    //const int incx = 1;
    //dscal_(&n, &alpha, x, &incx);
}

void Daxpy(const INT& n, const OCP_DBL& alpha, const OCP_DBL* x, OCP_DBL* y)
{

    OCP_axpy(n, alpha, x, y);

    //// y= ax +y
    //const int incx = 1, incy = 1;
    //daxpy_(&n, &alpha, x, &incx, y, &incy);
}

void DaABpbC(const INT&    m,
             const INT&    n,
             const INT&    k,
             const OCP_DBL& alpha,
             const OCP_DBL* A,
             const OCP_DBL* B,
             const OCP_DBL& beta,
             OCP_DBL*       C)
{
    /*  C' = alpha B'A' + beta C'
     *  A: m x k
     *  B: k x n
     *  C: m x n
     *  all column majored matrices, no tranpose
     *  A' in col-order in Fortran = A in row-order in C/Cpp
     */

    OCP_ABpC(m, n, k, A, B, C);

    //const char transa = 'N', transb = 'N';
    //dgemm_(&transa, &transb, &n, &m, &k, &alpha, B, &n, A, &k, &beta, C, &n);
}


void LUSolve(const INT& nrhs, const INT& N, OCP_DBL* A, OCP_DBL* b, INT* pivot)
{
    INT info;

#if OCPFLOATTYPEWIDTH == 64

    dgesv_(&N, &nrhs, A, &N, pivot, b, &N, &info);

    if (info < 0) {
        cout << "Wrong Input !" << endl;
    } else if (info > 0) {
        cout << "Singular Matrix !" << endl;
    }


#else
    OCP_ABORT("NOT AVAILABLE!");

#endif
}

INT SYSSolve(const INT& nrhs, const OCP_CHAR* uplo, const INT& N, OCP_DBL* A, OCP_DBL* b, INT* pivot, OCP_DBL* work, const INT& lwork)
{
    INT info;

#if OCPFLOATTYPEWIDTH == 64

    dsysv_(uplo, &N, &nrhs, A, &N, pivot, b, &N, work, &lwork, &info);
    if (info < 0) {
        cout << "Wrong Input !" << endl;
    } else if (info > 0) {
        cout << "Singular Matrix !" << endl;
    }

#else
    OCP_ABORT("NOT AVAILABLE!");

#endif

    return info;
}


void CalEigenSY(const INT& N, OCP_SIN* A, OCP_SIN* w, OCP_SIN* work, const INT& lwork)
{
    INT  info;
    INT  iwork[1] = { 0 };
    INT  liwork = 1;
    char uplo{ 'U' };
    char Nonly{ 'N' };

    ssyevd_(&Nonly, &uplo, &N, A, &N, w, work, &lwork, iwork, &liwork, &info);
    if (info > 0) {
        cout << "failed to compute eigenvalues!" << endl;
    }
}



void myDABpCp(const int& m,
    const int& n,
    const int& k,
    const double* A,
    const double* B,
    double* C,
    const int* flag,
    const int     N)
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            for (int p = 0; p < 3; p++) {
                if (flag[p] != 0) C[i * k + j] += A[i * n + p] * B[p * k + j];
            }
            for (int p = 0; p < 2; p++) {
                if (flag[p] != 0) {
                    for (int m = 0; m < N; m++) {
                        C[i * k + j] += A[i * n + 3 + p * (N + 1) + m] *
                            B[(3 + p * (N + 1) + m) * k + j];
                    }
                }
            }
        }
    }
}


void myDABpCp1(const int& m,
    const int& n,
    const int& k,
    const double* A,
    const double* B,
    double* C,
    const int* flag,
    const int     N)
{

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            int s = 0;
            for (int p = 0; p < 3; p++) {
                if (flag[p] != 0) {
                    C[i * k + j] += A[i * n + p] * B[s * k + j];
                    s++;
                }
            }
            for (int p = 0; p < 2; p++) {
                if (flag[p] != 0) {
                    for (int m = 0; m < N; m++) {
                        C[i * k + j] += A[i * n + 3 + p * (N + 1) + m] * B[s * k + j];
                        s++;
                    }
                }
            }
        }
    }
}

void myDABpCp2(const int& m,
    const int& n,
    const int& k,
    const double* A,
    const double* B,
    double* C,
    const int* flag,
    const int     N)
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            int s = 0;
            for (int p = 0; p < 3; p++) {
                if (flag[p] != 0) {
                    C[i * k + j] += A[i * n + s] * B[p * k + j];
                    s++;
                }
            }
            for (int p = 0; p < 2; p++) {
                if (flag[p] != 0) {
                    for (int m = 0; m < N; m++) {
                        C[i * k + j] += A[i * n + s] * B[(3 + p * (N + 1) + m) * k + j];
                        s++;
                    }
                }
            }
        }
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/21/2021      Create file                          */
/*----------------------------------------------------------------------------*/
