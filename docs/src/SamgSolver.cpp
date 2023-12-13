/*! \file    SamgSolver.cpp
 *  \brief   API of SANG Solver
 *  \author  Shizhe Li
 *  \date    Apr/07/2023
 *
 *  \note    use SAMG as linear solver
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifdef WITH_SAMG

#include "SamgSolver.hpp"


/// Set parameters.
void SamgSolver::SetupParam(const string& dir, const string& file)
{
    InitParam();
}


/// Initialize the Params for linear solver.
void SamgSolver::InitParam()
{

}


/// Allocate memoery for pardiso solver
void SamgSolver::Allocate(const OCPMatrix& mat)
{
    iA.resize(mat.maxDim * nb + 1);
    jA.resize(mat.max_nnz * nb * nb);
    A.resize(mat.max_nnz * nb * nb);

    if (nsys > 1) {
        iu.resize(mat.maxDim * nb);
        ip.resize(mat.maxDim * nb);
        for (OCP_INT n = 0; n < mat.maxDim; n++) {
            copy(iu_tmp.begin(), iu_tmp.end(), &iu[n * nb]);
            fill(ip.data() + n * nb, ip.data() + n * nb + nb, n + 1);
        }
    }
}


/// Calculate terms used in communication
void SamgSolver::CalCommTerm(const USI& actWellNum, const Domain* domain)
{
    // Two methods to calculate communication term

    // First, let SAMG calculate it automatically
    global_index = domain->CalGlobalIndex(actWellNum);

    // Or Second, give it directly
    // for send
    //const vector<vector<OCP_USI>>& sel = domain->send_element_loc;
    //npsnd = sel.size();
    //iranksnd.resize(npsnd);
    //ipts.resize(npsnd + 1, 1);
    //nshalo = 0;
    //for (USI s = 0; s < npsnd; s++) {
    //    iranksnd[s] = sel[s][0]; 
    //    nshalo += (sel[s].size() - 1) * nb;
    //    ipts[s + 1] = nshalo + 1;
    //}
    //isndlist.resize(nshalo);  // 1-based index
    //OCP_USI siter = 0;
    //for (USI s = 0; s < npsnd; s++) {
    //    for (USI i = 1; i < sel[s].size(); i++) {
    //        const OCP_USI bIds = sel[s][i] * nb + 1;
    //        for (USI j = 0; j < nb; j++) {
    //            isndlist[siter++] = bIds + j;
    //        }
    //    }
    //}
    //// for receive
    //const vector<vector<OCP_USI>>& rel = domain->recv_element_loc;
    //nprec = rel.size();
    //irankrec.resize(nprec);
    //iptr.resize(nprec + 1, 1);
    //nrhalo = 0;
    //for (USI r = 0; r < nprec; r++) {
    //    irankrec[r] = rel[r][0];
    //    nrhalo += (rel[r][2] - rel[r][1]) * nb;
    //    iptr[r + 1] = nrhalo + 1;
    //}
    //ireclist.resize(nrhalo);  // 1-based index
    //OCP_USI riter = 0;
    //for (USI r = 0; r < nprec; r++) {
    //    for (USI i = rel[r][1]; i < rel[r][2]; i++) {
    //        const OCP_USI bIdr = (i + actWellNum) * nb + 1;
    //        for (USI j = 0; j < nb; j++) {
    //            ireclist[riter++] = bIdr + j;
    //        }
    //    }
    //}
}


/// Assemble coefficient matrix.
void ScalarSamgSolver::AssembleMat(OCPMatrix& mat)
{
    nnu = mat.dim;
    b   = mat.b.data();
    x   = mat.u.data();

    iA[0] = 1; // 1-baesd index
    for (OCP_USI i = 1; i < mat.dim + 1; i++) {
        USI nnz_Row = mat.colId[i - 1].size();
        iA[i] = iA[i - 1] + nnz_Row;

        const OCP_USI bId = iA[i - 1] - 1;
        // the first entry is in diagnal line all right
        if (mat.val[i - 1][0] > 0) {
            for (USI j = 0; j < nnz_Row; j++) {
                jA[bId + j] = global_index->at(mat.colId[i - 1][j]) + 1;  // 1-baesd index
                A[bId + j]  = mat.val[i - 1][j];
            }
        }
        else {
            for (USI j = 0; j < nnz_Row; j++) {
                jA[bId + j] = global_index->at(mat.colId[i - 1][j]) + 1;  // 1-baesd index
                A[bId + j]  = -mat.val[i - 1][j];
            }
            b[i - 1] = -b[i - 1];
        }
    }

    nna = iA[nnu] - 1;
}


/// Solve the linear system.
OCP_INT SamgSolver::Solve()
{
    //exemplary demonstration of how to set secondary control parameters of SAMG:
     SAMG_INT levelx = 25; 
     SAMG_SET_LEVELX(&levelx);

    // SAMGP_OIL(&nnu, &nna, &nsys,
    //           &iA[0], &jA[0], &A[0], &b[0], &x[0], iu.data(), &ndiu, ip.data(), &ndip, &samg_matrix,
    //           &res_in, &res_out, &ncyc_done, &ierr,
    //           &ifirst, &eps, &ncyc, &iswtch,
    //           &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
    //           &chktol, &idump, &iout,
    //           nunknown_description.data(), &noil_approach,
    //           &noil_cyc, &noil_preparation,
    //           &nshalo, &npsnd, iranksnd.data(), ipts.data(), isndlist.data(),
    //           &nrhalo, &nprec, irankrec.data(), iptr.data(), ireclist.data(), &myComm);

	nrhalo = -1;
	SAMGP_PCRS_OIL(&nnu, &nna, &nsys,
		&iA[0], &jA[0], &A[0], &b[0], &x[0], iu.data(), &ndiu, ip.data(), &ndip, &samg_matrix,
		&res_in, &res_out, &ncyc_done, &ierr,
		&ifirst, &eps, &ncyc, &iswtch,
		&a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
		&chktol, &idump, &iout,
		nunknown_description.data(), &noil_approach,
		&noil_cyc, &noil_preparation,
		&nrhalo, &myComm);

    return ncyc_done;
}


/// Assemble coefficient matrix.
void VectorSamgSolver::AssembleMat(OCPMatrix& mat)
{
    nnu  = mat.dim * nb;
    ndiu = nnu;
    ndip = nnu;
    b    = mat.b.data();
    x    = mat.u.data();

    const USI blockSize = nb * nb;
    // Assemble iA, jA, A
    iA[0] = 1; // 1-baesd index
    for (OCP_USI i = 1; i < mat.dim + 1; i++) {
        const USI nnzR    = mat.colId[i - 1].size();
        const OCP_USI bId = (i - 1) * nb;

        // iA
        for (USI c = 0; c < nb; c++)
            iA[bId + c + 1] = iA[bId + c] + nnzR * nb;

        // iA and A
        // copy each block
        for (USI j = 0; j < nnzR; j++) {
            for (USI c = 0; c < nb; c++) {
                const OCP_USI bIdc = iA[bId + c] - 1;
                for (USI c1 = 0; c1 < nb; c1++)
                    jA[bIdc + j * nb + c1] = global_index->at(mat.colId[i - 1][j]) * nb + c1 + 1; // 1-baesd index

                const OCP_DBL* begin = &mat.val[i - 1][0] + j * blockSize + c * nb;
                const OCP_DBL* end = begin + nb;
                copy(begin, end, &A[bIdc + j * nb]);
            }
        }
        // the first entry should be in diagnal line and positive
        for (USI c = 0; c < nb; c++) {
            const OCP_USI bIdc = iA[bId + c] - 1;
            swap(jA[bIdc], jA[bIdc + c]);
            swap(A[bIdc], A[bIdc + c]);
            if (A[bIdc] < 0) {
                const OCP_USI eIdc = iA[bId + c + 1] - 1;
                for (OCP_USI j = bIdc; j < eIdc; j++)
                    A[j] = -A[j];
                b[bId + c] = -b[bId + c];
            }
        }
    }

    nna = iA[nnu] - 1;
}

#endif // WITH_SAMG


 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Apr/07/2023      Create file                          */
 /*----------------------------------------------------------------------------*/