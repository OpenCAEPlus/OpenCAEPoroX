/*! \file    SamgSolver.hpp
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
void SamgSolver::Allocate(const OCP_USI& max_nnz, const OCP_USI& maxDim)
{
    iA.resize(maxDim * blockdim + 1);
    jA.resize(max_nnz * blockdim * blockdim);
    A.resize(max_nnz * blockdim * blockdim);

    if (nsys > 1) {
        iu.resize(maxDim * blockdim);
        ip.resize(maxDim * blockdim);
        for (OCP_INT n = 0; n < maxDim; n++) {
            copy(iu_tmp.begin(), iu_tmp.end(), &iu[n * blockdim]);
            fill(ip.data() + n * blockdim, ip.data() + n * blockdim + blockdim, n + 1);
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
    //    nshalo += (sel[s].size() - 1) * blockdim;
    //    ipts[s + 1] = nshalo + 1;
    //}
    //isndlist.resize(nshalo);  // 1-based index
    //OCP_USI siter = 0;
    //for (USI s = 0; s < npsnd; s++) {
    //    for (USI i = 1; i < sel[s].size(); i++) {
    //        const OCP_USI bIds = sel[s][i] * blockdim + 1;
    //        for (USI j = 0; j < blockdim; j++) {
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
    //    nrhalo += (rel[r][2] - rel[r][1]) * blockdim;
    //    iptr[r + 1] = nrhalo + 1;
    //}
    //ireclist.resize(nrhalo);  // 1-based index
    //OCP_USI riter = 0;
    //for (USI r = 0; r < nprec; r++) {
    //    for (USI i = rel[r][1]; i < rel[r][2]; i++) {
    //        const OCP_USI bIdr = (i + actWellNum) * blockdim + 1;
    //        for (USI j = 0; j < blockdim; j++) {
    //            ireclist[riter++] = bIdr + j;
    //        }
    //    }
    //}
}


/// Assemble coefficient matrix.
void ScalarSamgSolver::AssembleMat(const vector<vector<USI>>& colId,
    const vector<vector<OCP_DBL>>& val,
    const OCP_USI& dim,
    vector<OCP_DBL>& rhs,
    vector<OCP_DBL>& u)
{
    nnu = dim;
    b   = rhs.data();
    x   = u.data();

    iA[0] = 1; // 1-baesd index
    for (OCP_USI i = 1; i < dim + 1; i++) {
        USI nnz_Row = colId[i - 1].size();
        iA[i] = iA[i - 1] + nnz_Row;

        const OCP_USI bId = iA[i - 1] - 1;
        // the first entry is in diagnal line all right
        if (val[i - 1][0] > 0) {           
            for (USI j = 0; j < nnz_Row; j++) {
                jA[bId + j] = global_index->at(colId[i - 1][j]) + 1;  // 1-baesd index
                A[bId + j]  = val[i - 1][j];
            }
        }
        else {
            for (USI j = 0; j < nnz_Row; j++) {
                jA[bId + j] = global_index->at(colId[i - 1][j]) + 1;  // 1-baesd index
                A[bId + j]  = -val[i - 1][j];
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
void VectorSamgSolver::AssembleMat(const vector<vector<USI>>& colId,
    const vector<vector<OCP_DBL>>& val,
    const OCP_USI& dim,
    vector<OCP_DBL>& rhs,
    vector<OCP_DBL>& u)
{
    nnu  = dim * blockdim;
    ndiu = nnu;
    ndip = nnu;
    b    = rhs.data();
    x    = u.data();

    const USI blockSize = blockdim * blockdim;
    // Assemble iA, jA, A
    iA[0] = 1; // 1-baesd index
    for (OCP_USI i = 1; i < dim + 1; i++) {
        const USI nnzR    = colId[i - 1].size();
        const OCP_USI bId = (i - 1) * blockdim;

        // iA
        for (USI c = 0; c < blockdim; c++)
            iA[bId + c + 1] = iA[bId + c] + nnzR * blockdim;

        // iA and A
        // copy each block
        for (USI j = 0; j < nnzR; j++) {
            for (USI c = 0; c < blockdim; c++) {
                const OCP_USI bIdc = iA[bId + c] - 1;
                for (USI c1 = 0; c1 < blockdim; c1++)
                    jA[bIdc + j * blockdim + c1] = global_index->at(colId[i - 1][j]) * blockdim + c1 + 1; // 1-baesd index

                const OCP_DBL* begin = &val[i - 1][0] + j * blockSize + c * blockdim;
                const OCP_DBL* end = begin + blockdim;
                copy(begin, end, &A[bIdc + j * blockdim]);
            }
        }
        // the first entry should be in diagnal line and positive
        for (USI c = 0; c < blockdim; c++) {
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