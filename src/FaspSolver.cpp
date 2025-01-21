/*! \file    FaspSolver.cpp
 *  \brief   FaspSolver class definition
 *  \author  Shizhe Li
 *  \date    Nov/22/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "../config/config.hpp"


#ifdef OCP_USE_FASP

#include <math.h>

#include "FaspSolver.hpp"

void FaspSolver::SetupParam(const string& dir, const string& file)
{
    solveDir  = dir;
    solveFile = file;

    string myfile = solveDir + solveFile;
    InitParam(); // Set default solver parameters
    ifstream ifs(myfile);
    if (!ifs.is_open()) {
        cout << "The input file " << myfile << " is missing!" << endl;
        myfile = solveDir + "../conf/csr.fasp";
        ifs.open(myfile);
        if (!ifs.is_open()) {
            cout << "The input file " << myfile << " is missing!" << endl;
            cout << "Using the default parameters of FASP" << endl;
        } else {
            ifs.close();
            cout << "Using the input file " << myfile << endl;
            fasp_param_input(myfile.data(), &inParam);
        }
    } else {
        ifs.close(); // if file has been opened, close it first
        fasp_param_input(myfile.data(), &inParam);
    }
    fasp_param_init(&inParam, &itsParam, &amgParam, &iluParam, &swzParam);
    //cout << endl << "OCP long double : " << sizeof(OCP_DBL) << endl;
    //cout << endl << sizeof(_Quad) << endl;
}


ScalarFaspSolver::ScalarFaspSolver(const string& dir, const string& file, const OCPMatrix& mat)
{
    SetupParam(dir, file);
    Allocate(mat);
}


void ScalarFaspSolver::AssembleMat(OCPMatrix& mat, const Domain* domain)
{
    // b & x
    b.row = mat.dim;
    b.val = mat.b.data();
    x.row = mat.dim;
    x.val = mat.u.data();
    // A
    OCP_USI nnz = 0;
    for (OCP_USI i = 0; i < mat.dim; i++) {
        nnz += mat.colId[i].size();
    }

    A.row = mat.dim;
    A.col = mat.dim;
    A.nnz = nnz;

    // IA
    A.IA[0] = 0;
    for (OCP_USI i = 1; i < mat.dim + 1; i++) {
        USI nnz_Row = mat.colId[i - 1].size();
        A.IA[i] = A.IA[i - 1] + nnz_Row;

        copy(mat.colId[i - 1].begin(), mat.colId[i - 1].end(), &A.JA[A.IA[i - 1]]);
        copy(mat.val[i - 1].begin(), mat.val[i - 1].end(), &A.val[A.IA[i - 1]]);
    }
}

OCP_INT ScalarFaspSolver::Solve()
{
    OCP_INT status = FASP_SUCCESS;

    const OCP_INT print_level = inParam.print_level;
    const OCP_INT solver_type = inParam.solver_type;
    const OCP_INT precond_type = inParam.precond_type;
    const OCP_INT output_type = inParam.output_type;

    if (output_type) {
        const char* outputfile = "out/Solver.out";
        printf("Redirecting outputs to file: %s ...\n", outputfile);
        freopen(outputfile, "w", stdout); // open a file for stdout
    }

    // Preconditioned Krylov methods
    if (solver_type >= 1 && solver_type <= 20) {

        // Using no preconditioner for Krylov iterative methods
        if (precond_type == PREC_NULL) {
            status = fasp_solver_dcsr_krylov(&A, &b, &x, &itsParam);
        }

        // Using diag(A) as preconditioner for Krylov iterative methods
        else if (precond_type == PREC_DIAG) {
            status = fasp_solver_dcsr_krylov_diag(&A, &b, &x, &itsParam);
        }

        // Using AMG as preconditioner for Krylov iterative methods
        else if (precond_type == PREC_AMG || precond_type == PREC_FMG) {
            if (print_level > PRINT_NONE) fasp_param_amg_print(&amgParam);
            status = fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itsParam, &amgParam);
        }

        // Using ILU as preconditioner for Krylov iterative methods
        else if (precond_type == PREC_ILU) {
            if (print_level > PRINT_NONE) fasp_param_ilu_print(&iluParam);
            status = fasp_solver_dcsr_krylov_ilu(&A, &b, &x, &itsParam, &iluParam);
        }

        // Undefined iterative methods
        else {
            printf("### ERROR: Wrong preconditioner type %d!!!\n", precond_type);
            status = ERROR_SOLVER_PRECTYPE;
        }

    }

    // AMG as the iterative solver
    else if (solver_type == SOLVER_AMG) {
        if (print_level > PRINT_NONE) fasp_param_amg_print(&amgParam);
        fasp_solver_amg(&A, &b, &x, &amgParam);
    }

    // Full AMG as the iterative solver
    else if (solver_type == SOLVER_FMG) {
        if (print_level > PRINT_NONE) fasp_param_amg_print(&amgParam);
        fasp_solver_famg(&A, &b, &x, &amgParam);
    }

#if WITH_MUMPS // use MUMPS directly
    else if (solver_type == SOLVER_MUMPS) {
        status = fasp_solver_mumps(&A, &b, &x, print_level);
        if (status >= 0) status = 1; // Direct solver returns 1
    }
#endif

#if WITH_SuperLU // use SuperLU directly
    else if (solver_type == SOLVER_SUPERLU) {
        status = fasp_solver_superlu(&A, &b, &x, print_level);
        if (status >= 0) status = 1; // Direct solver returns 1
    }
#endif

#if WITH_UMFPACK // use UMFPACK directly
    else if (solver_type == SOLVER_UMFPACK) {
        // Need to sort the matrix A for UMFPACK to work
        dCSRmat A_trans = fasp_dcsr_create(A.row, A.col, A.nnz);
        fasp_dcsr_transz(&A, NULL, &A_trans);
        fasp_dcsr_sort(&A_trans);
        status = fasp_solver_umfpack(&A_trans, &b, &x, print_level);
        fasp_dcsr_free(&A_trans);
        if (status >= 0) status = 1; // Direct solver returns 1
    }
#endif

#ifdef OCP_USE_PARDISO // use PARDISO directly
    else if (solver_type == SOLVER_PARDISO) {
        fasp_dcsr_sort(&A);
        status = fasp_solver_pardiso(&A, &b, &x, print_level);
        if (status >= 0) status = 1; // Direct solver returns 1
    }
#endif

    else {
        printf("### ERROR: Wrong solver type %d!!!\n", solver_type);
        status = ERROR_SOLVER_TYPE;
    }

    if (status < 0) {
        printf("### ERROR: Solver failed! Exit status = %d.\n\n", status);
    }

    if (output_type) fclose(stdout);
    return status;
}


void ScalarFaspSolver::InitParam()
{
    // Input/output
    inParam.print_level = PRINT_MIN;
    inParam.output_type = 0;

    // Problem information
    inParam.solver_type  = SOLVER_VFGMRES;
    inParam.precond_type = PREC_AMG;
    inParam.stop_type    = STOP_REL_RES;

    // LinearSolver parameters
    inParam.itsolver_tol   = 1e-4;
    inParam.itsolver_maxit = 100;
    inParam.restart        = 30;

    // ILU method parameters
    inParam.ILU_type    = ILUk;
    inParam.ILU_lfil    = 0;
    inParam.ILU_droptol = 0.001;
    inParam.ILU_relax   = 0;
    inParam.ILU_permtol = 0.0;

    // Schwarz method parameters
    inParam.SWZ_mmsize    = 200;
    inParam.SWZ_maxlvl    = 2;
    inParam.SWZ_type      = 1;
    inParam.SWZ_blksolver = SOLVER_DEFAULT;

    // AMG method parameters
    inParam.AMG_type                = CLASSIC_AMG;
    inParam.AMG_levels              = 20;
    inParam.AMG_cycle_type          = V_CYCLE;
    inParam.AMG_smoother            = SMOOTHER_GS;
    inParam.AMG_smooth_order        = CF_ORDER;
    inParam.AMG_presmooth_iter      = 1;
    inParam.AMG_postsmooth_iter     = 1;
    inParam.AMG_relaxation          = 1.0;
    inParam.AMG_coarse_dof          = 500;
    inParam.AMG_coarse_solver       = 0;
    inParam.AMG_tol                 = 1e-6;
    inParam.AMG_maxit               = 1;
    inParam.AMG_ILU_levels          = 0;
    inParam.AMG_SWZ_levels          = 0;
    inParam.AMG_coarse_scaling      = OFF;
    inParam.AMG_amli_degree         = 1;
    inParam.AMG_nl_amli_krylov_type = 2;

    // Classical AMG specific
    inParam.AMG_coarsening_type      = 1;
    inParam.AMG_interpolation_type   = 1;
    inParam.AMG_max_row_sum          = 0.9;
    inParam.AMG_strong_threshold     = 0.3;
    inParam.AMG_truncation_threshold = 0.2;
    inParam.AMG_aggressive_level     = 0;
    inParam.AMG_aggressive_path      = 1;

    // Aggregation AMG specific
    inParam.AMG_aggregation_type   = PAIRWISE;
    inParam.AMG_quality_bound      = 8.0;
    inParam.AMG_pair_number        = 2;
    inParam.AMG_strong_coupled     = 0.25;
    inParam.AMG_max_aggregation    = 9;
    inParam.AMG_tentative_smooth   = 0.67;
    inParam.AMG_smooth_filter      = ON;
    inParam.AMG_smooth_restriction = ON;
}


void ScalarFaspSolver::Allocate(const OCPMatrix& mat)
{
    A = fasp_dcsr_create(mat.maxDim, mat.maxDim, mat.max_nnz);
}


VectorFaspSolver::VectorFaspSolver(const string& dir, const string& file, const OCPMatrix& mat)
{
    SetupParam(dir, file);
    Allocate(mat);
}


void VectorFaspSolver::AssembleMat(OCPMatrix& mat, const Domain* domain)
{
    const OCP_USI nrow = mat.dim * mat.nb;
    // b & x
    b.row = nrow;
    b.val = mat.b.data();
    x.row = nrow;
    x.val = mat.u.data(); // x will be set to zero later

    // fsc & order
    fsc.row = nrow;
    order.row = nrow;

    // nnz
    OCP_USI nnz = 0;
    for (OCP_USI i = 0; i < mat.dim; i++) {
        nnz += mat.colId[i].size();
    }

    // Asc
    Asc.ROW = mat.dim;
    Asc.COL = mat.dim;
    Asc.nb = mat.nb;
    Asc.NNZ = nnz;

    // A
    A.ROW = mat.dim;
    A.COL = mat.dim;
    A.nb = mat.nb;
    A.NNZ = nnz;

    const USI block_size = mat.nb * mat.nb;
    A.IA[0] = 0;
    for (OCP_USI i = 1; i < mat.dim + 1; i++) {
        USI nnb_Row = mat.colId[i - 1].size();
        A.IA[i] = A.IA[i - 1] + nnb_Row;

        copy(mat.colId[i - 1].begin(), mat.colId[i - 1].end(), &A.JA[A.IA[i - 1]]);
        copy(mat.val[i - 1].begin(), mat.val[i - 1].end(), &A.val[A.IA[i - 1] * block_size]);
    }
}

OCP_INT VectorFaspSolver::Solve()
{
    OCP_INT status = FASP_SUCCESS;

    // Set local parameters
    const OCP_INT print_level = inParam.print_level;
    const OCP_INT solver_type = inParam.solver_type;
    const OCP_INT precond_type = inParam.precond_type;
    const OCP_INT output_type = inParam.output_type;

#ifdef OCP_USE_FASP4BLKOIL || WITH_FASPCPR // Currently, only fasp4blkoil requires decoupling
    const OCP_INT decoup_type = inParam.decoup_type;
#endif

    if (output_type) {
        const char* outputfile = "../output/test.out";
        printf("Redirecting outputs to file: %s ...\n", outputfile);
        freopen(outputfile, "w", stdout); // open a file for stdout
    }

    fasp_dvec_set(x.row, &x, 0);

    // Preconditioned Krylov methods
    if (solver_type >= 1 && solver_type <= 10) {

        // Preconditioned Krylov methods in BSR format
        switch (precond_type) {
        case PC_NULL:
            status = fasp_solver_dbsr_krylov(&A, &b, &x, &itsParam);
            break;
        case PC_DIAG:
            status = fasp_solver_dbsr_krylov_diag(&A, &b, &x, &itsParam);
            break;
        case PC_BILU:
            status = fasp_solver_dbsr_krylov_ilu(&A, &b, &x, &itsParam, &iluParam);
            break;

#if WITH_FASPCPR //! FASPCPR solver (i.e., CPR or ASCPR preconditioners) added by
            //! zhaoli, 2022.12.11
        case PC_FASP1:
            Decoupling(&A, &b, &Asc, &fsc, &order, Dmat.data(), decoup_type);
            status = FASP_BSRSOL_ASCPR(&Asc, &fsc, &x, &itsParam, &iluParam,
                &amgParam, 0);
            break;

        case PC_FASP1_SHARE:
            Decoupling(&A, &b, &Asc, &fsc, &order, Dmat.data(), decoup_type);
            status = FASP_BSRSOL_ASCPR(&Asc, &fsc, &x, &itsParam, &iluParam,
                &amgParam, RESET_CONST);
            break;
#endif

#ifdef OCP_USE_FASP4BLKOIL
        case PC_FASP1:
            Decoupling(&A, &b, &Asc, &fsc, &order, Dmat.data(), decoup_type);
#if WITH_FASP4CUDA // zhaoli 2022.04.04
            status = fasp_solver_dbsr_krylov_FASP1_cuda_interface(
                &Asc, &fsc, &x, &itsParam, &iluParam, &amgParam, NULL, &order);
#else
            status = fasp_solver_dbsr_krylov_FASP1a(
                &Asc, &fsc, &x, &itsParam, &iluParam, &amgParam, NULL, &order);
#endif
            break;
        case PC_FASP1_SHARE: // zhaoli 2021.03.24
            Decoupling(&A, &b, &Asc, &fsc, &order, Dmat.data(), decoup_type);
#if WITH_FASP4CUDA
            status = fasp_solver_dbsr_krylov_FASP1_cuda_share_interface(
                &Asc, &fsc, &x, &itsParam, &iluParam, &amgParam, NULL, &order,
                RESET_CONST);
#else
            status = fasp_solver_dbsr_krylov_FASP1a_share_interface(
                &Asc, &fsc, &x, &itsParam, &iluParam, &amgParam, NULL, &order,
                RESET_CONST);
#endif
            break;
        case PC_FASP2:
            Decoupling(&A, &b, &Asc, &fsc, &order, Dmat.data(), decoup_type);
            status = fasp_solver_dbsr_krylov_FASP2(
                &Asc, &fsc, &x, &itsParam, &iluParam, &amgParam, NULL, &order);
            break;
        case PC_FASP3:
            Decoupling(&A, &b, &Asc, &fsc, &order, Dmat.data(), decoup_type);
            status = fasp_solver_dbsr_krylov_FASP3(
                &Asc, &fsc, &x, &itsParam, &iluParam, &amgParam, NULL, &order);
            break;
        case PC_FASP4:
            Decoupling(&A, &b, &Asc, &fsc, &order, Dmat.data(), decoup_type);
#if WITH_FASP4CUDA
            status = fasp_solver_dbsr_krylov_FASP4_cuda(
                &Asc, &fsc, &x, &itsParam, &iluParam, &amgParam, NULL, &order);
#else
            status = fasp_solver_dbsr_krylov_FASP4(
                &Asc, &fsc, &x, &itsParam, &iluParam, &amgParam, NULL, &order);
#endif
            break;
        case PC_FASP4_SHARE: // zhaoli 2021.04.24
            Decoupling(&A, &b, &Asc, &fsc, &order, Dmat.data(), decoup_type);
#if WITH_FASP4CUDA // zhaoli 2022.08.03
            status = fasp_solver_dbsr_krylov_FASP4_cuda_share_interface(
                &Asc, &fsc, &x, &itsParam, &iluParam, &amgParam, NULL, &order,
                RESET_CONST);
#else
            status = fasp_solver_dbsr_krylov_FASP4_share_interface(
                &Asc, &fsc, &x, &itsParam, &iluParam, &amgParam, NULL, &order,
                RESET_CONST);
#endif
            break;
        case PC_FASP5:
            Decoupling(&A, &b, &Asc, &fsc, &order, Dmat.data(), decoup_type);
            status = fasp_solver_dbsr_krylov_FASP5(
                &Asc, &fsc, &x, &itsParam, &iluParam, &amgParam, NULL, &order);
            break;
#endif
        default:
            OCP_WARNING("fasp4blkoil was not linked correctly!");
            OCP_ABORT("Preconditioner type " + to_string(precond_type) +
                " not supported!");
        }
        fill(Dmat.begin(), Dmat.end(), 0.0);
    }

#if WITH_MUMPS // use MUMPS directly
    else if (solver_type == SOLVER_MUMPS) {
        dCSRmat Acsr = fasp_format_dbsr_dcsr(&A);
        status = fasp_solver_mumps(&Acsr, &b, &x, print_level);
        fasp_dcsr_free(&Acsr);
        if (status >= 0) status = 1; // Direct solver returns 1
    }
#endif

#if WITH_SuperLU // use SuperLU directly
    else if (solver_type == SOLVER_SUPERLU) {
        dCSRmat Acsr = fasp_format_dbsr_dcsr(&A);
        status = fasp_solver_superlu(&Acsr, &b, &x, print_level);
        fasp_dcsr_free(&Acsr);
        if (status >= 0) status = 1; // Direct solver returns 1
    }
#endif

#if WITH_UMFPACK // use UMFPACK directly
    else if (solver_type == SOLVER_UMFPACK) {
        // Need to sort the matrix A for UMFPACK to work
        dCSRmat Acsr = fasp_format_dbsr_dcsr(&A);
        dCSRmat A_trans = fasp_dcsr_create(Acsr.row, Acsr.col, Acsr.nnz);
        fasp_dcsr_transz(&Acsr, NULL, &A_trans);
        fasp_dcsr_sort(&A_trans);
        status = fasp_solver_umfpack(&A_trans, &b, &x, print_level);
        fasp_dcsr_free(&A_trans);
        fasp_dcsr_free(&Acsr);
        if (status >= 0) status = 1; // Direct solver returns 1
    }
#endif

#ifdef OCP_USE_PARDISO // use PARDISO directly
    else if (solver_type == SOLVER_PARDISO) {
        dCSRmat Acsr = fasp_format_dbsr_dcsr(&A);
        fasp_dcsr_sort(&Acsr);
        status = fasp_solver_pardiso(&Acsr, &b, &x, print_level);
        fasp_dcsr_free(&Acsr);
        if (status >= 0) status = 1; // Direct solver returns 1
    }
#endif

    else {
        printf("### ERROR: Wrong solver type %d!!!\n", solver_type);
        status = ERROR_SOLVER_TYPE;
    }

    if (print_level > PRINT_MIN) {
        if (status < 0) {
            cout << "\n### WARNING: Solver does not converge!\n" << endl;
        }
        else {
            cout << "\nSolver converges successfully!\n" << endl;
        }
    }

    if (output_type) fclose(stdout);

    return status;
}


void VectorFaspSolver::InitParam()
{
    // Input/output
    inParam.print_level = PRINT_MIN;
    inParam.output_type = 0;

    // Problem information
    inParam.solver_type  = SOLVER_VFGMRES;
    inParam.decoup_type  = 1;
    inParam.precond_type = 64;
    inParam.stop_type    = STOP_REL_RES;

    // Solver parameters
    inParam.itsolver_tol   = 1e-3;
    inParam.itsolver_maxit = 100;
    inParam.restart        = 30;

    // ILU method parameters
    inParam.ILU_type    = ILUk;
    inParam.ILU_lfil    = 0;
    inParam.ILU_droptol = 0.001;
    inParam.ILU_relax   = 0;
    inParam.ILU_permtol = 0.0;

    // Schwarz method parameters
    inParam.SWZ_mmsize    = 200;
    inParam.SWZ_maxlvl    = 2;
    inParam.SWZ_type      = 1;
    inParam.SWZ_blksolver = SOLVER_DEFAULT;

    // AMG method parameters
    inParam.AMG_type                = CLASSIC_AMG;
    inParam.AMG_levels              = 20;
    inParam.AMG_cycle_type          = V_CYCLE;
    inParam.AMG_smoother            = SMOOTHER_GS;
    inParam.AMG_smooth_order        = CF_ORDER;
    inParam.AMG_presmooth_iter      = 1;
    inParam.AMG_postsmooth_iter     = 1;
    inParam.AMG_relaxation          = 1.0;
    inParam.AMG_coarse_dof          = 500;
    inParam.AMG_coarse_solver       = 0;
    inParam.AMG_tol                 = 1e-6;
    inParam.AMG_maxit               = 1;
    inParam.AMG_ILU_levels          = 0;
    inParam.AMG_SWZ_levels          = 0;
    inParam.AMG_coarse_scaling      = OFF; // Require investigation --Chensong
    inParam.AMG_amli_degree         = 1;
    inParam.AMG_nl_amli_krylov_type = 2;

    // Classical AMG specific
    inParam.AMG_coarsening_type      = 1;
    inParam.AMG_interpolation_type   = 1;
    inParam.AMG_max_row_sum          = 0.9;
    inParam.AMG_strong_threshold     = 0.3;
    inParam.AMG_truncation_threshold = 0.2;
    inParam.AMG_aggressive_level     = 0;
    inParam.AMG_aggressive_path      = 1;

    // Aggregation AMG specific
    inParam.AMG_aggregation_type   = PAIRWISE;
    inParam.AMG_quality_bound      = 8.0;
    inParam.AMG_pair_number        = 2;
    inParam.AMG_strong_coupled     = 0.25;
    inParam.AMG_max_aggregation    = 9;
    inParam.AMG_tentative_smooth   = 0.67;
    inParam.AMG_smooth_filter      = ON;
    inParam.AMG_smooth_restriction = ON;
}


void VectorFaspSolver::Allocate(const OCPMatrix& mat)
{
    A = fasp_dbsr_create(mat.maxDim, mat.maxDim, mat.max_nnz, mat.nb, 0);
    Asc = fasp_dbsr_create(mat.maxDim, mat.maxDim, mat.max_nnz, mat.nb, 0);
    fsc = fasp_dvec_create(mat.maxDim * mat.nb);
    order = fasp_ivec_create(mat.maxDim);
    Dmat.resize(mat.maxDim * mat.nb * mat.nb);
}



#endif // OCP_USE_FASP

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/22/2021      Create file                          */
/*  Chensong Zhang      Jan/19/2022      Set FASP4BLKOIL as optional          */
/*  Li Zhao             Apr/04/2022      Set FASP4CUDA   as optional          */
/*----------------------------------------------------------------------------*/
