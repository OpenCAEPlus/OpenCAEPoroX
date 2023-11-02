/*! \file    LinearSystem.cpp
 *  \brief   Contains Internal Matrix structure, ptr to external linearsolver and
 *interface \author  Shizhe Li \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "LinearSystem.hpp"

void LinearSystem::AllocateRowMem(const USI& nb)
{
    blockSize = nb * nb;
    blockDim  = nb;
    maxDim    = domain->numElementLocal;
    dim       = 0;
    colId.resize(maxDim);
    val.resize(maxDim);
    b.resize((maxDim + domain->numGridGhost) * blockDim);
    u.resize((maxDim + domain->numGridGhost) * blockDim);  // ghost is included
}

void LinearSystem::AllocateColMem()
{
 
    USI maxPerfNum = 0;
    for (const auto& w : domain->wellWPB) {
        maxPerfNum = OCP_MAX(maxPerfNum, (w.size() - 1) / 2);
    }
    maxPerfNum++; // add well self

    max_nnz = 0;
    // bulk
    for (OCP_USI n = 0; n < domain->numGridInterior; n++) {
        max_nnz += domain->neighborNum[n];
        colId[n].reserve(domain->neighborNum[n]);
        val[n].reserve(domain->neighborNum[n] * blockSize);
    }
    // well
    for (OCP_USI n = domain->numGridInterior; n < maxDim; n++) {
        max_nnz += maxPerfNum;
        colId[n].reserve(maxPerfNum);
        val[n].reserve(maxPerfNum * blockSize);
    }
}

void LinearSystem::ClearData()
{
    for (OCP_USI i = 0; i < maxDim; i++) {
        colId[i].clear(); // actually, only parts of bulks needs to be clear
        val[i].clear();
    }
    // diagPtr.assign(maxDim, 0);
    fill(b.begin(), b.end(), 0.0);
    dim = 0;
    // In fact, for linear system the current solution is a good initial solution for
    // next step, so u will not be set to zero. u.assign(maxDim, 0);
}

void LinearSystem::AssembleRhsAccumulate(const vector<OCP_DBL>& rhs)
{
    OCP_USI nrow = dim * blockDim;
    for (OCP_USI i = 0; i < nrow; i++) {
        b[i] += rhs[i];
    }
}

void LinearSystem::AssembleRhsCopy(const vector<OCP_DBL>& rhs)
{
    copy(rhs.begin(), rhs.end(), b.data());
}

void LinearSystem::OutputLinearSystem(const string& fileA, const string& fileb) const
{
    string FileA = solveDir + fileA;
    string Fileb = solveDir + fileb;

    // out A
    // csr or bsr
    ofstream outA(FileA);
    if (!outA.is_open()) cout << "Can not open " << FileA << endl;

    ios::sync_with_stdio(false);
    outA.tie(0);

    outA << dim << "\n";
    if (blockDim != 1) {
        outA << blockDim << "\n";
    }
    // IA
    OCP_USI rowId = 0;
    for (OCP_USI i = 0; i < dim; i++) {
        outA << rowId << "\n";
        rowId += colId[i].size();
    }
    outA << rowId << "\n";
    // JA
    USI rowSize = 0;
    for (OCP_USI i = 0; i < dim; i++) {
        rowSize = colId[i].size();
        for (USI j = 0; j < rowSize; j++) {
            // outA << domain->global_index[colId[i][j]] << "\n";
            outA << colId[i][j] << "\n";
        }
    }
    // val
    for (OCP_USI i = 0; i < dim; i++) {
        rowSize = val[i].size();
        for (USI j = 0; j < rowSize; j++) {
            outA << val[i][j] << "\n";
        }
    }
    outA.close();

    // out b
    OCP_USI  nRow = dim * blockDim;
    ofstream outb(Fileb);
    if (!outb.is_open()) cout << "Can not open " << Fileb << endl;

    ios::sync_with_stdio(false);
    outb.tie(0);

    outb << dim << "\n";
    for (OCP_USI i = 0; i < nRow; i++) {
        outb << b[i] << "\n";
    }
    outb.close();
}

void LinearSystem::OutputSolution(const string& fileU) const
{
    string   FileU = solveDir + fileU;
    ofstream outu(FileU);
    if (!outu.is_open()) cout << "Can not open " << FileU << endl;
    outu << dim << "\n";
    OCP_USI nrow = dim * blockDim;
    for (OCP_USI i = 0; i < nrow; i++) outu << u[i] << "\n";
    outu.close();
}

void LinearSystem::CheckEquation() const
{
    // check A
    for (OCP_USI n = 0; n < dim; n++) {
        for (auto v : val[n]) {
            if (!isfinite(v)) {
                OCP_ABORT("NAN or INF in MAT");
            }
        }
    }
    // check b
    OCP_USI len = dim * blockDim;
    for (OCP_USI n = 0; n < len; n++) {
        if (!isfinite(b[n])) {
            OCP_ABORT("NAN or INF in rhs");
        }
    }
}

void LinearSystem::CheckSolution() const
{
    OCP_USI len = dim * blockDim;
    for (OCP_USI n = 0; n < len; n++) {
        if (!isfinite(u[n])) {
            OCP_ABORT("NAN or INF in u");
        }
    }
}

/// Setup LinearSolver
void LinearSystem::SetupLinearSolver(const OCPModel& model,
                                     const string&   dir,
                                     const string&   file)
{
    solveDir = dir;
    string lsMethod = file;
    auto pos = lsMethod.find_last_of('.');
    if (pos != string::npos) {
        // given params
        lsMethod = lsMethod.substr(pos + 1, lsMethod.size());
    }
    transform(lsMethod.begin(), lsMethod.end(), lsMethod.begin(), ::tolower);

    if (false) {}
#ifdef WITH_PARDISO
#if OCPFLOATTYPEWIDTH == 64
    else if (lsMethod == "pardiso") {
        if (blockDim > 1)    LS = new VectorPardisoSolver(blockDim);
        else                 LS = new PardisoSolver(blockDim);
    }
#endif
#endif // WITH_PARDISO
#ifdef WITH_SAMG
    else if (lsMethod == "samg") {
        if (blockDim > 1)    LS = new VectorSamgSolver(blockDim, model);
        else                 LS = new ScalarSamgSolver(model);
    }
#endif // WITH_SAMG
#if WITH_FASP
    else if (lsMethod == "fasp") {
        if (domain->numproc > 1)  OCP_ABORT("FASP is only available for single process now!");
        if (blockDim > 1)    LS = new VectorFaspSolver(blockDim);
        else                 LS = new ScalarFaspSolver();
}
#endif // WITH_FASP
#ifdef WITH_PETSCSOLVER
    else if (lsMethod == "petsc") {
        if (blockDim > 1)    LS = new VectorPetscSolver(blockDim, domain);
        else                 LS = new ScalarPetscSolver();
    }
#endif // WITH_PETSCSOLVER
    else {
        OCP_ABORT("Wrong Linear Solver Type " + lsMethod + " !");
    }

    LS->SetupParam(dir, file);
    LS->Allocate(max_nnz, maxDim);
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Shizhe Li           Nov/22/2021      renamed to LinearSystem              */
/*----------------------------------------------------------------------------*/