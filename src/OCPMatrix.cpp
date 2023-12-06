/*! \file    OCPMatrix.cpp
 *  \brief   OCPMatrix class declaration
 *  \author  Shizhe Li
 *  \date    Dec/06/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPMatrix.hpp"


void OCPMatrix::Allocate(const Domain& domain, const USI& blockdim)
{
    SetBlockDim(blockdim);
    AllocateRowMem(domain);
    AllocateColMem(domain);
}


void OCPMatrix::ClearData()
{
    // actually, only parts of bulks needs to be clear
    for (OCP_USI i = 0; i < maxDim; i++) {
        colId[i].clear(); 
        val[i].clear();
    }
    // diagPtr.assign(maxDim, 0);
    fill(b.begin(), b.end(), 0.0);
    dim = 0;
    // In fact, for linear system the current solution is a good initial solution for
    // next step, so u will not be set to zero. u.assign(maxDim, 0);
}


void OCPMatrix::SetBlockDim(const USI& blockdim)
{
    nb = blockdim;
    nb2 = nb * nb;
}


void OCPMatrix::AllocateRowMem(const Domain& domain)
{
    maxDim = domain.numElementLocal;
    dim = 0;
    colId.resize(maxDim);
    val.resize(maxDim);
    // ghost is included
    b.resize((maxDim + domain.numGridGhost) * nb);
    u.resize((maxDim + domain.numGridGhost) * nb);
}

void OCPMatrix::AllocateColMem(const Domain& domain)
{

    USI maxPerfNum = 0;
    for (const auto& w : domain.wellWPB) {
        maxPerfNum = OCP_MAX(maxPerfNum, (w.size() - 1) / 2);
    }
    maxPerfNum++; // add well self

    max_nnz = 0;
    // bulk
    for (OCP_USI n = 0; n < domain.numGridInterior; n++) {
        max_nnz += domain.neighborNum[n];
        colId[n].reserve(domain.neighborNum[n]);
        val[n].reserve(domain.neighborNum[n] * nb2);
    }
    // well
    for (OCP_USI n = domain.numGridInterior; n < maxDim; n++) {
        max_nnz += maxPerfNum;
        colId[n].reserve(maxPerfNum);
        val[n].reserve(maxPerfNum * nb2);
    }
}



OCP_USI OCPMatrix::AddDim(const OCP_USI& n)
{
    dim += n;
    return dim;
}


void OCPMatrix::NewDiag(const OCP_USI& n, const OCP_DBL& v)
{
    OCP_ASSERT(colId[n].size() == 0, "Wrong Diag");
    colId[n].push_back(n);
    val[n].push_back(v);
}


void OCPMatrix::NewDiag(const OCP_USI& n, const vector<OCP_DBL>& v)
{
    OCP_ASSERT(colId[n].size() == 0, "Wrong Diag");
    colId[n].push_back(n);
    val[n].insert(val[n].begin(), v.begin(), v.end());
}


void OCPMatrix::AddDiag(const OCP_USI& n, const OCP_DBL& v)
{
    OCP_ASSERT(colId[n].size() > 0, "Wrong Diag");
    val[n][0] += v;
}


void OCPMatrix::AddDiag(const OCP_USI& n, const vector<OCP_DBL>& v)
{
    OCP_ASSERT(colId[n].size() > 0, "Wrong Diag");
    for (USI i = 0; i < nb2; i++) {
        val[n][i] += v[i];
    }
}


void OCPMatrix::NewOffDiag(const OCP_USI& bId, const OCP_USI& eId, const OCP_DBL& v)
{
    OCP_ASSERT(colId[bId].size() > 0, "Wrong Diag");
    colId[bId].push_back(eId);
    val[bId].push_back(v);
}


void OCPMatrix::NewOffDiag(const OCP_USI& bId, const OCP_USI& eId, const vector<OCP_DBL>& v)
{
    OCP_ASSERT(colId[bId].size() > 0, "Wrong Diag");
    colId[bId].push_back(eId);
    val[bId].insert(val[bId].end(), v.begin(), v.end());
}


void OCPMatrix::AddRhs(const OCP_USI& n, const OCP_DBL& v) 
{
    b[n] += v;
}


void OCPMatrix::AddRhs(const OCP_USI& n, const vector<OCP_DBL>& v)
{
    for (USI i = 0; i < nb; i++) {
        b[n * nb + i] += v[i];
    }
}


void OCPMatrix::CopyRhs(const vector<OCP_DBL>& rhs)
{
    copy(rhs.begin(), rhs.end(), b.data());
}


void OCPMatrix::OutputLinearSystem(const string& dir, const string& fileA, const string& fileb) const
{
    string FileA = dir + fileA;
    string Fileb = dir + fileb;

    // out A
    // csr or bsr
    ofstream outA(FileA);
    if (!outA.is_open()) cout << "Can not open " << FileA << endl;

    ios::sync_with_stdio(false);
    outA.tie(0);

    outA << dim << "\n";
    if (nb != 1) {
        outA << nb << "\n";
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
    OCP_USI  nRow = dim * nb;
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


void OCPMatrix::OutputSolution(const string& dir, const string& fileU) const
{
    string   FileU = dir + fileU;
    ofstream outu(FileU);
    if (!outu.is_open()) cout << "Can not open " << FileU << endl;
    outu << dim << "\n";
    OCP_USI nrow = dim * nb;
    for (OCP_USI i = 0; i < nrow; i++) outu << u[i] << "\n";
    outu.close();
}


void OCPMatrix::CheckLinearSystem() const
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
    OCP_USI len = dim * nb;
    for (OCP_USI n = 0; n < len; n++) {
        if (!isfinite(b[n])) {
            OCP_ABORT("NAN or INF in rhs");
        }
    }
}

void OCPMatrix::CheckSolution() const
{
    OCP_USI len = dim * nb;
    for (OCP_USI n = 0; n < len; n++) {
        if (!isfinite(u[n])) {
            OCP_ABORT("NAN or INF in u");
        }
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/06/2023      Create file                          */
/*----------------------------------------------------------------------------*/