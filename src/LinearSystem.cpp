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


USI LinearSystem::Setup(const OCPModel& model, const string& dir, const string& file, const Domain& d, const USI& nb)
{
    solveDir  = dir;
    domain    = &d;
    nbs.push_back(nb);

    mat.Allocate(*domain, nb);
    SetupLinearSolver(model, file);

    return LS.size() - 1;
}


void LinearSystem::SetWorkLS(const USI& i) 
{ 
    wIndex = i; 
    mat.SetBlockDim(nbs[wIndex]);
};


/// Setup LinearSolver
void LinearSystem::SetupLinearSolver(const OCPModel& model,
                                     const string& lsFile)
{
    string lsMethod = lsFile;
    auto pos = lsMethod.find_last_of('.');
    if (pos != string::npos) {
        // given params
        lsMethod = lsMethod.substr(pos + 1, lsMethod.size());
    }
    transform(lsMethod.begin(), lsMethod.end(), lsMethod.begin(), ::tolower);

    if (false) {}
#ifdef WITH_PARDISO
#if    OCPFLOATTYPEWIDTH == 64
    else if (lsMethod == "pardiso") {
        if (mat.nb > 1)    LS.push_back(new VectorPardisoSolver(solveDir, lsFile, mat));
        else               LS.push_back(new PardisoSolver(solveDir, lsFile, mat));
        LStype.push_back(OCPLStype::pardiso);
    }
#endif // OCPFLOATTYPEWIDTH == 64
#endif // WITH_PARDISO
#ifdef WITH_SAMG
    else if (lsMethod == "samg") {
        if (mat.nb > 1)    LS.push_back(new VectorSamgSolver(solveDir, lsFile, mat, model));
        else               LS.push_back(new ScalarSamgSolver(solveDir, lsFile, mat, model));
        LStype.push_back(OCPLStype::samg);
    }
#endif // WITH_SAMG
#if WITH_FASP
    else if (lsMethod == "fasp") {
        // if (domain->global_numproc > 1)  OCP_ABORT("FASP is only available for single process now!");
        if (mat.nb > 1)    LS.push_back(new VectorFaspSolver(solveDir, lsFile, mat));
        else               LS.push_back(new ScalarFaspSolver(solveDir, lsFile, mat));
        LStype.push_back(OCPLStype::fasp);
    }
#endif // WITH_FASP
#ifdef WITH_PETSCSOLVER
    else if (lsMethod == "petsc") {
        if (mat.nb > 1)    LS.push_back(new VectorPetscSolver(solveDir, lsFile, mat, domain));
        else               LS.push_back(new ScalarPetscSolver(solveDir, lsFile, mat, domain));
        LStype.push_back(OCPLStype::petsc);
    }
#endif // WITH_PETSCSOLVER
    else {
        OCP_ABORT("Wrong Linear Solver Type " + lsMethod + " !");
    }
}


void LinearSystem::AssembleMatLinearSolver()
{ 
#ifdef DEBUG
    mat.CheckLinearSystem();
#endif
    LS[wIndex]->AssembleMat(mat, domain); 
}


OCP_INT LinearSystem::Solve()
{
    OCP_INT iters = LS[wIndex]->Solve();
    if (LStype[wIndex] == OCPLStype::fasp) {
        if (iters < 0) {
            iters = -LS[wIndex]->GetNumIters();
            OCP_WARNING("FASP failed! -- " + to_string(iters));
        }
    }

#ifdef DEBUG
    mat.CheckSolution();
#endif

    return iters;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Shizhe Li           Nov/22/2021      renamed to LinearSystem              */
/*----------------------------------------------------------------------------*/