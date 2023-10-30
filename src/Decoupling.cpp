/*! \file    Decoupling.cpp
 *  \brief   Decoupling methods for vector-value problems
 *  \author  Shizhe Li
 *  \date    Oct/15/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "DenseMat.hpp"
#include "FaspSolver.hpp"



// ABF decoupling method
static void decouple_abf(dBSRmat* A, OCP_DBL* diaginv, dBSRmat* B)
{
    // members of A
    const OCP_INT  ROW = A->ROW;
    const OCP_INT  NNZ = A->NNZ;
    const OCP_INT  nb  = A->nb;
    const OCP_INT  nb2 = nb * nb;
    const OCP_INT* IA  = A->IA;
    const OCP_INT* JA  = A->JA;
    OCP_DBL*      val = A->val;

    OCP_INT i, k, m;

    // Create a link to dBSRmat 'B'
    OCP_INT*  IAb  = B->IA;
    OCP_INT*  JAb  = B->JA;
    OCP_DBL* valb = B->val;
    memcpy(IAb, IA, (ROW + 1) * sizeof(OCP_INT));
    memcpy(JAb, JA, NNZ * sizeof(OCP_INT));

    for (i = 0; i < ROW; ++i) {
        // get the diagonal sub-blocks D and their inverse
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            if (JA[k] == i) {
                memcpy(diaginv + i * nb2, val + k * nb2, nb2 * sizeof(OCP_DBL));
                fasp_smat_inv(diaginv + i * nb2, nb);
            }
        }
        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop
}

// Analytical decoupling method
static void decouple_anl(dBSRmat* A, OCP_DBL* diaginv, dBSRmat* B)
{
    // members of A
    const OCP_INT  ROW = A->ROW;
    const OCP_INT  NNZ = A->NNZ;
    const OCP_INT  nb  = A->nb;
    const OCP_INT  nb2 = nb * nb;
    const OCP_INT* IA  = A->IA;
    const OCP_INT* JA  = A->JA;
    OCP_DBL*      val = A->val;

    OCP_INT i, k, m;

    // Create a dBSRmat 'B'
    OCP_INT*  IAb  = B->IA;
    OCP_INT*  JAb  = B->JA;
    OCP_DBL* valb = B->val;
    memcpy(IAb, IA, (ROW + 1) * sizeof(OCP_INT));
    memcpy(JAb, JA, NNZ * sizeof(OCP_INT));

    for (i = 0; i < ROW; ++i) {
        // form the diagonal sub-blocks for analytical decoupling
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            if (JA[k] == i) {
                m = k * nb2;
                fasp_smat_identity(diaginv + i * nb2, nb, nb2);
                for (int l = 0; l < nb - 1; l++)
                    diaginv[i * nb2 + 1 + l] = -val[m + 1 + l];
            }
        }

        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop
}

static void decouple_truetrans(dBSRmat* A, OCP_DBL* diaginv, dBSRmat* B)
{
    // members of A
    const OCP_INT  ROW = A->ROW;
    const OCP_INT  NNZ = A->NNZ;
    const OCP_INT  nb  = A->nb;
    const OCP_INT  nb2 = nb * nb;
    const OCP_INT* IA  = A->IA;
    const OCP_INT* JA  = A->JA;
    OCP_DBL*      val = A->val;

    OCP_INT i, k, l, m;

    // Create a dBSRmat 'B'
    OCP_INT*  IAb  = B->IA;
    OCP_INT*  JAb  = B->JA;
    OCP_DBL* valb = B->val;
    memcpy(IAb, IA, (ROW + 1) * sizeof(OCP_INT));
    memcpy(JAb, JA, NNZ * sizeof(OCP_INT));

    OCP_DBL *mat1, *mat2;
    mat1 = new OCP_DBL[nb2];
    mat2 = new OCP_DBL[nb2];

    // Dset0(nb2*ROW, diaginv);
    for (i = 0; i < ROW; ++i) {
        OCP_DBL Tt = 0.0;
        // get the diagonal sub-blocks
        // mat2 is OCP_TRUE-IMPES matrix.
        fasp_smat_identity(mat2, nb, nb2);
        // mat1 is the mobility ratio matrix
        fasp_smat_identity(mat1, nb, nb2);

        for (l = 1; l < nb; ++l) {
            mat2[l] = diaginv[i * nb2 + l];
            Tt += diaginv[i * nb2 + l] * diaginv[i * nb2 + l * nb];
        }
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            if (JA[k] == i) break;
        }
        for (l = 0; l < nb - 1; ++l) {
            if (val[k * nb2 + (l + 1) * nb] > 0)
                mat1[(l + 1) * nb] = -diaginv[i * nb2 + (l + 1) * nb] / Tt; // diag(l)
        }
        DaABpbC(nb, nb, nb, 1, mat1, mat2, 0, diaginv + i * nb2);

        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop

    delete mat1;
    delete mat2;
}

static void decouple_truetrans_alg(dBSRmat* A, OCP_DBL* diaginv, dBSRmat* B)
{
    // members of A
    const OCP_INT  ROW = A->ROW;
    const OCP_INT  NNZ = A->NNZ;
    const OCP_INT  nb  = A->nb;
    const OCP_INT  nb2 = nb * nb;
    const OCP_INT* IA  = A->IA;
    const OCP_INT* JA  = A->JA;
    OCP_DBL*      val = A->val;

    OCP_INT i, k, l, m;

    // Create a link to dBSRmat 'B'
    OCP_INT*  IAb  = B->IA;
    OCP_INT*  JAb  = B->JA;
    OCP_DBL* valb = B->val;
    memcpy(IAb, IA, (ROW + 1) * sizeof(OCP_INT));
    memcpy(JAb, JA, NNZ * sizeof(OCP_INT));

    OCP_DBL *mat1, *mat2;
    mat1 = new OCP_DBL[nb2];
    mat2 = new OCP_DBL[nb2];

    // Dset0(nb2*ROW, diaginv);
    for (i = 0; i < ROW; ++i) {
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            OCP_DBL Tt = 0.0;
            if (JA[k] == i) {
                m = k * nb2;
                fasp_smat_identity(mat2, nb, nb2);
                for (l = 1; l < nb; ++l) {
                    mat2[l] = -val[m + l];
                    if (val[m + l * nb] > 0) {
                        Tt += -val[m + l] * val[m + l * nb];
                    }
                }

                fasp_smat_identity(mat1, nb, nb2);
                for (l = 1; l < nb; ++l) {
                    if (val[m + l * nb] > 0) {
                        mat1[l * nb] = -val[m + l * nb] / Tt; // diag(l)
                    }
                }

                DaABpbC(nb, nb, nb, 1, mat1, mat2, 0, diaginv + i * nb2);
            }
        }

        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop

    delete mat1;
    delete mat2;
}

// Semi-analytical decoupling method
static void decouple_abftrue(dBSRmat* A, OCP_DBL* diaginv, dBSRmat* B)
{
    // members of A
    const OCP_INT  ROW = A->ROW;
    const OCP_INT  NNZ = A->NNZ;
    const OCP_INT  nb  = A->nb;
    const OCP_INT  nb2 = nb * nb;
    const OCP_INT* IA  = A->IA;
    const OCP_INT* JA  = A->JA;
    OCP_DBL*      val = A->val;

    OCP_INT i, k, m;

    // Create a dBSRmat 'B'
    OCP_INT*  IAb  = B->IA;
    OCP_INT*  JAb  = B->JA;
    OCP_DBL* valb = B->val;
    memcpy(IAb, IA, (ROW + 1) * sizeof(OCP_INT));
    memcpy(JAb, JA, NNZ * sizeof(OCP_INT));

    for (i = 0; i < ROW; ++i) {
        // get the diagonal sub-blocks
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            if (JA[k] == i) {
                // Form the ABF part first
                memcpy(diaginv + i * nb2, val + k * nb2, nb2 * sizeof(OCP_DBL));
                fasp_smat_inv(diaginv + i * nb2, nb);
                // Replace the first line with analytical decoupling
                m                = k * nb2;
                diaginv[i * nb2] = 1;
                for (int l = 0; l < nb - 1; l++)
                    diaginv[i * nb2 + 1 + l] = -val[m + 1 + l];
            }
        }

        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop
}

static void decouple_true_scale(dBSRmat* A, OCP_DBL* diaginv, dBSRmat* B)
{
    // members of A
    const OCP_INT  ROW = A->ROW;
    const OCP_INT  NNZ = A->NNZ;
    const OCP_INT  nb  = A->nb;
    const OCP_INT  nb2 = nb * nb;
    const OCP_INT* IA  = A->IA;
    const OCP_INT* JA  = A->JA;
    OCP_DBL*      val = A->val;

    OCP_INT i, k, m;

    // Create a dBSRmat 'B'
    OCP_INT*  IAb  = B->IA;
    OCP_INT*  JAb  = B->JA;
    OCP_DBL* valb = B->val;
    memcpy(IAb, IA, (ROW + 1) * sizeof(OCP_INT));
    memcpy(JAb, JA, NNZ * sizeof(OCP_INT));

    for (i = 0; i < ROW; ++i) {
        // get the diagonal sub-blocks
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            if (JA[k] == i) {
                m = k * nb2;
                fasp_smat_identity(diaginv + i * nb2, nb, nb2);
                diaginv[i * nb2] = 1 / val[m];
                for (int l = 0; l < nb - 1; l++)
                    diaginv[i * nb2 + 1 + l] = -val[m + 1 + l] / val[m];
            }
        }

        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop
}

static void decouple_rotate(dBSRmat* A, OCP_DBL* diaginv, dBSRmat* B)
{
    // members of A
    const OCP_INT  ROW = A->ROW;
    const OCP_INT  NNZ = A->NNZ;
    const OCP_INT  nb  = A->nb;
    const OCP_INT  nb2 = nb * nb;
    const OCP_INT* IA  = A->IA;
    const OCP_INT* JA  = A->JA;
    OCP_DBL*      val = A->val;

    OCP_INT i, k, m;

    // Create a dBSRmat 'B'
    OCP_INT*  IAb  = B->IA;
    OCP_INT*  JAb  = B->JA;
    OCP_DBL* valb = B->val;
    memcpy(IAb, IA, (ROW + 1) * sizeof(OCP_INT));
    memcpy(JAb, JA, NNZ * sizeof(OCP_INT));

    for (i = 0; i < ROW; ++i) {
        // get the diagonal sub-blocks
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            if (JA[k] == i) {
                m = k * nb2;
                for (int li = 0; li < nb; li++) {
                    for (int lj = 0; lj < nb; lj++) {
                        diaginv[i * nb2 + li * nb + lj] = 0;
                        if (lj - li == 1) {
                            diaginv[i * nb2 + li * nb + lj] = 1;
                        }
                        if (lj == 0 && li == nb - 1) {
                            diaginv[i * nb2 + li * nb + lj] = 1;
                        }
                    }
                }
            }
        }

        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop
}

// Quasi-IMPES decoupling method
static void decouple_quasi(dBSRmat* A, OCP_DBL* diaginv, dBSRmat* B)
{
    // members of A
    const OCP_INT  ROW = A->ROW;
    const OCP_INT  NNZ = A->NNZ;
    const OCP_INT  nb  = A->nb;
    const OCP_INT  nb2 = nb * nb;
    const OCP_INT* IA  = A->IA;
    const OCP_INT* JA  = A->JA;
    OCP_DBL*      val = A->val;

    OCP_INT i, k, m;

    OCP_INT*  IAb  = B->IA;
    OCP_INT*  JAb  = B->JA;
    OCP_DBL* valb = B->val;
    memcpy(IAb, IA, (ROW + 1) * sizeof(OCP_INT));
    memcpy(JAb, JA, NNZ * sizeof(OCP_INT));

    for (i = 0; i < ROW; ++i) {
        // get the diagonal sub-blocks
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            if (JA[k] == i) {
                m = k * nb2;
                fasp_smat_identity(diaginv + i * nb2, nb, nb2);
                for (int l = 0; l < nb - 1; l++)
                    diaginv[i * nb2 + 1 + l] =
                        -val[m + 1 + l] / val[m + (l + 1) * nb + l + 1];
            }
        }

        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop
}

// What's the difference with decouple_anl?
static void decouple_trueabf(dBSRmat* A, OCP_DBL* diaginv, dBSRmat* B)
{
    // members of A
    const OCP_INT  ROW = A->ROW;
    const OCP_INT  NNZ = A->NNZ;
    const OCP_INT  nb  = A->nb;
    const OCP_INT  nb2 = nb * nb;
    const OCP_INT* IA  = A->IA;
    const OCP_INT* JA  = A->JA;
    OCP_DBL*      val = A->val;

    OCP_INT i, k, m;

    OCP_INT*  IAb  = B->IA;
    OCP_INT*  JAb  = B->JA;
    OCP_DBL* valb = B->val;
    memcpy(IAb, IA, (ROW + 1) * sizeof(OCP_INT));
    memcpy(JAb, JA, NNZ * sizeof(OCP_INT));

    for (i = 0; i < ROW; ++i) {
        // get the diagonal sub-blocks
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            if (JA[k] == i) {
                m = k * nb2;
                fasp_smat_identity(diaginv + i * nb2, nb, nb2);
                for (int l = 0; l < nb; l++) diaginv[i * nb2 + l] = val[m + l];
                fasp_smat_inv(diaginv + i * nb2, nb);
            }
        }

        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop
}

/// Applying a decoupling algorithm for linear systems of FIM
void VectorFaspSolver::Decoupling(dBSRmat* Absr,
                                  dvector* b,
                                  dBSRmat* Asc,
                                  dvector* fsc,
                                  ivector* order,
                                  OCP_DBL*  Dmatvec,
                                  int      decoupleType)
{
    int              nrow = Absr->ROW;
    int              nb   = Absr->nb;
    OCP_DBL*         Dmat = Dmatvec;
    precond_diag_bsr diagA;

    // Natural ordering
    for (int i = 0; i < nrow; ++i) order->val[i] = i;

    // Without decoupling
    if (decoupleType == 0) {
        fasp_dbsr_cp(Absr, Asc); // Asc = Absr;
        fasp_dvec_cp(b, fsc);    // fsc = b;
        return;
    }

    // With decoupling
    switch (decoupleType) {
        case 2:
            decouple_anl(Absr, Dmat, Asc);
            break;
        case 3:
            decouple_quasi(Absr, Dmat, Asc);
            break;
        case 4:
            decouple_trueabf(Absr, Dmat, Asc);
            break;
        case 5:
            decouple_abftrue(Absr, Dmat, Asc);
            break;
        case 6:
            decouple_abftrue(Absr, Dmat, Asc);
            break;
        case 7:
            decouple_truetrans_alg(Absr, Dmat, Asc);
            break;
        case 8:
            decouple_truetrans(Absr, Dmat, Asc);
            break;
        case 9:
            decouple_true_scale(Absr, Dmat, Asc);
            break;
        case 10:
            decouple_rotate(Absr, Dmat, Asc);
            break;
        default: // case 1:
            decouple_abf(Absr, Dmat, Asc);
    }

    diagA.diag.row = nrow * nb * nb;
    diagA.diag.val = Dmat;
    diagA.nb       = nb;
    fasp_precond_dbsr_diag(b->val, fsc->val, &diagA);
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/15/2021      Create file                          */
/*  Chensong Zhang      Nov/09/2021      Restruct decoupling methods          */
/*  Chensong Zhang      Nov/30/2021      Add null decoupling                  */
/*----------------------------------------------------------------------------*/