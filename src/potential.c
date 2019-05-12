/*
  potential.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Arham Amouye Foumani

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

// all functions for loading potential files and preparing or managing data structures related to
// potentials come here.

#include "potential.h"
#include "base.h"
#include "cspline.h"
#include "array.h"
#include "list.h"

static void EAM_convert_r_to_r2(eam_t *EAMp, double *source, double *dest)
/* Consider two functions f1 and f2 with the relation f2(r^2)=f1(r)
 * between them. Given the array source[0..EAM.Nr-1] containing a table
 * of the function f1, i.e. source[i]=f1(r_i), with r_i=i*EAM.dr , this
 * routine returns an array dest[0..EAM.Nr2-1] that contains the
 * values of the function f2 at points j*EAM.dr2 . */
{
    double *sourceDD;
    int i;

    sourceDD = (double *)malloc(EAMp->Nr * sizeof(double));
    spline_prepare(EAMp->dr, source, EAMp->Nr, sourceDD);

    for (i=0; i < EAMp->Nr2; i++)
        dest[i] = spline_val(EAMp->dr, source, sourceDD, sqrt(i * EAMp->dr2));
    free(sourceDD);
}

void fmd_pot_setCutoffRadius(fmd_sys_t *sysp, double cutoff)
{
    sysp->cutoffRadius = cutoff;
}

void fmd_pot_eam_free(fmd_sys_t *sysp)
{
/*
    int i, j;

    for (i=0; i < sysp->EAM.elementsNo; i++)
    {
        free(sysp->EAM.elements[i].F);
        free(sysp->EAM.elements[i].rho);
        for (j=0; j<=i; j++)
            free(sysp->EAM.elements[i].phi[j]);
        free(sysp->EAM.elements[i].phi);
#ifdef USE_CSPLINE
        free(sysp->EAM.elements[i].F_DD);
        free(sysp->EAM.elements[i].rhoDD);
        for (j=0; j<=i; j++)
            free(sysp->EAM.elements[i].phiDD[j]);
        free(sysp->EAM.elements[i].phiDD);
#endif
    }
    free(sysp->EAM.elements);
*/
}

// the following function will be deleted finally.

/*void fmd_pot_eam_init(fmd_sys_t *sysp, char *filePath)
{
    FILE *fp;
    char str[1024];
    int i, j, k;
    double *tempArray;

    if (sysp->subDomain.myrank == MAINPROCESS(sysp->subDomain.numprocs))
    {
        fp = fopen(filePath, "r");
        handleFileOpenError(fp, filePath);
        for (i=0; i<3; i++)
            fgets(str, 1024, fp);
        fscanf(fp, "%d", &sysp->EAM.elementsNo);
        sysp->EAM.elements = (TElement *)malloc(sysp->EAM.elementsNo * sizeof(TElement));
        for (i=0; i < sysp->EAM.elementsNo; i++)
            fscanf(fp, "%s", sysp->EAM.elements[i].name);
        fscanf(fp, "%d%lf%d%lf%lf", &sysp->EAM.Nrho, &sysp->EAM.drho, &sysp->EAM.Nr, &sysp->EAM.dr, &sysp->EAM.cutoff);
        sysp->EAM.Nr2 = (sysp->EAM.Nr += 2);
        assert( (sysp->EAM.Nr-1) * sysp->EAM.dr > sysp->EAM.cutoff );
        sysp->EAM.dr2 = SQR((sysp->EAM.Nr-1) * sysp->EAM.dr) / (sysp->EAM.Nr2-1);

        tempArray = (double *)malloc(sysp->EAM.Nr * sizeof(double));
        for (i=0; i < sysp->EAM.elementsNo; i++)
        {
            fscanf(fp, "%s%lf%lf%s", str, &sysp->EAM.elements[i].mass,
                &sysp->EAM.elements[i].latticeParameter, str);
            sysp->EAM.elements[i].mass /= MD_MASS_UNIT;
            sysp->EAM.elements[i].F = (double *)malloc(sysp->EAM.Nrho * sizeof(double));
            for (j=0; j < sysp->EAM.Nrho; j++)
                fscanf(fp, "%lf", &sysp->EAM.elements[i].F[j]);

            sysp->EAM.elements[i].rho = (double *)malloc(sysp->EAM.Nr2 * sizeof(double));
            for (j=0; j < sysp->EAM.Nr-2; j++)  // read rho(r) values from file
                fscanf(fp, "%lf", &tempArray[j]);
            tempArray[sysp->EAM.Nr-1] = tempArray[sysp->EAM.Nr-2] = 0.;
            EAM_convert_r_to_r2(&sysp->EAM, tempArray, sysp->EAM.elements[i].rho);

            sysp->EAM.elements[i].phi = (double **)malloc(sysp->EAM.elementsNo * sizeof(double *));
#ifdef USE_CSPLINE
            sysp->EAM.elements[i].F_DD = (double *)malloc(sysp->EAM.Nrho * sizeof(double));
            sysp->EAM.elements[i].rhoDD = (double *)malloc(sysp->EAM.Nr2 * sizeof(double));
            spline_prepare(sysp->EAM.drho, sysp->EAM.elements[i].F, sysp->EAM.Nrho, sysp->EAM.elements[i].F_DD);
            spline_prepare(sysp->EAM.dr2, sysp->EAM.elements[i].rho, sysp->EAM.Nr2, sysp->EAM.elements[i].rhoDD);
            sysp->EAM.elements[i].phiDD = (double **)malloc(sysp->EAM.elementsNo * sizeof(double *));
#endif
        }

        for (i=0; i < sysp->EAM.elementsNo; i++)
            for (j=0; j<=i; j++)
            {
                for (k=0; k < sysp->EAM.Nr-2; k++) // read r*phi values from file
                {
                    fscanf(fp, "%lf", &tempArray[k]);
                    if (k==0)
                        tempArray[k] = FLT_MAX;
                    else
                        tempArray[k] /= k * sysp->EAM.dr;
                }
                tempArray[sysp->EAM.Nr-1] = tempArray[sysp->EAM.Nr-2] = 0.;
                sysp->EAM.elements[i].phi[j] = (double *)malloc(sysp->EAM.Nr2 * sizeof(double));
                EAM_convert_r_to_r2(&sysp->EAM, tempArray, sysp->EAM.elements[i].phi[j]);
                sysp->EAM.elements[j].phi[i] = sysp->EAM.elements[i].phi[j];
#ifdef USE_CSPLINE
                sysp->EAM.elements[i].phiDD[j] = (double *)malloc(sysp->EAM.Nr2 * sizeof(double));
                spline_prepare(sysp->EAM.dr2, sysp->EAM.elements[i].phi[j], sysp->EAM.Nr2, sysp->EAM.elements[i].phiDD[j]);
                sysp->EAM.elements[j].phiDD[i] = sysp->EAM.elements[i].phiDD[j];
#endif
            }
        free(tempArray);
        fclose(fp);
    }

    MPI_Bcast(&sysp->EAM, sizeof(eam_t), MPI_CHAR, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
    if (sysp->subDomain.myrank != MAINPROCESS(sysp->subDomain.numprocs))
    {
        sysp->EAM.elements = (TElement *)malloc(sysp->EAM.elementsNo * sizeof(TElement));
        for (i=0; i < sysp->EAM.elementsNo; i++)
        {
            sysp->EAM.elements[i].F = (double *)malloc(sysp->EAM.Nrho * sizeof(double));
            sysp->EAM.elements[i].rho = (double *)malloc(sysp->EAM.Nr2 * sizeof(double));
            sysp->EAM.elements[i].phi = (double **)malloc(sysp->EAM.elementsNo * sizeof(double *));
            for (j=0; j<=i; j++)
            {
                sysp->EAM.elements[i].phi[j] = (double *)malloc(sysp->EAM.Nr2 * sizeof(double));
                sysp->EAM.elements[j].phi[i] = sysp->EAM.elements[i].phi[j];
            }
#ifdef USE_CSPLINE
            sysp->EAM.elements[i].F_DD = (double *)malloc(sysp->EAM.Nrho * sizeof(double));
            sysp->EAM.elements[i].rhoDD = (double *)malloc(sysp->EAM.Nr2 * sizeof(double));
            sysp->EAM.elements[i].phiDD = (double **)malloc(sysp->EAM.elementsNo * sizeof(double *));
            for (j=0; j<=i; j++)
            {
                sysp->EAM.elements[i].phiDD[j] = (double *)malloc(sysp->EAM.Nr2 * sizeof(double));
                sysp->EAM.elements[j].phiDD[i] = sysp->EAM.elements[i].phiDD[j];
            }
#endif
        }
    }

    for (i=0; i < sysp->EAM.elementsNo; i++)
    {
        MPI_Bcast(&sysp->EAM.elements[i].mass, 1, MPI_DOUBLE, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        MPI_Bcast(&sysp->EAM.elements[i].latticeParameter, 1, MPI_DOUBLE, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        MPI_Bcast(sysp->EAM.elements[i].name, 3, MPI_CHAR, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        MPI_Bcast(sysp->EAM.elements[i].F, sysp->EAM.Nrho, MPI_DOUBLE, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        MPI_Bcast(sysp->EAM.elements[i].rho, sysp->EAM.Nr2, MPI_DOUBLE, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        for (j=0; j<=i; j++)
            MPI_Bcast(sysp->EAM.elements[i].phi[j], sysp->EAM.Nr2, MPI_DOUBLE,
                      MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
#ifdef USE_CSPLINE
        MPI_Bcast(sysp->EAM.elements[i].F_DD, sysp->EAM.Nrho, MPI_DOUBLE, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        MPI_Bcast(sysp->EAM.elements[i].rhoDD, sysp->EAM.Nr2, MPI_DOUBLE, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        for (j=0; j<=i; j++)
            MPI_Bcast(sysp->EAM.elements[i].phiDD[j], sysp->EAM.Nr2, MPI_DOUBLE,
                      MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
#endif
    }
}*/

static eam_t *load_DYNAMOsetfl(char *filePath)
{
    FILE *fp;
    fp = fopen(filePath, "r");
    handleFileOpenError(fp, filePath);

    eam_t *eam;
    eam = (eam_t *)malloc(sizeof(eam_t));

    char str[1024];
    int i, j, k;
    for (i=0; i<3; i++)
        fgets(str, 1024, fp);
    fscanf(fp, "%d", &eam->elementsNo);
    eam->elements = (eam_element_t *)malloc(eam->elementsNo * sizeof(eam_element_t));
    for (i=0; i < eam->elementsNo; i++)
        fscanf(fp, "%s", eam->elements[i].name);
    fscanf(fp, "%d%lf%d%lf%lf", &eam->Nrho, &eam->drho, &eam->Nr, &eam->dr, &eam->cutoff);
    eam->Nr2 = (eam->Nr += 2);
    assert( (eam->Nr-1) * eam->dr > eam->cutoff );
    eam->dr2 = SQR((eam->Nr-1) * eam->dr) / (eam->Nr2-1);

    double *tempArray;
    tempArray = (double *)malloc(eam->Nr * sizeof(double));
    for (i=0; i < eam->elementsNo; i++)
    {
        fscanf(fp, "%s%lf%lf%s", str, &eam->elements[i].mass,
            &eam->elements[i].latticeParameter, str);
        eam->elements[i].mass /= MD_MASS_UNIT;
        eam->elements[i].F = (double *)malloc(eam->Nrho * sizeof(double));
        for (j=0; j < eam->Nrho; j++)
            fscanf(fp, "%lf", &eam->elements[i].F[j]);

        eam->elements[i].rho = (double *)malloc(eam->Nr2 * sizeof(double));
        for (j=0; j < eam->Nr-2; j++)  // read rho(r) values from file
            fscanf(fp, "%lf", &tempArray[j]);
        tempArray[eam->Nr-1] = tempArray[eam->Nr-2] = 0.;
        EAM_convert_r_to_r2(eam, tempArray, eam->elements[i].rho);

        eam->elements[i].phi = (double **)malloc(eam->elementsNo * sizeof(double *));
#ifdef USE_CSPLINE
        eam->elements[i].F_DD = (double *)malloc(eam->Nrho * sizeof(double));
        eam->elements[i].rhoDD = (double *)malloc(eam->Nr2 * sizeof(double));
        spline_prepare(eam->drho, eam->elements[i].F, eam->Nrho, eam->elements[i].F_DD);
        spline_prepare(eam->dr2, eam->elements[i].rho, eam->Nr2, eam->elements[i].rhoDD);
        eam->elements[i].phiDD = (double **)malloc(eam->elementsNo * sizeof(double *));
#endif
    }

    for (i=0; i < eam->elementsNo; i++)
        for (j=0; j<=i; j++)
        {
            for (k=0; k < eam->Nr-2; k++) // read r*phi values from file
            {
                fscanf(fp, "%lf", &tempArray[k]);
                if (k==0)
                    tempArray[k] = FLT_MAX;
                else
                    tempArray[k] /= k * eam->dr;
            }
            tempArray[eam->Nr-1] = tempArray[eam->Nr-2] = 0.;
            eam->elements[i].phi[j] = (double *)malloc(eam->Nr2 * sizeof(double));
            EAM_convert_r_to_r2(eam, tempArray, eam->elements[i].phi[j]);
            eam->elements[j].phi[i] = eam->elements[i].phi[j];
#ifdef USE_CSPLINE
            eam->elements[i].phiDD[j] = (double *)malloc(eam->Nr2 * sizeof(double));
            spline_prepare(eam->dr2, eam->elements[i].phi[j], eam->Nr2, eam->elements[i].phiDD[j]);
            eam->elements[j].phiDD[i] = eam->elements[i].phiDD[j];
#endif
        }
    free(tempArray);
    fclose(fp);

    return eam;
}

fmd_pot_t *fmd_pot_eam_alloy_load(fmd_sys_t *sysp, char *filePath)
{
    if (sysp->subDomain.myrank != MAINPROCESS(sysp->subDomain.numprocs))
        return;

    eam_t *eam = load_DYNAMOsetfl(filePath);

    fmd_pot_t *pot = (fmd_pot_t *)malloc(sizeof(fmd_pot_t));
    pot->kind = POTKIND_EAM_ALLOY;
    pot->data = eam;

    sysp->potsys.potlist = fmd_list_prepend(sysp->potsys.potlist, pot);

    return pot;
}

double fmd_pot_eam_getCutoffRadius(fmd_sys_t *sysp)
{
    //return sysp->EAM.cutoff;
}

double fmd_pot_eam_getLatticeParameter(fmd_sys_t *sysp, int element)
{
    //return sysp->EAM.elements[element].latticeParameter;
}

void fmd_pot_setAtomKinds(fmd_sys_t *sysp, unsigned number, fmd_atomkind_name_t *names, double *masses)
{
    sysp->potsys.atomkinds_num = number;
    sysp->potsys.atomkinds = (atomkind_t *)malloc(number * sizeof(atomkind_t));
    for (int i=0; i<number; i++)
    {
        sysp->potsys.atomkinds[i].mass = masses[i];
        strcpy(sysp->potsys.atomkinds[i].name, names[i]);
    }
}

void fmd_pot_free(fmd_sys_t *sysp)
{
    if (sysp->potsys.atomkinds != NULL)
    {
        free(sysp->potsys.atomkinds);
        sysp->potsys.atomkinds = NULL;
    }
    if (sysp->potsys.lj_6_12 != NULL)
    {
        fmd_array_neat2d_free((void **)sysp->potsys.lj_6_12);
        sysp->potsys.lj_6_12 = NULL;
    }
}

void fmd_pot_init(fmd_sys_t *sysp)
{
    sysp->potsys.atomkinds = NULL;
    sysp->potsys.lj_6_12 = NULL;
    sysp->potsys.potlist = NULL;
    sysp->potsys.pottable = NULL;
}

static void create_pottable(fmd_sys_t *sysp)
{
    sysp->potsys.pottable = (fmd_pot_t **)fmd_array_neat2d_create(sysp->potsys.atomkinds_num,
                                                                  sysp->potsys.atomkinds_num,
                                                                  sizeof(fmd_pot_t));
}

void fmd_pot_apply(fmd_sys_t *sysp, unsigned atomkind1, unsigned atomkind2, fmd_pot_t *pot)
{
    if (sysp->subDomain.myrank != MAINPROCESS(sysp->subDomain.numprocs))
        return;

    // create the pottable if doesn't exist
    if (sysp->potsys.pottable == NULL) create_pottable(sysp);

    // apply the pot to pottable
    sysp->potsys.pottable[atomkind1][atomkind2] =
      sysp->potsys.pottable[atomkind2][atomkind1] = (*pot);
}

fmd_pot_t *fmd_pot_lj_apply(fmd_sys_t *sysp, unsigned atomkind1, unsigned atomkind2,
                            double sigma, double epsilon, double cutoff)
{
    if (sysp->subDomain.myrank != MAINPROCESS(sysp->subDomain.numprocs))
        return;

    LJ_6_12_t *lj = (LJ_6_12_t *)malloc(sizeof(LJ_6_12_t));
    lj->sig = sigma;
    lj->eps = epsilon;
    lj->cutoff_sqr = SQR(cutoff);

    fmd_pot_t *pot = (fmd_pot_t *)malloc(sizeof(fmd_pot_t));
    pot->kind = POTKIND_LJ_6_12;
    pot->data = lj;

    // add the pot to potlist
    sysp->potsys.potlist = fmd_list_prepend(sysp->potsys.potlist, pot);

    // apply the pot
    fmd_pot_apply(sysp, atomkind1, atomkind2, pot);

    return pot;
}
