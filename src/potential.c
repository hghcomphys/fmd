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

void fmd_pot_setCutoffRadius(fmd_sys_t *sysp, double cutoff)
{
    sysp->cutoffRadius = cutoff;
}

void fmd_pot_eam_free(fmd_sys_t *sysp)
{
/*
    int i, j;

    for (i=0; i < eam->elementsNo; i++)
    {
        free(eam->elements[i].F);
        free(eam->elements[i].rho);
        for (j=0; j<=i; j++)
            free(eam->elements[i].phi[j]);
        free(eam->elements[i].phi);
#ifdef USE_CSPLINE
        free(eam->elements[i].F_DD);
        free(eam->elements[i].rhoDD);
        for (j=0; j<=i; j++)
            free(eam->elements[i].phiDD[j]);
        free(eam->elements[i].phiDD);
#endif
    }
    free(eam->elements);
*/
}

static void EAM_convert_r_to_r2(eam_t *eam, double *source, double *dest)
/* Consider two functions f1 and f2 with the relation f2(r^2)=f1(r)
 * between them. Given the array source[0..eam.Nr-1] containing a table
 * of the function f1, i.e. source[i]=f1(r_i), with r_i=i*eam.dr , this
 * routine returns an array dest[0..eam.Nr2-1] that contains the
 * values of the function f2 at points j*eam.dr2 . */
{
    double *sourceDD;
    int i;

    sourceDD = (double *)malloc(eam->Nr * sizeof(double));
    spline_prepare(eam->dr, source, eam->Nr, sourceDD);

    for (i=0; i < eam->Nr2; i++)
        dest[i] = spline_val(eam->dr, source, sourceDD, sqrt(i * eam->dr2));
    free(sourceDD);
}

static eam_t *load_DYNAMOsetfl(fmd_sys_t *sysp, char *filePath)
{
    int i, j, k;

    eam_t *eam = (eam_t *)malloc(sizeof(eam_t));

    if (sysp->subDomain.myrank == MAINPROCESS(sysp->subDomain.numprocs))
    {
        FILE *fp = fopen(filePath, "r");
        handleFileOpenError(fp, filePath);

        char str[1024];
        for (i=0; i<3; i++)
            fgets(str, 1024, fp);

        fscanf(fp, "%d", &eam->elementsNo);
        eam->elements = (eam_element_t *)malloc(eam->elementsNo * sizeof(eam_element_t));
        for (i=0; i < eam->elementsNo; i++)
            fscanf(fp, "%s", eam->elements[i].name);

        double cutoff;
        fscanf(fp, "%d%lf%d%lf%lf", &eam->Nrho, &eam->drho, &eam->Nr, &eam->dr, &cutoff);
        eam->Nr2 = (eam->Nr += 2);
        assert( (eam->Nr-1) * eam->dr > cutoff );
        eam->cutoff_sqr = SQR(cutoff);
        eam->dr2 = SQR((eam->Nr-1) * eam->dr) / (eam->Nr2-1);

        double *tempArray = (double *)malloc(eam->Nr * sizeof(double));
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
    }

    MPI_Bcast(eam, sizeof(eam_t), MPI_CHAR, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
    if (sysp->subDomain.myrank != MAINPROCESS(sysp->subDomain.numprocs))
    {
        eam->elements = (eam_element_t *)malloc(eam->elementsNo * sizeof(eam_element_t));
        for (i=0; i < eam->elementsNo; i++)
        {
            eam->elements[i].F = (double *)malloc(eam->Nrho * sizeof(double));
            eam->elements[i].rho = (double *)malloc(eam->Nr2 * sizeof(double));
            eam->elements[i].phi = (double **)malloc(eam->elementsNo * sizeof(double *));
            for (j=0; j<=i; j++)
            {
                eam->elements[i].phi[j] = (double *)malloc(eam->Nr2 * sizeof(double));
                eam->elements[j].phi[i] = eam->elements[i].phi[j];
            }
#ifdef USE_CSPLINE
            eam->elements[i].F_DD = (double *)malloc(eam->Nrho * sizeof(double));
            eam->elements[i].rhoDD = (double *)malloc(eam->Nr2 * sizeof(double));
            eam->elements[i].phiDD = (double **)malloc(eam->elementsNo * sizeof(double *));
            for (j=0; j<=i; j++)
            {
                eam->elements[i].phiDD[j] = (double *)malloc(eam->Nr2 * sizeof(double));
                eam->elements[j].phiDD[i] = eam->elements[i].phiDD[j];
            }
#endif
        }
    }

    for (i=0; i < eam->elementsNo; i++)
    {
        MPI_Bcast(&eam->elements[i].mass, 1, MPI_DOUBLE, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        MPI_Bcast(&eam->elements[i].latticeParameter, 1, MPI_DOUBLE, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        MPI_Bcast(eam->elements[i].name, 3, MPI_CHAR, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        MPI_Bcast(eam->elements[i].F, eam->Nrho, MPI_DOUBLE, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        MPI_Bcast(eam->elements[i].rho, eam->Nr2, MPI_DOUBLE, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        for (j=0; j<=i; j++)
            MPI_Bcast(eam->elements[i].phi[j], eam->Nr2, MPI_DOUBLE,
                      MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
#ifdef USE_CSPLINE
        MPI_Bcast(eam->elements[i].F_DD, eam->Nrho, MPI_DOUBLE, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        MPI_Bcast(eam->elements[i].rhoDD, eam->Nr2, MPI_DOUBLE, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        for (j=0; j<=i; j++)
            MPI_Bcast(eam->elements[i].phiDD[j], eam->Nr2, MPI_DOUBLE,
                      MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
#endif
    }

    return eam;
}

fmd_pot_t *fmd_pot_eam_alloy_load(fmd_sys_t *sysp, char *filePath)
{
    eam_t *eam = load_DYNAMOsetfl(sysp, filePath);

    fmd_pot_t *pot = (fmd_pot_t *)malloc(sizeof(fmd_pot_t));
    pot->kind = POTKIND_EAM_ALLOY;
    pot->data = eam;

    sysp->potsys.potlist = fmd_list_prepend(sysp->potsys.potlist, pot);

    return pot;
}

double fmd_pot_eam_getCutoffRadius(fmd_sys_t *sysp, fmd_pot_t *pot)
{
    // TO-DO: handle error
    assert(pot->kind == POTKIND_EAM_ALLOY);

    return sqrt(((eam_t *)pot->data)->cutoff_sqr);
}

double fmd_pot_eam_getLatticeParameter(fmd_sys_t *sysp, int element)
{
    //return eam->elements[element].latticeParameter;
}

void fmd_pot_setAtomKinds(fmd_sys_t *sysp, unsigned number, fmd_atomkind_name_t *names, double *masses)
{
    sysp->potsys.atomkinds_num = number;
    sysp->potsys.atomkinds = (atomkind_t *)malloc(number * sizeof(atomkind_t));
    for (int i=0; i<number; i++)
    {
        sysp->potsys.atomkinds[i].mass = masses[i];
        strcpy(sysp->potsys.atomkinds[i].name, names[i]);
        sysp->potsys.atomkinds[i].aux = (atomkind_aux_t *)malloc(sizeof(atomkind_aux_t));
        // TO-DO: initialize aux
    }
}

void fmd_pot_free(fmd_sys_t *sysp)
{
    if (sysp->potsys.atomkinds != NULL)
    {
        // TO-DO: free aux
        free(sysp->potsys.atomkinds);
        sysp->potsys.atomkinds = NULL;
    }

    if (sysp->potsys.pottable != NULL)
    {
        fmd_array_neat2d_free((void **)sysp->potsys.pottable);
        sysp->potsys.pottable = NULL;
    }

    if (sysp->potsys.potlist != NULL)
    {
        // TO-DO: FREE_POTLIST
        // sysp->potsys.potlist = NULL;
    }

    if (sysp->potsys.potkinds != NULL)
    {
        // TO-DO: FREE_POTKINDS
        // sysp->potsys.potkinds = NULL;
    }
}

void fmd_pot_init(fmd_sys_t *sysp)
{
    sysp->potsys.atomkinds = NULL;
    sysp->potsys.potlist = NULL;
    sysp->potsys.pottable = NULL;
    sysp->potsys.potkinds = NULL;
}

static void pottable_create(fmd_sys_t *sysp)
{
    sysp->potsys.pottable = (potpair_t **)fmd_array_neat2d_create(sysp->potsys.atomkinds_num,
                                                                  sysp->potsys.atomkinds_num,
                                                                  sizeof(potpair_t));
    for (unsigned i=0; i < sysp->potsys.atomkinds_num; i++)
        for (unsigned j=0; j <= i; j++)
            sysp->potsys.pottable[i][j].kind = POTKIND_NONE;
}

static unsigned pot_eam_find_iloc(fmd_sys_t *sysp, eam_t *eam, unsigned atomkind)
{
    for (int i=0; i < eam->elementsNo; i++)
        if (strcmp(sysp->potsys.atomkinds[atomkind].name, eam->elements[i].name) == 0)
            return i;

    // return -1 if the given eam potential doesn't include the specified atom kind
    return -1;
}

void fmd_pot_apply(fmd_sys_t *sysp, unsigned atomkind1, unsigned atomkind2, fmd_pot_t *pot)
{
    // create the pottable if doesn't exist
    if (sysp->potsys.pottable == NULL) pottable_create(sysp);

    // if potential is of eam alloy kind, find the local indexes
    if (pot->kind == POTKIND_EAM_ALLOY)
    {
        unsigned loc1, loc2;

        loc1 = pot_eam_find_iloc(sysp, pot->data, atomkind1);
        loc2 = pot_eam_find_iloc(sysp, pot->data, atomkind2);

        if (loc1 == -1 || loc2 == -1)
        {
            // error should be handled here
            assert(loc1 != -1 && loc2 != -1);
        }

        sysp->potsys.pottable[atomkind1][atomkind2].iloc = loc1;
        sysp->potsys.pottable[atomkind1][atomkind2].jloc = loc2;
        sysp->potsys.pottable[atomkind2][atomkind1].iloc = loc2;
        sysp->potsys.pottable[atomkind2][atomkind1].jloc = loc1;
    }

    //
    sysp->potsys.pottable[atomkind1][atomkind2].kind =
      sysp->potsys.pottable[atomkind2][atomkind1].kind = pot->kind;
    sysp->potsys.pottable[atomkind1][atomkind2].data =
      sysp->potsys.pottable[atomkind2][atomkind1].data = pot->data;
}

fmd_pot_t *fmd_pot_lj_apply(fmd_sys_t *sysp, unsigned atomkind1, unsigned atomkind2,
                            double sigma, double epsilon, double cutoff)
{
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

static int potkind_compare(const void *a, const void *b)
{
    if ( *( (potkind_t *)a ) == *( (potkind_t *)b ) )
        return 0;
    else
        return 1;
}

// TO-DO?: first, clean the list
void fmd_pot_potkinds_update(fmd_sys_t *sysp)
{
    // error should be handled here
    assert(sysp->potsys.pottable != NULL);

    sysp->potsys.potkinds_num = 0;

    for (unsigned i=0; i < sysp->potsys.atomkinds_num; i++)
        for (unsigned j=0; j <= i; j++)
        {
            potpair_t *potpair = &sysp->potsys.pottable[i][j];

            // error should be handled here
            assert(potpair->kind != POTKIND_NONE);

            // add the potkind to potkinds list, if isn't already included there
            if (fmd_list_find_custom(sysp->potsys.potkinds, &potpair->kind, potkind_compare) == NULL)
            {
                sysp->potsys.potkinds_num++;
                sysp->potsys.potkinds = fmd_list_prepend(sysp->potsys.potkinds, &potpair->kind);
            }
        }
}
