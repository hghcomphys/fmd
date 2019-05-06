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

#include "base.h"
#include "cspline.h"

static void EAM_convert_r_to_r2(TEAM *EAMp, double *source, double *dest)
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

void fmd_pot_setCutoffRadius(TSystem *sysp, double cutoff)
{
    sysp->cutoffRadius = cutoff;
}

void fmd_pot_eam_free(TSystem *sysp)
{
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
}

void fmd_pot_eam_init(TSystem *sysp, char *filePath)
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

    MPI_Bcast(&sysp->EAM, sizeof(TEAM), MPI_CHAR, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
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
}

double fmd_pot_eam_getCutoffRadius(TSystem *sysp)
{
    return sysp->EAM.cutoff;
}

double fmd_pot_eam_getLatticeParameter(TSystem *sysp, int element)
{
    return sysp->EAM.elements[element].latticeParameter;
}
