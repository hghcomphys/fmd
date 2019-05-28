/*
  potential.h: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Arham Amouye Foumani, Hossein Ghorbanfekr Kalashami

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

#ifndef POTENTIAL_H
#define POTENTIAL_H

typedef char *fmd_string_t;

typedef struct
{
    double mass;
    double latticeParameter;
    double *F;
    double *F_DD;
    double *rho;
    double *rhoDD;
    double **phi;
    double **phiDD;
    fmd_string_t name;
} eam_element_t;

typedef struct
{
    eam_element_t *elements;
    double drho, dr, dr2, cutoff_sqr;
    int elementsNo;
    int Nrho, Nr, Nr2;
} eam_t;

typedef struct
{
    double eps;
    double sig;
    double cutoff_sqr;
} LJ_6_12_t;

typedef struct
{
    double D0;
    double alpha;
    double r0;
    double cutoff_sqr;
} morse_t;

typedef enum
{
    POTKIND_NONE,
    POTKIND_LJ_6_12,
    POTKIND_EAM_ALLOY,
    POTKIND_MORSE
} potkind_t;

typedef struct
{
    potkind_t kind;
    void *data;
} fmd_pot_t;

typedef struct
{
    potkind_t kind;
    void *data;
    unsigned iloc, jloc;    // local indexes inside the potential, used in potentials like EAM
} potpair_t;

typedef struct
{
    double eam_F0;
} atomkind_aux_t;

typedef struct
{
    double mass;
    fmd_string_t name;
    atomkind_aux_t *aux;
} atomkind_t;

typedef struct list_t list_t;

typedef struct
{
    unsigned atomkinds_num;
    atomkind_t *atomkinds;
    potpair_t **pottable;       // table of applied pots
    list_t *potkinds;           // list of pot kinds that are present in pottable
    unsigned potkinds_num;
    list_t *potlist;            // list of all pots, whether applied or not
    int hybridpasses[2];
} potsys_t;

typedef struct fmd_sys_t fmd_sys_t;

void fmd_pot_free(fmd_sys_t *sysp);
void fmd_pot_init(fmd_sys_t *sysp);
void fmd_pot_potkinds_update(fmd_sys_t *sysp);
void fmd_pot_hybridpasses_update(fmd_sys_t *sysp);

#endif /* POTENTIAL_H */
