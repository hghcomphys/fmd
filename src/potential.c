/*
  potential.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Arham Amouye Foumani, Hossein Ghorbanfekr

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
#include "array.h"
#include "list.h"
#include "eam.h"

void fmd_pot_setCutoffRadius(fmd_sys_t *sysp, double cutoff)
{
    sysp->cutoffRadius = cutoff;
}

void fmd_pot_setAtomKinds(fmd_sys_t *sysp, unsigned number, const fmd_string_t names[], const double masses[])
{
    // TO-DO: error should be handled here
    assert(number > 0);

    sysp->potsys.atomkinds_num = number;
    sysp->potsys.atomkinds = (atomkind_t *)malloc(number * sizeof(atomkind_t));
    for (unsigned i=0; i<number; i++)
    {
        sysp->potsys.atomkinds[i].mass = masses[i] / MD_MASS_UNIT;  // convert from amu to internal mass unit
        size_t len = strlen(names[i]);
        sysp->potsys.atomkinds[i].name = (char *)malloc(len + 1);
        strcpy(sysp->potsys.atomkinds[i].name, names[i]);
    }
}

static void potlist_free(fmd_sys_t *sysp)
{
    list_t *potlist = sysp->potsys.potlist;

    while (potlist != NULL)
    {
        fmd_pot_t *pot = (fmd_pot_t *)(potlist->data);

        switch (pot->kind)
        {
            case POTKIND_EAM_ALLOY:
                fmd_pot_eam_free((eam_t *)(pot->data));
                break;

            default:
                free(pot->data);
        }

        potlist = potlist->next;
    }

    fmd_list_free(sysp->potsys.potlist);
}


void fmd_potsys_free(fmd_sys_t *sysp)
{
    if (sysp->potsys.atomkinds != NULL)
    {
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
        potlist_free(sysp);
        sysp->potsys.potlist = NULL;
    }

    if (sysp->potsys.potkinds != NULL)
    {
        fmd_list_free(sysp->potsys.potkinds);
        sysp->potsys.potkinds = NULL;
    }
}

void fmd_potsys_init(fmd_sys_t *sysp)
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

void fmd_pot_apply(fmd_sys_t *sysp, unsigned atomkind1, unsigned atomkind2, fmd_pot_t *pot)
{
    // create the pottable if doesn't exist
    if (sysp->potsys.pottable == NULL) pottable_create(sysp);

    if (pot->kind == POTKIND_EAM_ALLOY)
    {
        // find the local indices
        unsigned loc1, loc2;

        loc1 = fmd_pot_eam_find_iloc(sysp, pot->data, atomkind1);
        loc2 = fmd_pot_eam_find_iloc(sysp, pot->data, atomkind2);

        if (loc1 == -1 || loc2 == -1)
        {
            // TO-DO: error should be handled here
            assert(loc1 != -1 && loc2 != -1);
        }

        sysp->potsys.pottable[atomkind1][atomkind2].iloc = loc1;
        sysp->potsys.pottable[atomkind1][atomkind2].jloc = loc2;
        sysp->potsys.pottable[atomkind2][atomkind1].iloc = loc2;
        sysp->potsys.pottable[atomkind2][atomkind1].jloc = loc1;

        // set eam_element
        sysp->potsys.atomkinds[atomkind1].eam_element = &((eam_t *)pot->data)->elements[loc1];
        sysp->potsys.atomkinds[atomkind2].eam_element = &((eam_t *)pot->data)->elements[loc2];
    }

    //
    sysp->potsys.pottable[atomkind1][atomkind2].kind =
      sysp->potsys.pottable[atomkind2][atomkind1].kind = pot->kind;
    sysp->potsys.pottable[atomkind1][atomkind2].data =
      sysp->potsys.pottable[atomkind2][atomkind1].data = pot->data;
}

static int potkind_compare(const void *a, const void *b)
{
    if ( *( (potkind_t *)a ) == *( (potkind_t *)b ) )
        return 0;
    else
        return 1;
}

static void pot_hybridpasses_update(fmd_sys_t *sysp)
{
    // fill hybridpasses with zeros
    memset(sysp->potsys.hybridpasses, 0, sizeof(sysp->potsys.hybridpasses));

    list_t *potkinds = sysp->potsys.potkinds;
    while (potkinds != NULL)
    {
        switch (*(potkind_t *)(potkinds->data))
        {
            case POTKIND_MORSE:
                sysp->potsys.hybridpasses[0] = 1;
                break;

            case POTKIND_LJ_6_12:
                sysp->potsys.hybridpasses[0] = 1;
                break;

            case POTKIND_EAM_ALLOY:
                sysp->potsys.hybridpasses[0] = 1;
                sysp->potsys.hybridpasses[1] = 1;
                break;
        }

        potkinds = potkinds->next;
    }
}

// TO-DO?: first, clean the potkinds list
void fmd_pot_prepareForForceComp(fmd_sys_t *sysp)
{
    // TO-DO: error should be handled here
    assert(sysp->potsys.pottable != NULL);

    sysp->potsys.potkinds_num = 0;

    for (unsigned i=0; i < sysp->potsys.atomkinds_num; i++)
        for (unsigned j=0; j <= i; j++)
        {
            potpair_t *potpair = &sysp->potsys.pottable[i][j];

            // TO-DO: error should be handled here
            assert(potpair->kind != POTKIND_NONE);

            // add the potkind to potkinds list, if isn't already included there
            if (fmd_list_find_custom(sysp->potsys.potkinds, &potpair->kind, potkind_compare) == NULL)
            {
                sysp->potsys.potkinds_num++;
                sysp->potsys.potkinds = fmd_list_prepend(sysp->potsys.potkinds, &potpair->kind);
            }

            // do these two atomkinds use EAM?
            if (potpair->kind != POTKIND_EAM_ALLOY)
                sysp->potsys.atomkinds[i].eam_element =
                sysp->potsys.atomkinds[j].eam_element = NULL;
        }

    pot_hybridpasses_update(sysp);
}
