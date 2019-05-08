/*
  potential.h: This file is part of Free Molecular Dynamics

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

#ifndef POTENTIAL_H
#define POTENTIAL_H

typedef char fmd_particle_name_t[17];

typedef struct
{
    double mass;
    fmd_particle_name_t name;
} particle_form_t;          // particle form

typedef struct
{
    unsigned pforms_num;
    particle_form_t *pforms;
} potential_t;
#endif /* POTENTIAL_H */
