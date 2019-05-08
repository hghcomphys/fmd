/*
  md_ghost.h: This file is part of Free Molecular Dynamics

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

#ifndef MD_GHOST_H
#define MD_GHOST_H

void fmd_ghostparticles_init(fmd_sys_t *sysp);
void fmd_ghostparticles_update_Femb(fmd_sys_t *sysp);
void fmd_ghostparticles_update_LocOrdParam(fmd_sys_t *sysp);
void fmd_ghostparticles_delete(fmd_sys_t *sysp);
void fmd_particles_migrate(fmd_sys_t *sysp);

#endif /* MD_GHOST_H */
