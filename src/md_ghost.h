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

void fmd_ghostparticles_init(TSystem *sysp);
void fmd_ghostparticles_update_Femb(TSystem *sysp);
void fmd_ghostparticles_update_LocOrdParam(TSystem *sysp);
void fmd_ghostparticles_delete(TSystem *sysp);
void fmd_particles_migrate(TSystem *sysp);
