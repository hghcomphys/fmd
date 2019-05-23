/*
  fmd.h: This file is part of Free Molecular Dynamics

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

#ifndef FMD_H
#define FMD_H

// symbolic constants

#define FMD_PHYS_AMU            1.036426957207970e-04       // atomic mass unit

// types

typedef struct fmd_sys_t fmd_sys_t;
typedef struct fmd_pot_t fmd_pot_t;
typedef enum
    {scmXYZParticlesNum, scmXYZSeparate, scmCSV, scmVTF} fmd_SaveConfigMode_t;
typedef char fmd_atomkind_name_t[16];

// functions

void fmd_matt_addVelocity(fmd_sys_t *system, int groupID, double vx, double vy, double vz);
void fmd_matt_setActiveGroup(fmd_sys_t *system, int groupID);
void fmd_matt_setDesiredTemperature(fmd_sys_t *system, double DesiredTemperature);
void fmd_matt_makeCuboidFCC(fmd_sys_t *system, double x, double y, double z,
  int dimx, int dimy, int dimz, double latticeParameter, int elementID, int groupID);
void fmd_matt_makeCuboidFCC_alloy(fmd_sys_t *system, double x, double y, double z,
  int dimx, int dimy, int dimz, double latticeParameter, double *proportions, int groupID);
void fmd_matt_saveConfiguration(fmd_sys_t *system);
double fmd_matt_getTotalEnergy(fmd_sys_t *system);
double fmd_matt_getGlobalTemperature(fmd_sys_t *system);
void fmd_matt_distribute(fmd_sys_t *system);
void fmd_matt_giveTemperature(fmd_sys_t *system, int groupID);

void fmd_box_setPBC(fmd_sys_t *system, int PBCx, int PBCy, int PBCz);
void fmd_box_setSize(fmd_sys_t *system, double sx, double sy, double sz);
void fmd_box_setSubDomains(fmd_sys_t *system, int dimx, int dimy, int dimz);
void fmd_box_createGrid(fmd_sys_t *system, double cutoff);

void fmd_io_setSaveDirectory(fmd_sys_t *system, char *directory);
void fmd_io_setSaveConfigMode(fmd_sys_t *system, fmd_SaveConfigMode_t mode);
void fmd_io_printf(fmd_sys_t *system, const char * restrict format, ...);
void fmd_io_loadState(fmd_sys_t *system, char *file, int useTime);
void fmd_io_saveState(fmd_sys_t *system, char *filename);

fmd_pot_t *fmd_pot_eam_alloy_load(fmd_sys_t *system, char *filePath);
double fmd_pot_eam_getLatticeParameter(fmd_sys_t *system, int element);
double fmd_pot_eam_getCutoffRadius(fmd_sys_t *system, fmd_pot_t *pot);
void fmd_pot_eam_free(fmd_sys_t *system);
void fmd_pot_setCutoffRadius(fmd_sys_t *system, double cutoff);
void fmd_pot_setAtomKinds(fmd_sys_t *system, unsigned number, fmd_atomkind_name_t *names, double *masses);
fmd_pot_t *fmd_pot_lj_apply(fmd_sys_t *system, unsigned atomkind1, unsigned atomkind2,
  double sigma, double epsilon, double cutoff);
fmd_pot_t *fmd_pot_morse_apply(fmd_sys_t *sysp, unsigned atomkind1, unsigned atomkind2,
                              double D0, double alpha, double r0, double cutoff);
void fmd_pot_apply(fmd_sys_t *system, unsigned atomkind1, unsigned atomkind2, fmd_pot_t *pot);


void fmd_subd_init(fmd_sys_t *system);
void fmd_subd_free(fmd_sys_t *system);

double fmd_proc_getWallTime(fmd_sys_t *system);
int fmd_proc_isMD(fmd_sys_t *system);

fmd_sys_t *fmd_sys_create();
void fmd_sys_free(fmd_sys_t *system, int finalizeMPI);

double fmd_dync_getTimeStep(fmd_sys_t *system);
void fmd_dync_setTimeStep(fmd_sys_t *system, double timeStep);
double fmd_dync_getTime(fmd_sys_t *system);
void fmd_dync_updateForces(fmd_sys_t *system);
void fmd_dync_updateForcesLJ(fmd_sys_t *system);
void fmd_dync_incTime(fmd_sys_t *system);
void fmd_dync_setBerendsenThermostatParameter(fmd_sys_t *system, double parameter);
void fmd_dync_velocityVerlet_takeFirstStep(fmd_sys_t *system, int useThermostat);
int fmd_dync_velocityVerlet_takeLastStep(fmd_sys_t *system);
void fmd_dync_equilibrate(fmd_sys_t *system, int groupID, double duration, double strength);

#endif /* FMD_H */
