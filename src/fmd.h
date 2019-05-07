/*
  fmd.h: This file is part of Free Molecular Dynamics

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

// types

typedef struct fmdt_sys fmdt_sys;
typedef enum
    {scmXYZParticlesNum, scmXYZSeparate, scmCSV, scmVTF} fmdt_SaveConfigMode;

// functions

void fmd_matt_addVelocity(fmdt_sys *system, int groupID, double vx, double vy, double vz);
void fmd_matt_setActiveGroup(fmdt_sys *system, int groupID);
void fmd_matt_setDesiredTemperature(fmdt_sys *system, double DesiredTemperature);
void fmd_matt_makeCuboidFCC(fmdt_sys *system, double x, double y, double z,
  int dimx, int dimy, int dimz, double latticeParameter, int elementID, int groupID);
void fmd_matt_saveConfiguration(fmdt_sys *system);
double fmd_matt_getTotalEnergy(fmdt_sys *system);
double fmd_matt_getGlobalTemperature(fmdt_sys *system);
void fmd_matt_distribute(fmdt_sys *system);
void fmd_matt_giveTemperature(fmdt_sys *system, int groupID);

void fmd_box_setPBC(fmdt_sys *system, int PBCx, int PBCy, int PBCz);
void fmd_box_setSize(fmdt_sys *system, double sx, double sy, double sz);
void fmd_box_setSubDomains(fmdt_sys *system, int dimx, int dimy, int dimz);
void fmd_box_createGrid(fmdt_sys *system, double cutoff);

void fmd_io_setSaveDirectory(fmdt_sys *system, char *directory);
void fmd_io_setSaveConfigMode(fmdt_sys *system, fmdt_SaveConfigMode mode);
void fmd_io_printf(fmdt_sys *system, const char * restrict format, ...);
void fmd_io_loadState(fmdt_sys *system, char *file, int useTime);
void fmd_io_saveState(fmdt_sys *system, char *filename);

void fmd_pot_eam_init(fmdt_sys *system, char *filePath);
double fmd_pot_eam_getLatticeParameter(fmdt_sys *system, int element);
double fmd_pot_eam_getCutoffRadius(fmdt_sys *system);
void fmd_pot_eam_free(fmdt_sys *system);
void fmd_pot_setCutoffRadius(fmdt_sys *system, double cutoff);

void fmd_subd_init(fmdt_sys *system);
void fmd_subd_free(fmdt_sys *system);

double fmd_proc_getWallTime(fmdt_sys *system);
int fmd_proc_isMD(fmdt_sys *system);

fmdt_sys *fmd_sys_create();
void fmd_sys_free(fmdt_sys *system, int finalizeMPI);

double fmd_dync_getTimeStep(fmdt_sys *system);
void fmd_dync_setTimeStep(fmdt_sys *system, double timeStep);
double fmd_dync_getTime(fmdt_sys *system);
void fmd_dync_updateForces(fmdt_sys *system);
void fmd_dync_incTime(fmdt_sys *system);
void fmd_dync_setBerendsenThermostatParameter(fmdt_sys *system, double parameter);
void fmd_dync_velocityVerlet_takeFirstStep(fmdt_sys *system, int useThermostat);
int fmd_dync_velocityVerlet_takeLastStep(fmdt_sys *system);
void fmd_dync_equilibrate(fmdt_sys *system, int groupID, double duration, double strength);
