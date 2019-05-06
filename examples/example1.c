/* Assuming that FMD is already compiled in src directory, this example
   can be compiled with the following command:

   $ gcc example1.c -L../src/  -Wl,-R../src/ -lfmd -lm -O3 -o example1.x

   and can be executed by

   $ mpirun -n 2 ./example1.x
*/

#include <math.h>
#include "../src/fmd.h"

int main(int argc, char *argv[])
{
    fmdt_sys *sys;

    // create an fmd-system instance
    sys = fmd_sys_create();

    // set size of the simulation box (in Angstrom)
    double latticeParameter = 3.6316;
    fmd_box_setSize(sys, 10*latticeParameter, 10*latticeParameter, 10*latticeParameter);

    // set periodic boundary conditions in three dimensions (0 = no PBC)
    fmd_box_setPBC(sys, 1, 1, 1);

    // partition the simulation box into subdomains for MPI-based parallel computation
    fmd_box_setSubDomains(sys, 1, 2, 1);

    /* sometimes the user launches more processes than the chosen number of subdomains; they're not needed here!
       the function fmd_proc_isMD() can be called only after fmd_box_setSubDomains() */
    if (! fmd_proc_isMD(sys))
    {
        fmd_sys_free(sys, 1);
        return 0;
    }

    // load the EAM file into memory; can be called only after fmd_box_setSubDomains()
    fmd_pot_eam_init(sys, "../potentials/Cu01.eam.alloy");

    // get the EAM potential cutoff radius and create the box grid by using it
    double cutoff = fmd_pot_eam_getCutoffRadius(sys);
    fmd_box_createGrid(sys, cutoff);

    // set the desired temperature (in Kelvin)
    fmd_matt_setDesiredTemperature(sys, 350.0);

    // make an fcc cuboid at a given position and with a given size
    fmd_matt_makeCuboidFCC(sys, 0.0, 0.0, 0.0, 10, 10, 10, latticeParameter, 0, 0);

    // distribute the matter among subdomains
    fmd_matt_distribute(sys);

    // set time step to 2 femtoseconds
    fmd_dync_setTimeStep(sys, 2e-3);

    // let us simulate for 2 picoseconds
    double final_time = 2.0;

    // set where to save output files (default = current directory)
    //fmd_io_setSaveDirectory(sys, "output/");

    // let configurations be saved as XYZ files
    fmd_io_setSaveConfigMode(sys, scmXYZParticlesNum);

    // set Berendsen thermostat parameter
    fmd_dync_setBerendsenThermostatParameter(sys, 2e-2);

    // compute forces for the first time
    fmd_dync_updateForces(sys);

    // the time loop starts here
    // fmd_dync_getTime() returns current internal time of the fmd-system instance
    while (fmd_dync_getTime(sys) < final_time)
    {
        // save configuration every 40 femtoseconds
        if (fmod(fmd_dync_getTime(sys), 0.04) < fmd_dync_getTimeStep(sys))
            fmd_matt_saveConfiguration(sys);

        // report some quantities every time step
        fmd_io_printf(sys, "%f\t%f\t%e\n", fmd_dync_getTime(sys),
                                           fmd_matt_getGlobalTemperature(sys),
                                           fmd_matt_getTotalEnergy(sys));

        // take first step of velocity Verlet integrator
        fmd_dync_velocityVerlet_takeFirstStep(sys, 1);

        // compute forces
        fmd_dync_updateForces(sys);

        // take last step of velocity Verlet integrator
        fmd_dync_velocityVerlet_takeLastStep(sys);

        // increase internal time by one time step
        fmd_dync_incTime(sys);
    }
    // end of the time loop

    // save system's final state in a file
    fmd_io_saveState(sys, "state0.stt");

    // release memory taken for EAM potential
    fmd_pot_eam_free(sys);

    // another report
    fmd_io_printf(sys, "The run took about %.3f seconds to finish.\n", fmd_proc_getWallTime(sys));

    // release memory taken for the fmd-system instance (including subdomain and all particles)
    // also finalize MPI
    fmd_sys_free(sys, 1);

    return 0;
}