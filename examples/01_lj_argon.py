"""PyFMD

This example demonstrates an MD simulation of liquid Argon using Lennard-Jones potential.
Beforehand, make sure that PyFMD is installed and properly linked to the FMD library file.

How to run:

    $ mpirun -n 2 python 01_lj_argon.py
"""

import sys
from pyfmd import Calculator

# create an fmd-system instance
md = Calculator()

# set size of the simulation box (in Angstrom) and periodic boundary conditions
lattice_parameter = 5.26
md.box_size = 10*lattice_parameter, 10*lattice_parameter, 10*lattice_parameter

# set periodic boundary conditions in three dimensions (False = no PBC)
md.box_pbc = True, True, True

# partition the simulation box into subdomains for MPI-based parallel computation
md.subdomains = (1, 2, 1)

# sometimes the user launches more processes than the chosen number of subdomains; they're not needed here!
# the function fmd_proc_isMD() can be called only after setting subdomains.
if not md.is_process_md:
    md.free_system(finalize_mpi=True)
    sys.exit()

# have only argon atoms
name = ["Ar"]
mass = [39.948]
md.set_potential_atom_kinds(name, mass)

# use a 12-6 Lennard-Jones potential for Argon atoms
sigma = 3.4
epsilon = 0.0104
cutoff = 2.5 * sigma
md.apply_potential_lj(0, 0, sigma, epsilon, cutoff)

# create the box grid
md.create_box_grid(cutoff)

# set the desired temperature (in Kelvin)
md.desired_temperature = 100.0

# make an fcc cuboid at a given position and with a given size
md.make_cuboid_fcc((0, 0, 0), (10, 10, 10), lattice_parameter, 0, 0)

# distribute the matter among subdomains
md.material_distribute()

# set time step to 2 femtoseconds
md.time_step = 0.002

# set where to save output files (default = current directory)
# md.io_directory = "../"

# let configurations be saved as XYZ files
md.set_io_config_mode(0)

# set Berendsen thermostat parameter
md.set_berendsen_thermostat_parameter(0.02)

# compute forces for the first time
md.update_force()

final_step = 2.0
while md.time < final_step:

    # save configuration every 40 femtoseconds
    if md.get_time() % 0.04 < md.get_time_step():
        md.save_io_config()

    # report some quantities every time step
    print (md.time, md.temperature, md.total_energy)

    # take first step of velocity Verlet integrator
    md.velocity_verlet_first_step(use_thermostat=True)  # NVT ensemble

    # compute forces
    md.update_force()

    # take last step of velocity Verlet integrator
    md.velocity_verlet_second_step()

    # Or, simply call velocity_verlet method
    # md.update_position_and_velocity()

    # increase internal time by one time step
    md.next_time_step()

# save system's final state in a file
md.save_io_state("state0.stt")

# another report
print ("The run took about %.3f seconds to finish." % md.process_wall_time)

# release memory taken for potential and fmd-system
del md
