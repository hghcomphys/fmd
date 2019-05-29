"""PyFMD example

This example simulates collision between two identical fcc-copper crystals in NVE ensemble.

How to execute:

    $ mpirun -n 2 python 02_collision.py
"""

import sys
from pyfmd import Calculator

# create an fmd-system instance
md = Calculator()

# set size of the simulation box (in Angstrom)
latticeParameter = 3.6316
boxx = 40 * latticeParameter
boxy = 20 * latticeParameter
boxz = 20 * latticeParameter
md.box_size = boxx, boxy, boxz

# set periodic boundary conditions in three dimensions (False = no PBC)
md.box_pbc = False, False, False

# partition the simulation box into subdomains for MPI-based parallel computation
md.subdomains = (1, 1, 2)

# sometimes the user launches more processes than the chosen number of subdomains; they're not needed here!
# the function fmd_proc_isMD() can be called only after setting subdomains.
if not md.is_process_md:
    md.free_system(finalize_mpi=True)
    sys.exit()

# add copper atoms
md.set_potential_atom_kinds({'Cu': 63.546})  # symbols and masses (a.m.u)

# load the EAM file into memory; can be called only after fmd_box_setSubDomains()
pot = md.load_potential_eam_alloy("../potentials/Cu01.eam.alloy")

# apply the potential
md.apply_potential(0, 0, pot)

# get the EAM potential cutoff radius and define the box grid by using it
cutoff = md.get_potential_cutoff(pot)
md.create_box_grid(cutoff)

# set the desired temperature (in Kelvin)
# md.desired_temperature = 300.0  # default

# make an fcc cuboid at a given position and with a given size
cusize = 7
posx = 10.0
posy = (boxy - cusize*latticeParameter) / 2.0
posz = (boxz - cusize*latticeParameter) / 2.0
md.make_cuboid_fcc((posx, posy+cusize*latticeParameter*.35, posz), (cusize, cusize, cusize), latticeParameter, 0, 0)

# add another fcc cuboid with a different groupID
posx = boxx - cusize*latticeParameter - posx
md.make_cuboid_fcc((posx, posy-cusize*latticeParameter*.35, posz), (cusize, cusize, cusize), latticeParameter, 0, 1)

# distribute the matter among subdomains
md.material_distribute()

# set time step to 2 femtoseconds
md.time_step = 0.002

# set where to save output files (default = current directory)
md.io_directory = "../"

# let configurations be saved as VTF files
md.set_io_config_mode(0)

# set Berendsen thermostat parameter
md.set_berendsen_thermostat_parameter(0.02)

# equilibrate the two clusters
print("equilibrating the first cluster...")
md.equilibrate(0, 3.0, 2e-2)
print("equilibrating the second cluster...")
md.equilibrate(1, 3.0, 2e-2)

# add some center-of-mass velocity to atoms of groups 0 and 1
md.add_velocity(0, (+15., 0., 0.))
md.add_velocity(1, (-15., 0., 0.))

# activate all groups for dynamics; -1 as a groupID means all groups
md.set_active_group(-1)

# compute forces for the first time
md.update_force()

final_step = 8.0
while md.time < final_step:

    # save configuration every 40 femtoseconds
    if md.get_time() % 0.06 < md.get_time_step():
        md.save_io_config()

    # report some quantities every time step
    # print (md.time, md.temperature, md.total_energy)

    # time integrator
    md.velocity_verlet(use_thermostat=False)  # NVE

    # increase internal time by one time step
    md.next_time_step()

# save system's final state in a file
md.save_io_state("state0.stt")

# another report
print("The run took about %.3f seconds to finish." % md.process_wall_time)

# release memory taken for potential and fmd-system
# del md
