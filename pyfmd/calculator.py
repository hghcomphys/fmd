#   calculator.py: This file is part of Free Molecular Dynamics
#
#   Copyright (C) 2019 Hossein Ghorbanfekr
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""FMD core library wrapper"""

import ctypes as ct

from .error import Error

__all__ = ["Calculator"]


class Calculator:
    """This base class provides an interface to FMD core library.
    Basically a direct molecular dynamics simulation can be performed via this class."""

    def __init__(self):
        """Initialize Calculator class by loading library file and creating an instance of system."""
        try:
            lib_file = 'libfmd.so'
            self._lib_file = str(lib_file)
            self._lib = ct.CDLL(self._lib_file)  # load library file
            self._sys = self._create_system()
            self._is_potential = False  # check a potential loaded is loaded
            self._dummy = 0  # dummy variable
        except:
            raise

    def __del__(self):
        """Release memory taken by fmd-system instance."""
        self.free_potential()  # release memory taken by potential
        self.free_system()  # release memory taken by system

    # FMD-system instance pointer --------------------
    @property
    def _sys_none(self):
        """Assert error if pointer to fmd-system is None."""
        if self._sys is None:
            raise AssertionError("No fmd-system is found!")
        return self._sys

    # FMD-system instance ----------------------------
    def _create_system(self):
        """Create an fmd-system instance."""
        self._lib.fmd_sys_create.argtype = ()
        self._lib.fmd_sys_create.restype = ct.c_void_p
        return self._lib.fmd_sys_create()

    def free_system(self, finalize_mpi=True):
        """Release memory taken for the fmd-system instance (including subdomain and all particles)."""
        if self._sys is not None:
            self._lib.fmd_sys_free.argtypes = (ct.c_void_p, ct.c_int)
            self._lib.fmd_sys_free(self._sys, int(finalize_mpi))
            self._sys = None  # set None flag for fmd-system pointer

    # Potential -----------------------------------------
    def init_potential(self, file_name):
        """Load the EAM file into memory."""
        self._lib.fmd_pot_eam_init.argtypes = (ct.c_void_p, ct.c_char_p)
        self._lib.fmd_pot_eam_init(self._sys_none, str(file_name).encode('utf-8'))
        self._is_potential = True  # check a potential is loaded
        return self

    def free_potential(self):
        """Release memory taken for potential."""
        if self._is_potential:
            self._lib.fmd_pot_eam_free.argtypes = (ct.c_void_p, ct.c_double)
            self._lib.fmd_pot_eam_free(self._sys_none, self._dummy)
            self._is_potential = False
        return self

    @property
    def potential_cutoff(self):
        """Get the potential cutoff radius."""
        return self.get_potential_cutoff()

    def get_potential_cutoff(self):
        """Get the potential cutoff radius."""
        self._lib.fmd_pot_eam_getCutoffRadius.argtypes = (ct.c_void_p, ct.c_double)
        self._lib.fmd_pot_eam_getCutoffRadius.restype = ct.c_double
        return self._lib.fmd_pot_eam_getCutoffRadius(self._sys_none, self._dummy)

    # NEW
    def set_potential_atom_kinds(self, number, names, masses):
        self._lib.fmd_pot_setAtomKinds.argtypes = (ct.c_void_p, ct.c_int, ct.c_void_p, ct.c_void_p)
        self._lib.fmd_pot_setAtomKinds(self._sys_none, number, [str(name).encode('utf-8') for name in names], masses)
        return self

    # NEW
    def apply_potential_lj(self, atom_kind1, atom_kind2, sigma, epsilon, cutoff):
        self._lib.fmd_pot_lj_apply.argtypes = (ct.c_void_p, ct.c_int, ct.c_int, ct.c_double, ct.c_double, ct.c_double)
        self._lib.fmd_pot_lj_apply.restype = ct.c_void_p
        return self._lib.fmd_pot_lj_apply(self._sys_none, atom_kind1, atom_kind2, sigma, epsilon, cutoff)

    # Box ----------------------------------------
    @property
    def box_size(self):
        # TODO: not implemented yet!
        return None

    @box_size.setter
    def box_size(self, box):
        """Set size of the simulation box (in Angstrom)."""
        self.set_box_size(box)

    def set_box_size(self, box):
        """Set size of the simulation box (in Angstrom)."""
        box = Error.check_tuple(box, float)
        self._lib.fmd_box_setSize.argtypes = (ct.c_void_p, ct.c_double, ct.c_double, ct.c_double)
        self._lib.fmd_box_setSize(self._sys_none, box[0], box[1], box[2])
        return self

    @property
    def box_pbc(self):
        # TODO: not implemented yet!
        return None

    @box_pbc.setter
    def box_pbc(self, pbc):
        """Set periodic boundary conditions in three dimensions."""
        self.set_box_pbc(pbc)

    def set_box_pbc(self, pbc):
        """Set periodic boundary conditions in three dimensions."""
        pbc = Error.check_tuple(pbc, int)
        self._lib.fmd_box_setPBC.argtypes = (ct.c_void_p, ct.c_int, ct.c_int, ct.c_int)
        self._lib.fmd_box_setPBC(self._sys_none, pbc[0], pbc[1], pbc[2])
        return self

    def set_box_grid(self, cutoff):
        """Create box grid."""
        self._lib.fmd_box_createGrid.argtypes = (ct.c_void_p, ct.c_double)
        self._lib.fmd_box_createGrid(self._sys_none, float(cutoff))
        return self

    # Sub-domains -------------------------------------
    def init_subdomain(self):
        """Prepare subdomain before adding some matter."""
        self._lib.fmd_subd_init.argtypes = (ct.c_void_p, ct.c_double)
        self._lib.fmd_subd_init(self._sys_none, self._dummy)

    @property
    def subdomains(self):
        # TODO: not implemented yet!
        return None

    @subdomains.setter
    def subdomains(self, dim):
        """Partition the simulation box into subdomains for MPI-based parallel computation."""
        self.set_subdomains(dim)

    def set_subdomains(self, dim):
        """Partition the simulation box into subdomains for MPI-based parallel computation."""
        dim = Error.check_tuple(dim, Error.int_gt_zero)
        self._lib.fmd_box_setSubDomains.argtypes = (ct.c_void_p, ct.c_int, ct.c_int, ct.c_int)
        self._lib.fmd_box_setSubDomains(self._sys_none, dim[0], dim[1], dim[2])
        return self

    # Process -----------------------------------------
    @property
    def is_process_md(self):
        """Sometimes the user launches more processes than the chosen number of subdomains; they're not needed here!"""
        return self.get_is_process_md()

    def get_is_process_md(self):
        """Sometimes the user launches more processes than the chosen number of subdomains; they're not needed here!"""
        self._lib.fmd_proc_isMD.argtypes = (ct.c_void_p, ct.c_double)
        self._lib.fmd_proc_isMD.restype = ct.c_int
        return bool(self._lib.fmd_proc_isMD(self._sys_none, self._dummy))

    @property
    def process_wall_time(self):
        """Get the wall time for each process."""
        return self.get_process_wall_time()

    def get_process_wall_time(self):
        """Get the wall time for each process."""
        # TODO: io_fprintf
        # void fmd_io_printf.argtypes= (fmdt_sys *system, const char * restrict format, ...);
        self._lib.fmd_proc_getWallTime.argtypes = (ct.c_void_p, ct.c_double)
        self._lib.fmd_proc_getWallTime.restype = ct.c_double
        return self._lib.fmd_proc_getWallTime(self._sys_none, self._dummy)

    # Temperature ---------------------------------------
    @property
    def desired_temperature(self):
        # TODO: not implemented yet!
        return None

    @desired_temperature.setter
    def desired_temperature(self, temperature):
        """Set the desired temperature (in Kelvin)."""
        self.set_desired_temperature(temperature)

    def set_desired_temperature(self, temperature):
        """Set the desired temperature (in Kelvin)."""
        self._lib.fmd_matt_setDesiredTemperature.argtypes = (ct.c_void_p, ct.c_double)
        self._lib.fmd_matt_setDesiredTemperature(self._sys_none, Error.float_gt_zero(temperature))
        return self

    @property
    def temperature(self):
        """Return global Temperature."""
        return self.get_temperature()

    def get_temperature(self):
        """Return global Temperature."""
        self._lib.fmd_matt_getGlobalTemperature.argtypes = (ct.c_void_p, ct.c_double)
        self._lib.fmd_matt_getGlobalTemperature.restype = ct.c_double
        return self._lib.fmd_matt_getGlobalTemperature(self._sys_none, self._dummy)

    # Material builder (demo) ----------------------------
    def make_cuboid_fcc(self, pos, dim, lattice_parameter, element_id, group_id):
        """Make an fcc cuboid at a given position and with a given size."""
        dim = Error.check_tuple(dim, Error.int_gt_zero)
        pos = Error.check_tuple(pos, float)
        self._lib.fmd_matt_makeCuboidFCC.argtypes = (ct.c_void_p, ct.c_double, ct.c_double, ct.c_double,
                                               ct.c_int, ct.c_int, ct.c_int, ct.c_double, ct.c_int, ct.c_int)
        # TODO: input arguments error handel
        self._lib.fmd_matt_makeCuboidFCC(self._sys_none, pos[0], pos[1], pos[2], dim[0], dim[1], dim[2],
                                         Error.float_gt_zero(lattice_parameter), element_id, group_id)

    def material_distribute(self):
        """Distribute the matter among subdomains."""
        self._lib.fmd_matt_distribute.argtypes = (ct.c_void_p, ct.c_double)
        self._lib.fmd_matt_distribute(self._sys_none, self._dummy)
        return self

    # Time (step) ------------------------------------------
    @property
    def time_step(self):
        """Get time step for dynamics."""
        return self.get_time()

    @time_step.setter
    def time_step(self, time_step):
        """Set time step (in picoseconds)."""
        self.set_time_step(time_step)

    def get_time_step(self):
        """Get time step for dynamics."""
        self._lib.fmd_dync_getTimeStep.argtypes = (ct.c_void_p, ct.c_double)
        self._lib.fmd_dync_getTimeStep.restype = ct.c_double
        return self._lib.fmd_dync_getTimeStep(self._sys_none, self._dummy)

    def set_time_step(self, time_step=0.001):
        """Set time step (in picoseconds)."""
        self._lib.fmd_dync_setTimeStep.argtypes = (ct.c_void_p, ct.c_double)
        self._lib.fmd_dync_setTimeStep(self._sys_none, Error.float_gt_zero(time_step))
        return self

    @property
    def time(self):
        return self.get_time()

    def get_time(self):
        """Get time for dynamics."""
        self._lib.fmd_dync_getTime.argtypes = (ct.c_void_p, ct.c_double)
        self._lib.fmd_dync_getTime.restype = ct.c_double
        return self._lib.fmd_dync_getTime(self._sys_none, self._dummy)

    def next_time_step(self):
        """Time increment for time integration."""
        self._lib.fmd_dync_incTime.argtypes = (ct.c_void_p, ct.c_double)
        self._lib.fmd_dync_incTime(self._sys_none, self._dummy)
        return self

    # Thermostat ------------------------------------------
    def set_berendsen_thermostat_parameter(self, parameter=0.01):
        """Set Berendsen thermostat parameter (in picoseconds)."""
        self._lib.fmd_dync_setBerendsenThermostatParameter.argtypes = (ct.c_void_p, ct.c_double)
        self._lib.fmd_dync_setBerendsenThermostatParameter(self._sys_none, Error.float_gt_zero(parameter))
        return self

    # I/O -------------------------------------------------
    @property
    def io_directory(self):
        # TODO: not implemented yet!
        return None

    @io_directory.setter
    def io_directory(self, path):
        """Set where to save output files (default = current directory)."""
        self.set_io_directory(path)

    def set_io_directory(self, path):
        """Set where to save output files (default = current directory)."""
        self._lib.fmd_io_setSaveDirectory.argtypes = (ct.c_void_p, ct.c_char_p)
        self._lib.fmd_io_setSaveDirectory(self._sys_none, str(path).encode('utf-8'))
        return self

    def set_io_config_mode(self, mode=0):
        """Set saving configuration mode."""
        self._lib.fmd_io_setSaveConfigMode.argtypes = (ct.c_void_p, ct.c_int)
        self._lib.fmd_io_setSaveConfigMode(self._sys_none, mode)  # TODO: enum fmdt_SaveConfigMode
        return self

    def save_io_config(self):
        """Save configuration file."""
        self._lib.fmd_matt_saveConfiguration.argtypes = (ct.c_void_p, ct.c_double)
        self._lib.fmd_matt_saveConfiguration(self._sys_none, self._dummy)
        return self

    def save_io_state(self, file_name):
        """Save state in to a file (in .stt format)."""
        self._lib.fmd_io_saveState.argtypes = (ct.c_void_p, ct.c_char_p)
        self._lib.fmd_io_saveState(self._sys_none, str(file_name).encode('utf-8'))
        return self

    def load_io_state(self, file_name, use_time=True):
        """Load state file (in .stt format)."""
        self._lib.fmd_io_loadState.argtypes = (ct.c_void_p, ct.c_char_p, ct.c_int)
        self._lib.fmd_io_loadState(self._sys_none, str(file_name).encode('utf-8'), int(use_time))
        return self

    # Dynamics -----------------------------------------------
    def update_force(self):
        """Compute forces."""
        self._lib.fmd_dync_updateForces.argtypes = (ct.c_void_p, ct.c_double)
        self._lib.fmd_dync_updateForces(self._sys_none, self._dummy)
        return self

    def velocity_verlet_first_step(self, use_thermostat=False):
        """Velocity-Verlet (first step)."""
        self._lib.fmd_dync_velocityVerlet_takeFirstStep.argtypes = (ct.c_void_p, ct.c_int)
        self._lib.fmd_dync_velocityVerlet_takeFirstStep(self._sys_none, int(use_thermostat))
        return self

    def velocity_verlet_second_step(self):
        """Velocity-Verlet (last step)."""
        self._lib.fmd_dync_velocityVerlet_takeLastStep.argtypes = (ct.c_void_p, ct.c_double)
        self._lib.fmd_dync_velocityVerlet_takeLastStep.restype = ct.c_int
        self._lib.fmd_dync_velocityVerlet_takeLastStep(self._sys_none, self._dummy)
        return self

    def velocity_verlet(self, use_thermostat=False):
        """Apply a complete velocity-Verlet algorithm."""
        self.velocity_verlet_first_step(use_thermostat)
        self.update_force()
        self.velocity_verlet_second_step()
        return self

    def equilibrate(self, group_id, duration, strength):
        """Equilibrate a group of atoms for specified time period."""
        self._lib.fmd_dync_equilibrate.argtypes = (ct.c_void_p, ct.c_int, ct.c_double, ct.c_double)
        self._lib.fmd_dync_equilibrate(self._sys_none, group_id, Error.float_gt_zero(duration), strength)
        return self

    # Velocity --------------------------------------------------
    def add_velocity(self, group_id, velocity_vec=(0.0, 0.0, 0.0)):
        """Add some center-of-mass velocity to a specified group of atoms."""
        vel = Error.check_tuple(velocity_vec, float)
        self._lib.fmd_matt_addVelocity.argtypes = (ct.c_void_p, ct.c_int, ct.c_double, ct.c_double, ct.c_double)
        self._lib.fmd_matt_addVelocity(self._sys_none, group_id, vel[0], vel[1], vel[2])
        return self

    def set_velocity(self, group_id):
        """Set velocity of atoms based on desired temperature."""
        self._lib.fmd_matt_giveTemperature.argtypes = (ct.c_void_p, ct.c_int)
        self._lib.fmd_matt_giveTemperature(self._sys_none, group_id)
        return self

    # Group-ID --------------------------------------------------
    def set_active_group(self, group_id):
        """Apply a complete velocity-Verlet algorithm."""
        self._lib.fmd_matt_setActiveGroup.argtypes = (ct.c_void_p, ct.c_int)
        self._lib.fmd_matt_setActiveGroup(self._sys_none, group_id)
        return self

    # Energy ----------------------------------------------------
    @property
    def total_energy(self):
        """Return the total Energy."""
        return self.get_total_energy()

    def get_total_energy(self):
        """Return the total Energy."""
        self._lib.fmd_matt_getTotalEnergy.argtypes = (ct.c_void_p, ct.c_double)
        self._lib.fmd_matt_getTotalEnergy.restype = ct.c_double
        return self._lib.fmd_matt_getTotalEnergy(self._sys_none, self._dummy)

