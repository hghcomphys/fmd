#   atom.py: This file is part of Free Molecular Dynamics
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

"""Atom and Atoms"""

import numpy as np
import periodictable as pt

from .error import Error
from .section import ListSection
from .unit import *

__all__ = ['Atom', 'Atoms']


class Atom:
    """A base class that contains attributes for an atom and relevant methods.
    It includes positions, group_id, velocity, etc."""

    def __init__(self, symbol='X', position=(0.0, 0.0, 0.0), group_id=0, velocity=(0.0, 0.0, 0.0)):
        # charge = 0.0, atom_id = 0, molecule_id = 0, atom_type = 0, imx = 0, imy = 0, imz = 0
        try:
            self.symbol = symbol  # atom symbol (string)
            self.position = position  # set atom position
            self.velocity = velocity  # set atom velocity
            self.group_id = group_id  # set group id for atom

            # TODO: below attributes are not implemented yet!
            # self.__charge = float(charge)  # atom charge
            # self.__atom_id = int_ge_zero(atom_id)  # atom id
            # self.__molecule_id = int_ge_zero(molecule_id)  # molecule id
            # self.__atom_type = int_ge_zero(atom_type)  # atom type
            # self.__imx = int(imx)  # image index along-x
            # self.__imy = int(imy)  # image index along-y
            # self.__imz = int(imz)  # image index along-z

        except (ValueError, IndexError, TypeError):
            raise AssertionError("Unexpected value for initializing %s!" % self.__class__.__name__)

    # Symbol -----------------------------
    @property
    def symbol(self):
        """Return symbol of atom."""
        return self.__symbol

    @symbol.setter
    def symbol(self, symbol):
        """Set a given symbol for atom."""
        self.__symbol = str(symbol)

    def set_symbol(self, symbol):
        """Set a given symbol for atom."""
        self.symbol = symbol

    # Position ---------------------------
    @property
    def x(self):
        """Return x position of atom."""
        return self.__x

    @x.setter
    def x(self, x):
        """Set x position of atom."""
        self.__x = float(x)

    @property
    def y(self):
        """Return y position of atom."""
        return self.__y

    @y.setter
    def y(self, y):
        """Set y position of atom."""
        self.__y = float(y)

    @property
    def z(self):
        """Return z position of atom."""
        return self.__z

    @z.setter
    def z(self, z):
        """Set z position of atom."""
        self.__z = float(z)

    @property
    def position(self):
        """Return vector position of the atom."""
        return self.x, self.y, self.z

    @position.setter
    def position(self, position):
        """Set vector position of the atom."""
        self.x, self.y, self.z = Error.check_tuple_float(position)

    def get_position(self):
        """Return vector position of the atom."""
        return self.position

    def set_position(self, position):
        """Set vector position of the atom."""
        self.position = position
        return self

    def translate(self, translation_vec=(0.0, 0.0, 0.0)):
        """Translate atom along a given vector."""
        translation_vec = Error.check_tuple_float(translation_vec)
        self.position = self.x+translation_vec[0], self.y+translation_vec[1], self.z+translation_vec[2]
        return self

    # Velocity -----------------------------
    @property
    def vx(self):
        """Return x component of atom's velocity."""
        return self.__vx

    @vx.setter
    def vx(self, vx):
        """Set x component of atom's velocity."""
        self.__vx = float(vx)

    @property
    def vy(self):
        """Return y component of atom's velocity."""
        return self.__vy

    @vy.setter
    def vy(self, vy):
        """Set x component of atom's velocity."""
        self.__vy = float(vy)

    @property
    def vz(self):
        """Return z component of atom's velocity."""
        return self.__vz

    @vz.setter
    def vz(self, vz):
        """Set z component of atom's velocity."""
        self.__vz = float(vz)

    @property
    def velocity(self):
        """Return velocity vector of atom."""
        return self.vx, self.vy, self.vz

    @velocity.setter
    def velocity(self, velocity):
        """Set velocity vector of atom"""
        self.vx, self.vy, self.vz = Error.check_tuple_float(velocity)

    def get_velocity(self):
        """Return velocity vector of atom."""
        return self.velocity

    def set_velocity(self, velocity=(0.0, 0.0, 0.0)):
        """Set velocity vector of atom"""
        self.velocity = velocity
        return self

    def add_velocity(self, velocity_vec=(0.0, 0.0, 0.0)):
        """A constant velocity vector added to the atom's velocity."""
        velocity_vec = Error.check_tuple_float(velocity_vec)
        self.velocity = self.vx + velocity_vec[0], self.vy + velocity_vec[1], self.vz + velocity_vec[2]
        return self

    # Group-id ----------------------------
    @property
    def group_id(self):
        """Returns the given group-id of atom."""
        return self.__group_id

    @group_id.setter
    def group_id(self, group_id):
        """Sets a given group-id for atom."""
        self.__group_id = int(group_id)

    # String representation -------------------------
    def __str__(self):
        """String conversion for Atom."""
        out = ''
        # for attribute in [self.__atom_id, self.__molecule_id, self.__atom_type, self.__q, self.__x, self.__y,
        #                   self.__z, self.__imx, self.__imy, self.__imz, self.__label]:
        for attribute in [self.symbol, self.x, self.y, self.z, self.vx, self.vy, self.vz, self.group_id]:
            out += str(attribute) + ' '
        return out

    # Mass ---------------------------------------
    @property
    def mass(self):
        """A method that returns mass value from its atomic symbol (in atomic mass unit)."""
        try:
            # TODO: This method depends on 'periodictable' module and might be limiting.
            return  eval("pt.%s.mass" % self.symbol)
        except AttributeError:
            raise AssertionError("Cannot find mass value for '%s!'" % self.symbol)

    def get_mass(self):
        """A method that returns mass value from its atomic symbol (in atomic mass unit)."""
        return self.mass

    # Kinetic energy -----------------------------
    @property
    def kinetic(self):
        """This method returns kinetic energy of the atom (in eV unit)."""
        v = self.velocity
        return 0.5*self.mass*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])*1.0365E-4  # amu*(A/ps)^2 * 1.0365E-4 --> eV

    def get_kinetic(self):
        """This method returns kinetic energy of the atom (in eV unit)."""
        return self.kinetic


class Atoms(ListSection):
    """A derived class that contains particularly list of atoms with relevant methods.
    I includes methods that relates to the group of atoms such as number of atoms, selecting group of atoms, etc."""

    def __init__(self, atoms=None):
        """Initialize atoms with either empty, atom, or list of atoms."""

        # Initial base class
        ListSection.__init__(self, 'Atoms')  # set the name for the list section

        # Initialize accordingly
        if atoms is None:
            pass
        elif isinstance(atoms, Atom) or isinstance(atoms, list):
            self.append(atoms)
        else:
            raise AssertionError("Unexpected argument for %s!" % self.__class__.__name__)

    # Items ------------------------------------
    @property
    def atoms(self):
        """Return list of atoms."""
        return self.items

    @property
    def number_of_atoms(self):
        """Return number of atoms."""
        return self.number_of_items

    def get_number_of_atoms(self):
        return self.number_of_atoms

    def __append_atom(self, atom):
        """Append an Atom to the list of atoms."""
        if not isinstance(atom, Atom):
            raise AssertionError("Expected Atom type for %s method!" % self.__append_atom.__name__)
        self.atoms.append(atom)

    def append(self, atoms):
        """Append atoms."""
        if isinstance(atoms, list):
            for atom in atoms:
                self.__append_atom(atom)
            return self
        else:
            self.__append_atom(atoms)
        return self

    def translate(self, translation_vec=(0.0, 0.0, 0.0)):
        """Translate all atoms along a given direction."""
        for atom in self.atoms:
            atom.translate(translation_vec)
        return self

    def add_velocities(self, velocity_vec=(0.0, 0.0, 0.0)):
        """Add a given vector to all atoms's velocity."""
        for atom in self.atoms:
            atom.add_velocity(velocity_vec)
        return self

    # Operators ----------------------------------
    def __add__(self, other):
        """ Define integration '+' between two Atoms instances."""
        # TODO: simply joining two list of atoms
        new_atom_section = Atoms()
        new_atom_section.append(self.atoms)
        new_atom_section.append(other.atoms)
        return new_atom_section

    # Select -------------------------------------
    def select(self, select_obj=None):
        """Return a new Atoms based on a input Select instance."""
        # TODO: handle unexpected input select_obj
        atoms_sel = Atoms()
        for atom in self.atoms:
            if select_obj.is_selected(atom):
                atoms_sel.append(atom)
        return atoms_sel

    # Group-id -----------------------------------
    def set_group_id(self, group_id):
        """A method that sets a given group id for atoms."""
        for atom in self.atoms:
            atom.group_id = group_id
        return self

    # Mass ----------------------------------------
    def get_total_mass(self):
        """This methods return total mass."""
        mass = 0.0
        for atom in self.atoms:
            mass += atom.mass
        return mass

    @property
    def total_mass(self):
        """This methods return total mass."""
        return self.get_total_mass()

    # Center of mass ------------------------------
    def get_center_of_mass(self):
        """This method returns center of mass vector."""
        com = np.array([0.0, 0.0, 0.0])
        for atom in self.atoms:
            com += np.array(atom.position) * atom.mass
        com /= self.total_mass
        return tuple(com)

    @property
    def center_of_mass(self):
        """This method returns center of mass vector."""
        return self.get_center_of_mass()

    def set_center_of_mass(self, position=(0.0, 0.0, 0.0)):
        """This method translates the center of mass to a given vector position."""
        com = self.center_of_mass
        pos = Error.check_tuple_float(position)
        for atom in self.atoms:
            atom.translate((pos[0] - com[0], pos[1] - com[1], pos[1] - com[1]))
        return self

    @center_of_mass.setter
    def center_of_mass(self, position):
        """This method translates the center of mass to a given vector position."""
        self.set_center_of_mass(position)

    def get_center_of_mass_velocity(self):
        """This method returns center of mass vector."""
        vel_com = np.array([0.0, 0.0, 0.0])
        for atom in self.atoms:
            vel_com += np.array(atom.velocity) * atom.mass
        vel_com /= self.total_mass
        return tuple(vel_com)

    @property
    def center_of_mass_velocity(self):
        """This method returns center of mass vector."""
        return self.get_center_of_mass_velocity()

    def set_center_of_mass_velocity(self, velocity=(0.0, 0.0, 0.0)):
        """This method translates the center of mass to a given vector position."""
        vel_com = self.center_of_mass_velocity
        vel = Error.check_tuple_float(velocity)
        for atom in self.atoms:
            atom.add_velocity((vel[0] - vel_com[0], vel[1] - vel_com[1], vel[2] - vel_com[2]))
        return self

    @center_of_mass_velocity.setter
    def center_of_mass_velocity(self, velocity):
        """This method translates the center of mass to a given vector position."""
        self.set_center_of_mass_velocity(velocity)

    # Temperature -------------------------------
    def get_temperature(self):
        """A method that returns the temperature (in Kelvin)."""
        total_kinetic = 0.0
        for atom in self.atoms:
            total_kinetic += atom.kinetic  # eV
        return total_kinetic / (1.5 * self.number_of_atoms) * EV_TO_KELVIN  # Kelvin

    @property
    def temperature(self):
        """A method that returns the temperature (in Kelvin)."""
        return self.get_temperature()

    def set_temperature(self, temperature, seed=None):
        """Set the temperature via assigning random velocities to the atoms."""
        assert temperature > 0.0, "Unexpected value for temperature!"
        # Set random seed
        if seed is not None:
              np.random.seed(seed)
        # Assign random velocities to the atoms (Maxwell-Boltzmann distributions)
        for atom in self.atoms:
            scale = np.sqrt(temperature*KELVIN_TO_JOULE/(atom.mass*AMU_TO_KILOGRAM)) * 0.01  # m/s * 0.01 -- > A/ps
            atom.velocity = scale * np.random.normal(0, 1, 3)  # Normal distribution
        # Set COM velocity to zero
        self.center_of_mass_velocity = (0.0, 0.0, 0.0)
        # Rescale velocity based on the given temperature
        temperature0 = self.temperature
        for atom in self.atoms:
            atom.velocity = np.array(atom.velocity) * np.sqrt(temperature/temperature0)
        return self

    @temperature.setter
    def temperature(self, temperature):
        """Set the temperature via assigning random velocities to the atoms."""
        self.set_temperature(temperature)

    # get/set operator ------------------------------
    def __getitem__(self, key):
        """A method that gets specified atom (index) from Atoms."""
        try:
            return self.atoms[int(key)]
        except:
            raise ValueError("Atom index not found!")

    def __setitem__(self, key, item):
        """This method set atom (using given index and Atom object) in Atoms."""
        if not isinstance(item, Atom):
            raise AssertionError("Unexpected Atom type!")
        try:
            self.atoms[int(key)] = item
        except:
            raise AssertionError("Atom index not found!")


