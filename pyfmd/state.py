#   state.py: This file is part of Free Molecular Dynamics
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

"""State header section"""

from .section import FixedSection
from .error import Error
from .box import Box

__all__ = ["StateHeader"]


class StateHeader(FixedSection):
    """A class that contains header information for .stt file format."""

    def __init__(self, simulation_time=0.0, number_of_atoms=0, box=(0.0, 0.0, 0.0), pbc=(True, True, True)):
        FixedSection.__init__(self, 'StateHeader')  # set name of the section
        try:
            self.simulation_time = simulation_time
            self.number_of_atoms = number_of_atoms
            self.box_obj = box  # create a Box instance from given box sizes
            self.pbc = pbc  # periodic boundary conditions along x,y,z directions

        except ValueError:
            raise AssertionError("Unexpected value for %s!" % self.__class__.__name__)

    # State -------------------------------------
    @property
    def simulation_time(self):
        """Set simulation time in pico seconds."""
        return self.__simulation_time

    @simulation_time.setter
    def simulation_time(self, simulation_time):
        """Set simulation time in pico seconds."""
        self.__simulation_time = Error.float_ge_zero(simulation_time)

    @property
    def number_of_atoms(self):
        """Return the number of atoms."""
        return self.__number_of_atoms

    @number_of_atoms.setter
    def number_of_atoms(self, number_of_atoms):
        """Set the number of atoms."""
        self.__number_of_atoms = Error.int_ge_zero(number_of_atoms)

    @property
    def pbc(self):
        """Return periodic boundary condition for the simulation box."""
        return self.__pbc

    @pbc.setter
    def pbc(self, pbc):
        """Set periodic boundary condition for the simulation box."""
        self.__pbc = Error.check_tuple(pbc, bool)

    # Box -----------------------------------------
    @property
    def box_obj(self):
        """Return Box object in StateHeader."""
        return self.__box_obj

    @box_obj.setter
    def box_obj(self, box):
        """Set Box object of StateHeader."""
        if isinstance(box, Box):
            self.__box_obj = box
        else:
            self.__box_obj = Box(box)

    @property
    def box(self):
        """Return box."""
        return self.box_obj.box

    @box.setter
    def box(self, box):
        """Set box."""
        self.box_obj.box = box

    # String representation -------------------------------------
    def __str__(self):
        """String representation of the StateHeader section."""

        out = self.name + "\n\n"
        # out = ''
        out += str(self.simulation_time) + '\n'
        out += str(self.number_of_atoms) + '\n'
        for i in range(3):
            out += str(int(self.box[i])) + ' '
        out += '\n'
        for i in range(3):
            out += str(int(self.pbc[i])) + ' '
        out += '\n'
        return out

    # Operators ------------------------------------
    def __add__(self, other):
        """This method defines '+' between box instances."""

        if not isinstance(other, Box):
            AssertionError("Expected %s for '=' operator!" % self.__class__.__name__)

        # TODO: it picks a box with a larger volume
        if self.box_obj.volume >= other.box_obj.volume:
            return self
        else:
            return other
