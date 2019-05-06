#   adaptor.py: This file is part of Free Molecular Dynamics
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

"""Adaptor for an external package"""

from ase import Atom as AtomASE
from ase import Atoms as AtomsASE
import numpy as np
# TODO: import ASE module error handling

from .atom import Atom, Atoms
from .box import Box

__all__ = []


class Adaptor:
    """A generic class adaptor for working with an external package."""

    def __init__(self, package_instance=None):
        self.__package_instance = package_instance  # external package object including structural info
        self.__adaptor_sections = dict()  # sections of the molecular frame

    @property
    def package_instance(self):
        """Return external package name"""
        return self.__package_instance

    @package_instance.setter
    def package_instance(self, package_instance):
        """Set external package instance"""
        self.__package_instance = package_instance

    @property
    def adaptor_sections(self):
        """Return adaptor molecular sections"""
        return self.__adaptor_sections

    def make(self, package_name):
        """A method that creates an appropriate adaptor instance."""
        try:
            # Make a drived class using eval, the name for subclass has to start with "Adaptor"!
            formatter = eval("Adaptor" + str(package_name))(self.package_instance)

        except (SyntaxError, NameError, TypeError):
            raise AssertionError("Unexpected type for Formatter!")

        # Return the created instance of the adaptor
        return formatter


class AdaptorASE(Adaptor):
    """Adaptor for ASE package that imports and exports the molecular frame (structure)."""

    def get_sections(self):
        """A method that returns sections imported though ASE adaptor."""

        # Box section
        cell_info = self.package_instance.get_cell().diagonal()
        # TODO: only works for orthogonal cell!
        # TODO: what if package instance has not box info?
        box_section = Box(cell_info)
        self.adaptor_sections[box_section.name] = box_section

        # Atoms section
        atom_positions = self.package_instance.get_positions()
        atom_labels = self.package_instance.get_chemical_symbols()
        # TODO: what if package instance has not atoms/symbol info?
        atoms_section = Atoms()
        for aid, pos, symbol in zip(range(len(atom_labels)), atom_positions, atom_labels):
            atom = Atom(symbol=symbol, position=(pos[0], pos[1], pos[2]))
            # TODO: more info can be added!
            atoms_section.append(atom)
            self.adaptor_sections[atoms_section.name] = atoms_section

        # Return the created dictionary of molecular frame sections
        return self.adaptor_sections

    def set_sections(self, molecular_frame):
        """This method adds atoms from a molecular frame to the package instance"""

        self.package_instance = AtomsASE()
        # Add atoms
        if Atoms().name in molecular_frame.get_sections_name():
            # TODO: limited to .xyz format!
            for atom in molecular_frame.atoms_section.atoms:
                self.package_instance.append(AtomASE(atom.symbol,(atom.x, atom.y, atom.z)))
        # Add box
        if Box().name in molecular_frame.get_sections_name():
            box = molecular_frame.box_section
            # TODO: limited to orthogonal box!
            cell = np.array([[box.lx, 0., 0.],
                            [0., box.ly, 0.],
                            [0., 0., box.lz]])
            self.package_instance.set_cell(cell)

        # Returns the created package instance
        return self.package_instance
