#   molframe.py: This file is part of Free Molecular Dynamics
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

""" Molecular Frame"""

from copy import deepcopy

from .section import ListSection, FixedSection
from .adaptor import Adaptor
from .formatter import Formatter
from .state import StateHeader

__all__ = ['MolecularFrame']


class MolecularFrame:
    """A class for preparing molecular structure that contains molecular topologies and simulation box.
    It contains in fact a collection of sections that basically includes list of atoms, box, masses, etc.
    It has several features such as read and write files in several formats (such as .xyz or .stt),
    import (export) molecular structure from (to) an external packages (i.e. ASE),
    selecting specific group of atoms, and integrating two molecular frames."""

    def __init__(self, name='Molecular Frame', sections=None):
        """Make a new instance of molecular frame."""
        self.name = name  # assign a given name for molecular frame
        self.__sections = {}  # initialize sections as a dictionary

        if sections is None:
            pass
        else:
            self.set_sections(sections)

    @property
    def name(self):
        """Returns molecular frame name."""
        return self.__name

    # Name ---------------------------------------
    @name.setter
    def name(self, name):
        self.__name = str(name)

    def get_name(self):
        return self.name

    def set_name(self, name):
        self.name = name
        return self

    # Sections -------------------------------------
    @property
    def sections(self):
        """This method returns all sections in form of a dictionary."""
        return self.__sections

    @sections.setter
    def sections(self, sections):
        """This method sets given molecular sections."""
        if isinstance(sections, dict):
            for sec in sections.values():
                self.__set_section(sec)
        elif isinstance(sections, list):
            for sec in sections:
                self.__set_section(sec)
        else:
            self.__set_section(sections)

    def __set_section(self, section):
        """This method sets a specified molecular section into molecular frame."""
        if not (isinstance(section, ListSection) or isinstance(section, FixedSection)):
            raise AssertionError("Unexpected given molecular section for %s!" % self.__class__.__name__)
        self.__sections[section.name] = section

    def get_sections(self):
        return self.sections

    def set_sections(self, sections):
        self.sections = sections
        return self

    def set_section(self, section):
        self.sections = section
        return self

    @property
    def sections_name(self):
        """Returns list of name of sections."""
        return self.sections.keys()

    def get_sections_name(self):
        """Returns list of name of sections."""
        return self.sections_name

    def get_section(self, section_name):
        """A method that returns an specified molecular section by giving its name (key)."""
        if section_name not in self.sections.keys():
            raise AssertionError("Cannot find section name for %s!" % self.get_section.__name__)
        return self.sections[section_name]

    # String representation ------------------------------
    def __str__(self):
        """A method that defines string conversion for a molecular frame instance."""
        out = self.__name + '\n\n'
        for mol_sec in self.sections.values():
            out += str(mol_sec) + '\n'
        return out

    # Clone ------------------------------------------
    def clone(self):
        """A method that makes a deep copy of the molecular frame."""
        return deepcopy(self)

    # Clear ------------------------------------------
    def clear(self):
        """A method that clears sections of the molecular frame."""
        self.sections.clear()
        return self

    # Atoms ------------------------------------------
    @property
    def atoms_section(self):
        """Returns list of atoms"""
        return self.get_section("Atoms")

    def get_atoms_section(self):
        """Return atom sections"""
        return self.atoms_section

    @property
    def number_of_atoms(self):
        """Returns the number of atoms"""
        return self.atoms_section.number_of_atoms

    def get_number_of_atoms(self):
        return self.number_of_atoms

    def translate_atoms(self, translation_vec=(0.0, 0.0, 0.0), inplace=True):
        """Translates all atoms along a given direction."""

        if inplace:
            self_new = self  # apply inplace translation
        else:
            self_new = self.clone()  # make a deep copy the self and then applying the translation

        self_new.atoms_section.translate(translation_vec)
        return self_new

    def add_velocities(self, volocity_vec=(0.0, 0.0, 0.0), inplace=True):
        """Add a constant value to all atoms's velocities."""

        if inplace:
            self_new = self  # apply inplace velocity translation
        else:
            self_new = self.clone()  # make a deep copy the self and then applying the velocity translation

        self_new.atoms_section.add_velocities(volocity_vec)
        return self_new

    def set_group_id(self, group_id):
        self.atoms_section.set_group_id(group_id)
        return self

    # Box ----------------------------------------------
    @property
    def box_section(self):
        """Returns box section."""
        return self.get_section("Box")

    def get_box_section(self):
        """Returns box section."""
        return self.box_section

    def translate_box(self, translation_vec=(0.0, 0.0, 0.0)):
        """A method that translate box, it has no effect on box sizes and only changes box min/max."""
        self.box_section.translate(translation_vec)

    # Formatter ----------------------------------------
    def read(self, file_name, file_format=None):
        """This method reads a file with specified format into the molecular frame."""
        if file_format is None:
            file_format = file_name.split('.')[-1]
        self.clear().set_sections(Formatter(self).make(file_format.upper()).read(file_name))
        return self

    def write(self, file_name, file_format=None):
        """A method that writes the molecular frame into a file with given format."""
        if file_format is None:
            file_format = file_name.split('.')[-1]
        Formatter(self).make(file_format.upper()).write(file_name)
        # return self

    # Adaptor --------------------------------------------
    def import_from(self, package_instance, package_name="ASE"):
        """This method imports the molecular frame from an external module (i.e. ASE)."""
        # TODO: default package name is set to "ASE"!
        # TODO: error handling for a given input package
        self.set_sections(Adaptor(package_instance).make(package_name).get_sections())
        return self

    def export(self, package_name="ASE"):
        """This method exports the  molecular frame as an specified external package instance (i.e. ASE)."""
        return Adaptor().make(package_name).set_sections(self)

    # Perators  ---------------------------------------
    def __eq__(self, other):
        """This method defines equal '=' operator between two molecular frame instances."""
        if not isinstance(other, MolecularFrame):
            raise AssertionError("Expected %s for '=' operator!" % self.__class__.__name__)

        # Set name and molecular sections from 'other' instance
        self.name = other.name
        self.set_sections(other.sections)
        return self

    def __add__(self, other):
        """This method defines '+' operator between two molecular frame instances."""
        if not isinstance(other, MolecularFrame):
            raise AssertionError("Expected %s for '+' operator!" % self.__class__.__name__)

        # make a new molecular frame by integrating two given ones
        new_mf = MolecularFrame(self.name + ' + ' + other.name)
        # TODO: it only picks the shared sections
        for sec_name in self.sections_name:
            if sec_name in other.sections_name:
                new_mf.set_sections(self.get_section(sec_name) + other.get_section(sec_name))
        return new_mf

    # Select --------------------------------------------
    def select(self, select_obj=None):
        """This method returns a molecular frame determined the Select object."""
        mf_sel = MolecularFrame(self.name, self.sections)
        mf_sel.set_section(self.atoms_section.select(select_obj))
        # TODO: there is a better way to do it (observer pattern design)!
        if StateHeader().name in self.sections_name:
            self.get_section(StateHeader().name).number_of_atoms = mf_sel.number_of_atoms
        return mf_sel

    # Center of mass ------------------------------------
    @property
    def center_of_mass(self):
        """This method returns the center of mass vector position."""
        return self.atoms_section.center_of_mass

    @center_of_mass.setter
    def center_of_mass(self, position):
        """This method translates the center of mass to the given vector position."""
        self.atoms_section.center_of_mass = position

    @property
    def center_of_mass_velocity(self):
        """This method returns center of mass vector."""
        return self.atoms_section.center_of_mass_velocity

    @center_of_mass_velocity.setter
    def center_of_mass_velocity(self, velocity):
        """This method translates the center of mass to a given vector position."""
        self.atoms_section.center_of_mass_velocity = velocity

    # Temperature ----------------------------------------
    def get_temperature(self):
        """A method that returns the temperature (in Kelvin)."""
        return self.atoms_section.temperature

    @property
    def temperature(self):
        """A method that returns the temperature (in Kelvin)."""
        return self.get_temperature()

    def set_temperature(self, temperature, seed=None):
        """Set the temperature via assigning random velocities to the atoms."""
        self.atoms_section.set_temperature(temperature, seed)

    @temperature.setter
    def temperature(self, temperature):
        """Set the temperature via assigning random velocities to the atoms."""
        self.set_temperature(temperature)


