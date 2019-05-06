#   section.py: This file is part of Free Molecular Dynamics
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

"""Section classes for Molecular Frame"""

__all__ = []


class SectionName:
    """A class that holds name for a section."""

    def __init__(self, name):
        self.name = name  # setting section's name

    @property
    def name(self):
        """This method returns given name of the section."""
        return self.__name

    @name.setter
    def name(self, name):
        """This method sets a given name for the section."""
        self.__name = str(name)

    def get_name(self):
        return self.name

    def set_name(self, name):
        self.name = name
        return self


class FixedSection(SectionName):
    """A base class that contains a single item such as simulation box."""

    def __init__(self, name):
        SectionName.__init__(self, name)


class ListSection(SectionName):
    """A base class that contains a list of items such as list of atoms, bonds, etc."""

    def __init__(self, name):
        SectionName.__init__(self, name)
        self.__items = []  # empty list of items

    # Items -------------------------------------------------
    @property
    def items(self):
        """This method returns list of items in the list section."""
        return self.__items

    @property
    def number_of_items(self):
        """This method returns number of items in the list."""
        return len(self.items)

    # String representation -------------------------------------
    def __str__(self):
        """This method defines string representation of the list section."""
        out = self.name + "\n\n"
        for item in self.items:
            out += str(item) + '\n'
        return out

    # Operators ----------------------------------------------
    def __eq__(self, other):
        """This method defines equal '=' operator between two list sections."""

        if not isinstance(other, ListSection):
            raise AssertionError("Expected %s for '=' operator!" % self.__class__.__name__)

        self.__name = other.name
        self.__items = other.items
        return self

    def __add__(self, other):
        """A dummy method for adding two box objects."""
        return self
