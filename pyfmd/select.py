#   select.py: This file is part of Free Molecular Dynamics
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

"""Select (atoms)"""

__all__ = ["Select", "SelectRegion", "SelectGroupId"]


class SelectBase:
    """A base class for generic section of atoms."""

    def is_selected(self, atom):
        """An abstract method that has to be defined in any Select derived class."""
        return True

    def AND(self, other):
        """Implement AND between two Select instances."""
        return SelectOperator(self, other, lambda x, y: x and y)

    def OR(self, other):
        """Implement OR between two Select instances."""
        return SelectOperator(self, other, lambda x, y: x or y)

    def NOT(self):
        """Implement NOT for a Select instance."""
        return SelectNOT(self)


class Select:
    """A base class for an easier interface of making select objects."""

    def __init__(self):
        pass

    def region(self, region_fn):
        """Make a SelectRegion instance."""
        return SelectRegion(region_fn)

    def group_id(self, group_id):
        """Make a SelectTag instance."""
        return SelectGroupId(group_id)


class SelectOperator(SelectBase):
    """A derived class that defines logical operator between two select instances."""

    def __init__(self, sel_obj1, sel_obj2, operator):
        self.sel_obj1 = sel_obj1
        self.sel_obj2 = sel_obj2
        self.operator = operator

    def is_selected(self, atom):
        return self.operator(self.sel_obj1.is_selected(atom), self.sel_obj2.is_selected(atom))


class SelectNOT(SelectBase):
    """A class that invert the selection within a molecular frame."""

    def __init__(self, sel_obj):
        self.sel_obj = sel_obj

    def is_selected(self, atom):
        return not self.sel_obj.is_selected(atom)


class SelectRegion(SelectBase):
    """This class selects atoms within a specified region."""

    def __init__(self, region_fn=None):
        """Region function has to be in form of f(x,y,z) where input arguments refer to coordinates of an atom."""
        if region_fn is None:
            raise AssertionError("Expected region function f(x,y,z) for %s!" % self.__class__.__name__)

        # Set a region function
        self.__region_fn = region_fn

    def is_selected(self, atom):
        """Define criteria for selecting atoms."""
        try:
            if bool(self.__region_fn(atom.x, atom.y, atom.z)):
                return True
            else:
                return False
        except TypeError:
            raise


class SelectGroupId(SelectBase):
    """A class that selects atoms with certain group-id(s)."""

    def __init__(self, group_id=None):
        self.__group_id = []  # make a empty list of group-id

        try:
            if isinstance(group_id, list):
                for each in group_id:
                    self.__group_id.append(int(each))
            else:
                self.__group_id.append(int(group_id))

        except TypeError:
            raise AssertionError("Expected a group-id for %s!" % self.__class__.__name__)

    def is_selected(self, atom):
        """Define a criteria for selecting the atoms."""

        if atom.group_id in self.__group_id:
            return True
        else:
            return False
