#   box.py: This file is part of Free Molecular Dynamics
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

"""Box"""

from .section import FixedSection
from .error import Error

__all__ = ['Box']


class Box(FixedSection):
    """This derived class contains data for simulation box and relevant methods.
    It includes box sizes, volume, etc."""

    def __init__(self, box=(0.0, 0.0, 0.0)): # xy = 0.0, xz = 0.0, yz = 0.0
        """Initialize box"""
        FixedSection.__init__(self, 'Box')  # set name of the section
        try:
            self.box = box # set box sizes
            # TODO: non-orthogonal not implemented yet!

        except ValueError:
            raise AssertionError("Unexpected value for %s!" % self.__class__.__name__)

    # Lx -----------------------------------
    @property
    def xlo(self):
        return self.__xlo

    @xlo.setter
    def xlo(self, xlo):
        self.__xlo = float(xlo)

    @property
    def xhi(self):
        return self.__xhi

    @xhi.setter
    def xhi(self, xhi):
        self.__xhi = float(xhi)

    @property
    def lx(self):
        return self.xhi - self.xlo

    @lx.setter
    def lx(self, lx):
        self.xlo = 0.0
        self.xhi = Error.float_ge_zero(lx)

    # Ly -----------------------------------
    @property
    def ylo(self):
        return self.__ylo

    @ylo.setter
    def ylo(self, ylo):
        self.__ylo = float(ylo)

    @property
    def yhi(self):
        return self.__yhi

    @yhi.setter
    def yhi(self, yhi):
        self.__yhi = float(yhi)

    @property
    def ly(self):
        return self.yhi - self.ylo

    @ly.setter
    def ly(self, ly):
        self.ylo = 0.0
        self.yhi = Error.float_ge_zero(ly)

    # Lz -----------------------------------
    @property
    def zlo(self):
        return self.__zlo

    @zlo.setter
    def zlo(self, zlo):
        self.__zlo = float(zlo)

    @property
    def zhi(self):
        return self.__zhi

    @zhi.setter
    def zhi(self, zhi):
        self.__zhi = float(zhi)

    @property
    def lz(self):
        return self.zhi - self.zlo

    @lz.setter
    def lz(self, lz):
        self.zlo = 0.0
        self.zhi = Error.float_ge_zero(lz)

    # Box -------------------------------------
    @property
    def box(self):
        return self.lx, self.ly, self.lz

    @box.setter
    def box(self, box=(0.0, 0.0, 0.0)):
        self.lx, self.ly, self.lz = Error.check_tuple_float(box)

    def get_box(self):
        return self.box

    def set_box(self, box=(0.0, 0.0, 0.0)):
        self.box = box
        return self

    def translate(self, translation_vec=(0.0, 0.0, 0.0)):
        """Translate box along a given vector."""
        vec = Error.check_tuple_float(translation_vec)
        self.xlo += vec[0]
        self.xhi += vec[0]
        self.ylo += vec[1]
        self.yhi += vec[1]
        self.zlo += vec[2]
        self.zhi += vec[2]
        return self

    # String representation -----------------------
    def __str__(self):
        """String representation of the box section."""
        out = self.name + '\n\n'
        for lo, hi, l in zip([self.xlo, self.ylo, self.zlo],
                             [self.xhi, self.yhi, self.zhi], [self.lx, self.ly, self.lz]):
            out += str(lo) + ' ' + str(hi) + ' ' + str(l)
            out += '\n'
        return out

    # Volume ---------------------------------------
    @property
    def volume(self):
        """A method that returns box volume."""
        # TODO: only works for orthogonal box
        return self.lx * self.ly * self.lz

    def get_volume(self):
        return self.volume

    # Operators ------------------------------------
    def __add__(self, other):
        """This method defines '+' between box instances."""

        if not isinstance(other, Box):
            AssertionError("Expected %s for '=' operator!" % self.__class__.__name__)

        # TODO: it picks a box with a larger volume
        if self.volume >= other.volume:
            return self
        else:
            return other
