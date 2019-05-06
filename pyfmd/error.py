#   error.py: This file is part of Free Molecular Dynamics
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

"""This script contains functions for a unified error handling of classes mainly through getters and setters."""

__all__ = []


class Error(Exception):
    """A Base class that provides static methods for errors handling."""

    @staticmethod
    def ge_zero(x, func, msg=None):
        """A function for unexpected negative value of generic input argument."""
        if msg is None:
            msg = "Unexpected negative value!"
        x = func(x)
        if x < 0:
            raise ValueError(msg)
        return x

    @staticmethod
    def int_ge_zero(n, msg=None):
        """A function for unexpected negative value of integer input argument."""
        return Error.ge_zero(n, int, msg)

    @staticmethod
    def float_ge_zero(x, msg=None):
        """A function for unexpected negative value of floating-point input argument."""
        return Error.ge_zero(x, float, msg)

    @staticmethod
    def gt_zero(x, func, msg=None):
        """A function for unexpected negative or zero value of generic input argument."""
        if msg is None:
            msg = "Unexpected negative or zero value!"
        x = func(x)
        if x <= 0:
            raise ValueError(msg)
        return x

    @staticmethod
    def int_gt_zero(n, msg=None):
        """A function for unexpected negative or zero value of integer input argument."""
        return Error.gt_zero(n, int, msg)

    @staticmethod
    def float_gt_zero(x, msg=None):
        """A function for unexpected negative or zero value of float input argument."""
        return Error.gt_zero(x, float, msg)

    @staticmethod
    def check_tuple(x, func, size=3, msg=None):
        """A method that checks type and size of generic input tuple argument."""
        if msg is None:
            msg = "Unexpected tuple!"
        try:
            if not len(x) == size:
                raise ValueError("Unexpected tuple size! (%d instead of %d)" % (len(x), size))
            return tuple([func(_) for _ in x])

        except (TypeError, ValueError):
            raise ValueError(msg)

    @staticmethod
    def check_tuple_float(x, size=3, msg=None):
        """A method that checks type and size of floating-point input tuple argument."""
        return Error.check_tuple(x, float, size, msg)

    @staticmethod
    def check_tuple_int(x, size=3, msg=None):
        """A method that checks type and size of integer input tuple argument."""
        return Error.check_tuple(x, float, size, msg)


