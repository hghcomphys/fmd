#   unit.py: This file is part of Free Molecular Dynamics
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

"""Unit and physical constants"""

AMU_TO_KILOGRAM = 1.6605E-27
JOULE_TO_EV = 6.242E+18
KELVIN_TO_EV = 8.617332478E-5

KILOGRAM_TO_AMU = 1.0 / AMU_TO_KILOGRAM
EV_TO_JOULE = 1.0 / JOULE_TO_EV
EV_TO_KELVIN = 1.0 / KELVIN_TO_EV
KELVIN_TO_JOULE = KELVIN_TO_EV * EV_TO_JOULE
