#   formatter.py: This file is part of Free Molecular Dynamics
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

""" I/O file formatter"""

from .atom import Atom, Atoms
from .box import Box
from .state import StateHeader

__all__ = []


class Formatter:
    """A base class for implementing a generic file format read and write into molecular frame.
    The required methods have to be implemented in the derived class."""

    def __init__(self, molecular_frame):
        # TODO: handle TypeError for input argument
        self.__formatter_molecular_frame = molecular_frame  # initialize formatter molecular frame

    @property
    def formatter_molecular_frame(self):
        """Returns formatter molecular frame"""
        return self.__formatter_molecular_frame

    # Static method, it is implicitly a class method
    def make(self, format):
        try:
            # Make a derived class using eval, the name for subclass has to start with "Formatter"!
            formatter = eval("Formatter"+str(format))(self.formatter_molecular_frame)

        except (SyntaxError, NameError, TypeError):
            raise AssertionError("Unexpected type for Formatter!")
        return formatter


class FormatterXYZ(Formatter):
    """A derived class reads and writes (via implemented methods) .xyz file format."""

    def write(self, file_name):
        """A method that writes molecular frame into a given file in xyz format."""
        atoms_section = self.formatter_molecular_frame.get_section('Atoms')
        with open(file_name, "w") as out_file:
            out_file.write('%d\n\n'%atoms_section.get_number_of_atoms())
            for atom in atoms_section.atoms:
                out_file.write("%s %f %f %f\n"%(atom.symbol, atom.x, atom.y, atom.z))

    def read(self, file_name, frame=-1):
        """A method that reads specific frame of a *.xyz file and returns molecular frame."""

        n_frame = 0  # initializing frame counter to zero
        atoms = []  # list of read frames

        with open(file_name, 'r') as in_file:

            # Loop over lines in file
            for line in in_file:

                n_atoms = int(line)  # read number of atoms at the begging of the frame
                next(in_file)  # skip one line
                n_frame += 1  # increment frame
                if n_frame >= frame:
                    for index in range(n_atoms):
                        line = next(in_file)
                        line = line.rstrip("/n").split()
                        atom = Atom(symbol=line[0], position=(line[1], line[2], line[3]))
                        # TODO: more info can be read!
                        atoms.append(atom)
                else:
                    # skipping the frame
                    for i in range(n_atoms):
                        next(in_file)

                # check either desired frame or last frame is reached
                if (frame > 0) and (n_frame >= frame):
                    break

        self.formatter_molecular_frame.set_sections(Atoms(atoms))

        # returning sections of created molecular frame
        return self.formatter_molecular_frame.sections


class FormatterSTT(Formatter):
    """A derived class reads and writes (via implemented methods) state (.stt) file format.
    The state format is specific for FMD code."""

    def write(self, file_name):
        """A method that writes molecular frame into .stt (state) format."""

        # TODO: what if there is no such section in MF?
        state = self.formatter_molecular_frame.get_section("StateHeader")
        atoms_section = self.formatter_molecular_frame.get_section('Atoms')

        with open(file_name, "w") as out_file:
            out_file.write("%f\n" % state.simulation_time)
            out_file.write("%d\n" % state.number_of_atoms)
            out_file.write("%f %f %f\n" % state.box)
            out_file.write("%d %d %d\n" % state.pbc)
            for atom in atoms_section.atoms:
                out_file.write("%s %d\n%f %f %f\n%f %f %f\n" % (atom.symbol, atom.group_id, atom.x, atom.y, atom.z,
                                                             atom.vx, atom.vy, atom.vz))

    def read(self, file_name, frame=-1):
        """A method that reads .stt (state) format into molecular frame sections."""

        # Atoms
        atoms_section = Atoms()

        with open(file_name, 'r') as in_file:

            # loop over lines in file
            for line in in_file:
                simulation_time = float(line.split()[0])
                line = next(in_file)
                number_of_atoms = int(line.split()[0])
                line = next(in_file)
                line = line.rstrip("/n").split()
                box = (float(line[0]), float(line[1]), float(line[2]))
                line = next(in_file)
                line = line.rstrip("/n").split()
                # TODO: pbc is boolean!
                pbc = (int(line[0]), int(line[1]), int(line[2]))
                # StateHeader
                state_header = StateHeader(simulation_time, number_of_atoms, box, pbc)
                self.formatter_molecular_frame.set_sections(state_header)
                # Box section
                # TODO: there shoul be a better way!
                self.formatter_molecular_frame.set_sections(Box(box))

                # reading atoms info
                for n in range(number_of_atoms):
                    line = next(in_file)
                    line = line.rstrip("/n").split()
                    symbol, group_id = line[0], int(line[1])
                    line = next(in_file)
                    line = line.rstrip("/n").split()
                    position = (float(line[0]), float(line[1]), float(line[2]))
                    line = next(in_file)
                    line = line.rstrip("/n").split()
                    velocity = (float(line[0]), float(line[1]), float(line[2]))
                    # add an atom to the Atoms
                    atoms_section.append(Atom(symbol=symbol, position=position, group_id=group_id, velocity=velocity))
                break  # stop reading file
        self.formatter_molecular_frame.set_sections(atoms_section)

        # returns all molecular frame's sections
        return self.formatter_molecular_frame.sections



