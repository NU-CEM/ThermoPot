"""
Module to parse and store data from electronic structure calculations.
Contains the parent class Calculation to store data from a variety of sources.
Contains the child class AimsCalculation to read and store data from a FHI-aims calculation.
"""

import re


class Calculation:
    """
    Parent class for parsing and storing data from electronic structure calculations.

    Example:

        BaS_calc = Calculation(volume=63.2552, energy=-235926.586148547, xc='pbesol', num_atoms=2)

    Attributes:

        volume (float): volume of the periodic unit cell in Angstrom^3
        filepath (str): path to the calculation output files
        energy (float): DFT total energy in eV
        xc (str): XC functional used to calculate the total energy
        num_atoms (int): number of atoms in the periodic unit cell

    Note:

        If gas is True then no volume attribute is required.
    """

    def __init__(
        self,
        energy=None,
        xc=None,
        num_atoms=None,
        volume=None,
        filepath=None,
        gas=False,
    ):
        """
        Note:

            All attributes are None until set by derived classes or specified by user.

        Args:

            volume (float): volume of the periodic unit cell in Angstrom^3
            filepath (str, optional): path to the calculation output files
            energy (float): DFT total energy in eV
            xc (str): XC functional used to calculate the total energy
            num_atoms (int): number of atoms in the periodic unit cell
            gas (bool): True if gas species, False otherwise

        Note:

            If gas is True then volume is None
        """
        if not gas:
            self.volume = volume
        self.filepath = filepath
        self.energy = energy
        self.xc = xc
        self.num_atoms = num_atoms

        # self.check_attributes()

    def check_attributes(self):
        """Check that the Calculation class attributes make basic sense."""

        assert (
            type(self.filepath) == str or self.filepath is None
        ), "filepath must be a string"
        assert type(self.energy) == float, "energy must be a float"
        assert type(self.xc) == str, "xc must be a string"
        assert (
            type(self.num_atoms) == int
        ) and self.num_atoms >= 1, "num_atoms must be an integer >= 1"
        assert (
            type(self.volume) == float
        ) and self.volume > 0, "volume must be a float > 0"


class AimsCalculation(Calculation):
    """Class for parsing and storing data from a FHI-AIMS total energy calculation.

    Example:

       BaS_calc = AimsCalculation("./aims_output/output.aims")

    Attributes:

        volume (float): volume of the periodic unit cell in Angstrom^3
        filepath (str): path to the calculation output files
        energy (float): DFT total energy in eV
        xc (str): XC functional used to calculate the total energy
        num_atoms (int): number of atoms in the periodic unit cell
    """

    def __init__(self, filepath="./calculation.out", gas=False):
        """
        Args:

            filepath (str): path to the calculation output files
            gas (bool): True if gas species, False otherwise

        Note:

            If gas is True then volume is None
        """
        super().__init__()
        self.filepath = filepath
        if not gas:
            self.volume = self.get_volume()
        self.energy = self.get_energy()
        self.xc = self.get_xc()
        self.num_atoms = self.get_num_atoms()

    def get_volume(self):
        """
        Returns:

            (float): volume of the periodic unit cell in Angstrom^3
        """
        with open(self.filepath) as contents:
            return float(
                re.findall("Unit cell volume\s+:\s*(.*)\sA", contents.read())[-1]
            )

    def get_energy(self):
        """
        Returns:

            (float): DFT total energy in eV
        """
        with open(self.filepath) as contents:
            return float(
                re.findall(
                    "Total energy of the DFT[^0-9]+(-\d*\.?\d*) eV", contents.read()
                )[-1]
            )

    def get_xc(self):
        """
        Returns:

            (str): XC functional used to calculate the total energy
        """
        with open(self.filepath) as contents:
            return re.findall("xc               (\S+)", contents.read())[0]

    def get_num_atoms(self):
        """
        Returns:

            (int): number of atoms in the periodic unit cell
        """
        with open(self.filepath) as contents:
            return int(
                re.findall("Number of atoms\s +:\s + (\S+)", contents.read())[-1]
            )


class QHACalculation:
    def __init__(self, ev_filepath, filepath):
        self.ev_filepath = ev_filepath
        self.filepath = filepath
        self.volumes = self.get_volumes()
        self.energies = self.get_energies()
        self.xc = self.get_xc()
        self.num_atoms = self.get_num_atoms()

    def get_volumes(self):
        volumes = []
        with open(self.ev_filepath) as contents:
            for line in contents:
                line = line.split()
                volumes.append(float(line[0]))
        return volumes

    def get_energies(self):
        energies = []
        with open(self.ev_filepath) as contents:
            for line in contents:
                line = line.split()
                energies.append(float(line[1]))
        return energies

    def get_xc(self):
        """
        Returns:

            (str): XC functional used to calculate the total energy
        """
        with open(self.filepath) as contents:
            return re.findall("xc               (\S+)", contents.read())[0]

    def get_num_atoms(self):
        """
        Returns:

            (int): number of atoms in the periodic unit cell
        """
        with open(self.filepath) as contents:
            return int(
                re.findall("Number of atoms\s +:\s + (\S+)", contents.read())[-1]
            )
