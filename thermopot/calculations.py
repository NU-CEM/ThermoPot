import re


class Calculation():
    """Parent class for parsing and storing data from electronic structure calculations."""

    def __init__(self, energy=None, xc=None, NAtoms=None, volume=None, filepath=None):
        """All attributes are None until set by derived classes.

        Attributes:
            volume (float): volume of the periodic unit cell
            filepath (str): path to the calculation output files
            energy (float): DFT total energy in eV
            xc (str): XC functional used to calculate the total energy. Options are "hse06" or "pbesol".
            NAtoms (int): number of atoms in the periodice unit cell

        Returns:
            None.
        """

        self.volume = volume
        self.filepath = filepath
        self.energy = energy
        self.xc = xc
        self.NAtoms = NAtoms

    def check_attributes(self):
        """Check that the Calculation class attributes make basic sense."""

        assert type(self.filepath) == str, "filepath must be a string"
        assert type(self.energy) == float, "energy must be a float"
        assert type(self.xc) == str, "xc must be a string"
        assert (
            type(self.NAtoms) == int
        ) and self.NAtoms >= 1, "NAtoms must be an integer >= 1"
        assert (
            type(self.volume) == float
        ) and self.volume > 0, "volume must be a float > 0"


class AimsCalculation(Calculation):
    """Class for parsing and storing data from a FHI-AIMS total energy calculation."""

    def __init__(self, filepath="./calculation.out"):
        """
        Args:
            filepath: path to the calculation output files
        """
        super().__init__()
        self.filepath = filepath
        self.volume = self.get_volume()
        self.energy = self.get_energy()
        self.xc = self.get_xc()
        self.NAtoms = self.get_NAtoms()

    def get_volume(self):
        with open(self.filepath) as contents:
            return float(
                re.findall("Unit cell volume\s+:\s*(.*)\sA", contents.read())[-1]
            )

    def get_energy(self):
        with open(self.filepath) as contents:
            return float(
                re.findall(
                    "Total energy of the DFT[^0-9]+(-\d*\.?\d*) eV", contents.read()
                )[-1]
            )

    def get_xc(self):
        with open(self.filepath) as contents:
            return re.findall("xc\s+(\S+)", contents.read())[-1]

    def get_NAtoms(self):
        with open(self.filepath) as contents:
            return int(
                re.findall("Number of atoms\s +:\s + (\S+)", contents.read())[-1]
            )
