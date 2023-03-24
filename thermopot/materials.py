"""
Module contains the classes Solid and IdealGas to store basic material data.
Each class provides methods for calculating various thermodynamic properties.
"""

import numpy as np
from scipy import constants
from thermopot import interpolate


import os  # get correct path for datafiles when called from another directory

materials_directory = os.path.dirname(__file__)
# Append a trailing slash to make coherent directory name - this would select the
# root directory in the case of no prefix, so we need to check
if materials_directory:
    materials_directory = materials_directory + "/"


class Material(object):
    """
    Parent class for storing materials properties.

    Attributes:

        name (str): Identifying string
        stoichiometry (dict): relates element to the number of atoms in a single formula unit
        energies (dict): relates xc functional to DFT total energy in eV
        N (int): number of atoms per formula unit
    """

    def __init__(self, name, stoichiometry, energies):
        self.name = name
        self.stoichiometry = stoichiometry
        self.energies = energies
        self.N = sum(self.stoichiometry.values())


class Solid(Material):
    """
    Class for solid material data.

    Note:

        The material is assumed to be incompressible and without thermal expansion.

    Example:

        BaS = Solid('BaS',{'Ba':1,'S':1},"./phonon_data/Ba_S",calculation=BaS_calc)

    Attributes:

        name (str): Identifying string
        stoichiometry (dict): relates element to the number of atoms in a single formula unit
        energies (dict): relates xc functional to DFT total energy in eV
        N (int): number of atoms per formula unit
        fu_cell (int): number of formula units in periodic unit cell
        volume (float): volume of unit cell in Angstroms^3
        phonon_filepath (str): path to the phonon output data
        NAtoms (int): number of atoms in periodic unit cell
    """

    def __init__(
        self,
        name,
        stoichiometry,
        phonon_filepath,
        calculation=False,
        QHACalculation=False,
        volume=False,
        energies=False,
        NAtoms=1,
    ):
        """
        Args:

           name (str): Identifying string
           stoichiometry (dict): relates element to the number of atoms in a single formula unit
           phonon_filepath (str): path to the phonon output data
           calculation (thermopot.calculation.Calculation, optional): instance of the thermopot.calculation.Calculation class
           volume (float, optional): volume of unit cell in Angstroms^3
           energies (dict, optional): relates xc functional to DFT total energy in eV
           NAtoms (int): number of atoms in periodic unit cell
        """
        if QHACalculation is not False:
            Material.__init__(
                self, name, stoichiometry, {QHACalculation.xc: QHACalculation.energies}
            )
            self.volume = QHACalculation.volumes
            self.NAtoms = QHACalculation.NAtoms

        elif calculation is not False:
            if type(calculation) is not list:
                Material.__init__(
                    self, name, stoichiometry, {calculation.xc: calculation.energy}
                )
                self.volume = calculation.volume
                self.NAtoms = calculation.NAtoms
            else:
                pass

        else:
            Material.__init__(self, name, stoichiometry, energies)

            self.NAtoms = NAtoms
            self.volume = volume

        self.fu_cell = self.NAtoms / self.N
        self.phonon_filepath = materials_directory + phonon_filepath

    def V(self, T, xc="pbesol", units="eV"):
        V_func = interpolate.get_potential_aims(self.phonon_filepath, "V")
        return V_func

    def U(self, T, xc="pbesol", units="eV"):
        """
        Calculates the internal energy of one formula unit of solid.

        Example:

            U = BaS.U(300,xc="pbesol",units="eV")

        Args:

            T (float/ndarray): 1D Numpy array containing temperature data (in Kelvin) as floats, or a single temperature as a float.
            xc (str, optional): DFT XC functional used to calculate the ground state energy
            units (str, optional):  specifies the units as "eV", "J" (J/mol) or "kJ" (kJ/mol)

        Returns:

            U (float/ndarray): 1D Numpy array (with the same dimensions as T) containing the internal energies of one formula unit of solid, or a single internal energy float when a single temperature is passed as an argument.

        """
        U_func = interpolate.get_potential_aims(self.phonon_filepath, "U")
        E_dft = self.energies[xc]

        U_eV = (E_dft + U_func(T)) / self.fu_cell

        if units == "eV":
            return U_eV

        elif units == "J":
            return (
                U_eV(T, xc=xc)
                * constants.physical_constants["electron volt-joule relationship"][0]
                * constants.N_A
            )

        elif units == "kJ":
            return (
                U_eV(T, xc=xc)
                * constants.physical_constants["electron volt-joule relationship"][0]
                * constants.N_A
                * 0.001
            )

    def H(self, T, P, xc="pbesol", units="eV"):
        """
        Calculates the Enthalpy (H = U + PV) of one formula unit of solid.

        Examples:

            H = BaS.H(300,1E3,xc="pbesol",units="eV")
            H = BaS.H(np.linspace(100,700,1000),np.array(np.logspace(1, 7, 100),ndmin=2).transpose())

        Args:

            T (float/ndarray): 1D Numpy array containing temperature data (in Kelvin) as floats, or a single temperature as a float.
            P (float/ndarray): 2D Numpy array with a single row containing pressure data (in Pa) as floats, or a single pressure as a float.
            xc (str, optional): DFT XC functional used to calculate the ground state energy
            units (str, optional):  specifies the units as "eV", "J" (J/mol) or "kJ" (kJ/mol)

        Note:

            T, P are orthogonal 2D arrays of length m and n, populated in one row/column: in this case H is an m x n matrix.
            Other T, P arrays will result in undefined behaviour.

        Returns:

            H (float/ndarray): Enthalpy of one formula unit of solid expressed as floats in a m x n Numpy array where T, P are orthogonal 2D arrays of length m and n
        """
        U_func = interpolate.get_potential_aims(self.phonon_filepath, "U")
        PV = (
            P
            * self.volume
            * 1e-30
            * constants.physical_constants["joule-electron volt relationship"][0]
            / constants.N_A
        )
        E_dft = self.energies[xc]
        H_eV = ((E_dft + U_func(T)) + PV) / self.fu_cell

        if units == "eV":
            return H_eV

        elif units == "J":
            return (
                H_eV
                * constants.physical_constants["electron volt-joule relationship"][0]
                * constants.N_A
            )

        elif units == "kJ":
            return (
                H_eV
                * constants.physical_constants["electron volt-joule relationship"][0]
                * constants.N_A
                * 0.001
            )

    def mu(self, T, P, xc="pbesol", units="eV"):
        """
        Calculates the Gibbs Free Energy (mu = U + PV - TS) of one formula unit of solid.

        Examples:

            mu = BaS.mu(300,1E3,xc="pbesol",units="eV")
            mu = BaS.mu(np.linspace(100,700,1000),np.array(np.logspace(1, 7, 100),ndmin=2).transpose())

        Args:

            T (float/ndarray): 1D Numpy array containing temperature data (in Kelvin) as floats, or a single temperature as a float.
            P (float/ndarray): 2D Numpy array with a single row containing pressure data (in Pa) as floats, or a single pressure as a float.
            xc (str, optional): DFT XC functional used to calculate the ground state energy
            units (str, optional):  specifies the units as "eV", "J" (J/mol) or "kJ" (kJ/mol)

        Note:

            T, P are orthogonal 2D arrays of length m and n, populated in one row/column: in this case mu is an m x n matrix.
            Other T, P arrays will result in undefined behaviour.

        Returns:

            mu (float/ndarray): Gibbs Free Energy of one formula unit of solid expressed as floats in a m x n Numpy array where T, P are orthogonal 2D arrays of length m and n
        """
        TS_func = interpolate.get_potential_aims(self.phonon_filepath, "TS")
        H = self.H(T, P, xc=xc)
        mu_eV = H - (TS_func(T)) / self.fu_cell

        if units == "eV":
            return mu_eV

        elif units == "J":
            return (
                mu_eV
                * constants.physical_constants["electron volt-joule relationship"][0]
                * constants.N_A
            )

        elif units == "kJ":
            return (
                mu_eV
                * constants.physical_constants["electron volt-joule relationship"][0]
                * constants.N_A
                * 0.001
            )

    def Cv(self, T, units="kB"):
        """
        Calculates the Constant-volume heat capacity of one formula unit of solid.

        Examples:

            Cv = BaS.mu(300,xc="pbesol",units="eV")
            Cv = BaS.mu(np.linspace(100,700,1000),xc="pbesol",units="kJ")

        Args:

            T (float/ndarray): 1D Numpy array containing temperature data (in Kelvin) as floats, or a single temperature as a float.
            xc (str, optional): DFT XC functional used to calculate the ground state energy
            units (str, optional):  specifies the units as "eV", "J" (J/mol) or "kJ" (kJ/mol)

        Returns:

            Cv (float/ndarray): 1D Numpy array (with the same dimensions as T) containing the Constant-volume heat capacity of one formula unit of solid, or a single heat capacity float when a single temperature is passed as an argument.
        """
        Cv_func = interpolate.get_potential_aims(self.phonon_filepath, "Cv")
        Cv_kB = Cv_func(T) / self.fu_cell

        if units == "kB":
            return Cv_kB

        elif units == "eV":
            return Cv_kB * constants.physical_constants["Boltzmann constant in eV/K"][0]

        elif units == "J":
            return (
                Cv_kB
                * constants.physical_constants["Boltzmann constant"][0]
                * constants.N_A
            )

        elif units == "kJ":
            return (
                Cv_kB
                * constants.physical_constants["Boltzmann constant"][0]
                * constants.N_A
                * 0.001
            )


class IdealGas(Material):
    """
    Class for ideal gas properties.

    Example:

        S2_gas = materials.IdealGas("S2", {"S":2}, "./thermo_data/S2",calculation=S2_gas_calc)

    Attributes:

       name (str): Identifying string
       stoichiometry (dict): relates element to the number of atoms in a single formula unit
       thermo_dile (str): path to the thermodynamics data
       calculation (thermopot.calculation.Calculation, optional): instance of the thermopot.calculation.Calculation class
       energies (dict, optional): relates xc functional to DFT total energy in eV
       zpe_pbesol (float, optional): zero point energy calculated using the pbesol XC-functional
       zpe_hse06 (float, optional): zero point energy calculated using the hse06 XC-functional
       zpe_lit (float, optional): zero point energy calculated using literature values

    Note:

        Ideal gas law PV=nRT is applied: specifically (dH/dP) at const. T = 0 and int(mu)^P2_P1 dP = kTln(P2/P1).
        Enthalpy has no P dependence as volume is not restricted / expansion step is defined as isothermal
    """

    def __init__(
        self,
        name,
        stoichiometry,
        thermo_file,
        calculation=False,
        energies=False,
        zpe_pbesol=0,
        zpe_hse06=0,
        zpe_lit=0,
    ):
        """
        Args:

        name (str): Identifying string
        stoichiometry (dict): relates element to the number of atoms in a single formula unit
        thermo_dile (str): path to the thermodynamics data
        calculation (thermopot.calculation.Calculation, optional): instance of the thermopot.calculation.Calculation class
        energies (dict, optional): relates xc functional to DFT total energy in eV
        zpe_pbesol (float, optional): zero point energy calculated using the pbesol XC-functional
        zpe_hse06 (float, optional): zero point energy calculated using the hse06 XC-functional
        zpe_lit (float, optional): zero point energy calculated using literature values
        """
        if calculation is not False:
            Material.__init__(
                self, name, stoichiometry, {calculation.xc: calculation.energy}
            )
        else:
            Material.__init__(self, name, stoichiometry, energies)
        self.thermo_file = materials_directory + thermo_file

        if zpe_hse06 > 0:
            self.zpe = zpe_pbesol
        elif zpe_pbesol > 0:
            self.zpe = zpe_pbesol
        elif zpe_lit > 0:
            self.zpe = zpe_lit
        else:
            self.zpe = 0

    def U(self, T, xc="pbesol", units="eV"):
        """
        Calculates the internal energy of one formula unit of ideal gas.

        Example:

            U = S2_gas.U(300,xc="pbesol",units="eV")

        Args:

            T (float/ndarray): 1D Numpy array containing temperature data (in Kelvin) as floats, or a single temperature as a float.
            xc (str, optional): DFT XC functional used to calculate the ground state energy
            units (str, optional):  specifies the units as "eV", "J" (J/mol) or "kJ" (kJ/mol)

        Returns:

            U (float/ndarray): 1D Numpy array (with the same dimensions as T) containing the internal energies of one formula unit of gas, or a single internal energy float when a single temperature is passed as an argument.

        """
        U_func = interpolate.get_potential_nist_table(self.thermo_file, "U")
        E_dft = self.energies[xc]
        U_eV = (
            E_dft
            + self.zpe
            + U_func(T)
            * constants.physical_constants["joule-electron volt relationship"][0]
            / constants.N_A
        )

        if units == "eV":
            return U_eV

        elif units == "J":
            return (
                U_eV
                * constants.physical_constants["electron volt-joule relationship"][0]
                * constants.N_A
            )

        elif units == "kJ":
            return (
                U_eV
                * constants.physical_constants["electron volt-joule relationship"][0]
                * constants.N_A
                * 0.001
            )

    def H(self, T, xc="pbesol", units="eV"):
        """
        Calculates the Enthalpy of one formula unit of ideal gas.

        Examples:

            H = S2_gas.H(300,xc="pbesol",units="eV")
            H = S2_gas.H(np.linspace(100,700,1000))

        Args:

            T (float/ndarray): 1D Numpy array containing temperature data (in Kelvin) as floats, or a single temperature as a float.
            xc (str, optional): DFT XC functional used to calculate the ground state energy
            units (str, optional):  specifies the units as "eV", "J" (J/mol) or "kJ" (kJ/mol)

        Returns:

            H (float/ndarray):  1D Numpy array (with the same dimensions as T) containing the enthalpy of one formula unit of gas, or a single enthalpy float when a single temperature is passed as an argument.
        """
        H_func = interpolate.get_potential_nist_table(self.thermo_file, "H")
        E_dft = self.energies[xc]
        H_eV = (
            E_dft
            + self.zpe
            + H_func(T)
            * constants.physical_constants["joule-electron volt relationship"][0]
            / constants.N_A
        )

        if units == "eV":
            return H_eV

        elif units == "J":
            return (
                H_eV
                * constants.physical_constants["electron volt-joule relationship"][0]
                * constants.N_A
            )

        elif units == "kJ":
            return (
                H_eV
                * constants.physical_constants["electron volt-joule relationship"][0]
                * constants.N_A
                * 0.001
            )

    def mu(self, T, P, xc="pbesol", units="eV"):
        """
        Calculates the Gibbs Free Energy of one formula unit of ideal gas.

        Examples:

            mu = S2_gas.mu(300,xc="pbesol",units="eV")
            mu = S2_gas.mu(np.linspace(100,700,1000))

        Args:

            T (float/ndarray): 1D Numpy array containing temperature data (in Kelvin) as floats, or a single temperature as a float.
            P (float/ndarray): 2D Numpy array with a single row containing pressure data (in Pa) as floats, or a single pressure as a float.
            xc (str, optional): DFT XC functional used to calculate the ground state energy
            units (str, optional):  specifies the units as "eV", "J" (J/mol) or "kJ" (kJ/mol)

        Note:

            T, P are orthogonal 2D arrays of length m and n, populated in one row/column: in this case mu is an m x n matrix.
            Other T, P arrays will result in undefined behaviour.

        Returns:

            mu (float/ndarray): Gibbs Free Energy of one formula unit of ideal gas expressed as floats in a m x n Numpy array where T, P are orthogonal 2D arrays of length m and n
        """

        S_func = interpolate.get_potential_nist_table(self.thermo_file, "S")
        S = (
            S_func(T)
            * constants.physical_constants["joule-electron volt relationship"][0]
            / constants.N_A
        )
        H = self.H(T, xc=xc)
        mu_eV = (
            H
            - T * S
            + constants.physical_constants["Boltzmann constant in eV/K"][0]
            * T
            * np.log(P / 1e5)
        )

        if units == "eV":
            return mu_eV

        elif units == "J":
            return (
                mu_eV
                * constants.physical_constants["electron volt-joule relationship"][0]
                * constants.N_A
            )

        elif units == "kJ":
            return (
                mu_eV
                * constants.physical_constants["electron volt-joule relationship"][0]
                * constants.N_A
                * 0.001
            )
