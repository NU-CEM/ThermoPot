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

        if calculation is not False:
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

    def U(self, T, xc="pbesol", units="eV"):
        """
        Calculates the internal energy of one formular unit of solid.
        
        Example:
        
            U = BaS.U(300,xc="pbesol",units="eV")
            
        Args:
        
            T (float/ndarray): 1D Numpy array containing temperature data (in Kelvin) as a float, or a single temperature as a float.
            xc (str, optional): DFT XC functional used to calculate the ground state energy
            units (str, optional):  specifies the units as "eV", "J" or "kJ"
        
        Returns:
        
            U (float/ndarray): 1D Numpy array (with the same dimensions as T) containing the internal energies of one formula unit of solid 
        
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
        Enthalpy of solid, H = U + PV.
        H = solid.H(T,P)

         Pressure is given in Pa.
        The units keyword specifies the units used.
        It defaults to "eV per formula unit". Other options are "J" (J/mol) and "kJ" (kJ/mol)

        The xc keyword specifies the DFT XC functional used to calculate the ground state energy.
        If not specified, it defaults to `pbesol`.

        T, P may be orthogonal 2D arrays of length m and n, populated in one row/column:
        in this case H is an m x n matrix.

        T, P may instead be equal-length non-orthogonal 1D arrays, in which case H is a vector
        of H values corresponding to T,P pairs.

        Other T, P arrays may result in undefined behaviour.
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
        Free energy of one formula unit of solid, mu = U + PV - TS
        mu = solid.mu(T,P)

        The units keyword specifies the units used.
        It defaults to "eV per formular unit". Other options are "J" (J/mol) and "kJ" (kJ/mol)


        The xc keyword specifies the DFT XC functional used to calculate the ground state energy.
        If not specified, it defaults to `pbesol`.

        T, P may be orthogonal 2D arrays of length m and n, populated in one row/column:
        in this case H is an m x n matrix.

        T, P may instead be equal-length non-orthogonal 1D arrays, in which case H is a vector
        of H values corresponding to T,P pairs.

        Other T, P arrays may result in undefined behaviour.
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
        Constant-volume heat capacity of one formula unit of solid
        Cv = solid.Cv(T)

        The "units" keyword specifies the units used.
        If unspecified, it defaults to units kB (Boltzmann constant).
        The other options are "eV" (eV/K), "J" (J/molK) and "kJ" (kJ/molK).

        T may be an array, in which case Cv will be an array of the same dimensions.
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

    Sets properties:
    -------------------
    ideal_gas.name             (string)
    ideal_gas.stoichiometry    (Dict relating element to number of atoms in a single formula unit)
    ideal_gas.energies         (Dict relating xc functional with DFT total energy)
    ideal_gas.thermo_data      (String containing path to NIST data table)
    ideal_gas.N                (Number of atoms per formula unit)

    Sets methods:
    -------------------
    ideal_gas.U(T) : Internal energy
    ideal_gas.H(T) : Enthalpy H = U + PV
    ideal_gas.mu(T,P) : Chemical potential mu = U + PV - TS


    Ideal gas law PV=nRT is applied: specifically (dH/dP) at const. T = 0 and int(mu)^P2_P1 dP = kTln(P2/P1)
    Enthalpy has no P dependence as volume is not restricted / expansion step is defined as isothermal

    Temperature is given in Kelvin, Pressure is given in Pa.
    """

    # TODO:script for calculating zpe values

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

        if calculation is not False:
            Material.__init__(
                self, name, stoichiometry, {calculation.xc: calculation.energy}
            )
        else:
            Material.__init__(self, name, stoichiometry, energies)
        self.thermo_file = materials_directory + thermo_file
        # Initialise ZPE to HSE06 value if provided.
        # This looks redundant at the moment: the intent is to implement
        # some kind of switch or heirarchy of methods further down the line.
        if zpe_hse06 > 0:
            self.zpe = zpe_pbesol
        elif zpe_pbesol > 0:
            self.zpe = zpe_pbesol
        elif zpe_lit > 0:
            self.zpe = zpe_lit
        else:
            self.zpe = 0

    def U(self, T, xc="pbesol", units="eV"):
        """Internal energy of one formula unit of ideal gas.
        U = ideal_gas.U_eV(T)

        The units keyword specifies the units used.
        If not specified, it  defaults to  "eV".
        Can also be specified as "J" for J/mol or "kJ" for kJ/mol"
        The xc keyword specifies the DFT XC functional used to calculate the ground state energy.
        If not specified, it defaults to `pbesol`.

        Returns a matrix with the same dimensions as T
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

    def H(self, T, *P, xc="pbesol", units="eV"):
        """Enthalpy of one formula unit of ideal gas
        H = ideal_gas.H(T)

        Default is expressed in eV. Can also be expressed in J/mol or kJ/mol.

        The xc keyword specifies the DFT XC functional used to calculate the ground state energy.
        If not specified, it defaults to `pbesol`.

        Returns an array with the same dimensions as T

        Accepts ideal_gas.H_eV(T,P): P is unused
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
        Free energy of one formula unit of ideal gas.
        mu = ideal_gas.mu(T,P)

        Default is expressed in eV.
        Can also been expressed in J/mol or kJ/mol.

        The xc keyword specifies the DFT XC functional used to calculate the ground state energy.
        If not specified, it defaults to `pbesol`.

        T, P may be orthogonal 2D arrays of length m and n, populated in one row/column:
        in this case H is an m x n matrix.

        T, P may instead be equal-length non-orthogonal 1D arrays, in which case H is a vector
        of H values corresponding to T,P pairs.

        Other T, P arrays may result in undefined behaviour.
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
