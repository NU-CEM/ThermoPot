################################################################################
#  Copyright Adam J. Jackson (2014)                                            #
#                                                                              #
#   This program is free software: you can redistribute it and/or modify       #
#   it under the terms of the GNU General Public License as published by       #
#   the Free Software Foundation, either version 3 of the License, or          #
#   (at your option) any later version.                                        #
#                                                                              #
#   This program is distributed in the hope that it will be useful,            #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#   GNU General Public License for more details.                               #
#                                                                              #
#   You should have received a copy of the GNU General Public License          #
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.      #
################################################################################

import numpy as np
from scipy import constants
from interpolate_thermal_property import get_potential_aims, get_potential_nist_table

import os  # get correct path for datafiles when called from another directory
materials_directory = os.path.dirname(__file__)
# Append a trailing slash to make coherent directory name - this would select the
#  root directory in the case of no prefix, so we need to check
if materials_directory:
    materials_directory = materials_directory + '/'

class material(object):
    """Parent class for materials properties"""
    def __init__(self,name,pbesol_energy_eV=False, N=1):
        self.name = name
        self.pbesol_energy_eV = pbesol_energy_eV
        self.N = N

class solid(material):
    """
    Class for solid material data. 

    Sets properties:
    -------------------
    solid.name             (Identifying string)
    solid.pbesol_energy_eV (DFT total energy in eV with PBEsol XC functional)
    solid.fu_cell          (Number of formula units in periodic unit cell)
    solid.volume           (Volume of unit cell in cubic angstroms (m3 * 10^30))
    solid.phonons          (String containing path to phonopy-FHI-aims output data file)
    solid.N                (Number of atoms per formula unit)

    Sets methods:
    -------------------
    solid.U_eV(T), solid.U_J(T), solid.U_kJ(T) : Internal energy 
    solid.H_eV(T,P), solid.H_J(T,P), solid.H_kJ(T,P) : Enthalpy H = U + PV
    solid.mu_eV(T,P), solid.mu_J(T,P), solid.mu_kJ(T,P) : Chemical potential mu = U + PV - TS

    The material is assumed to be incompressible and without thermal expansion
    """
    def __init__(self, name, pbesol_energy_eV, fu_cell, volume, phonons=False, N=1):
        material.__init__(self,name,pbesol_energy_eV,N)

        self.fu_cell = fu_cell
        self.volume = volume
        self.phonons = materials_directory + phonons

    def U_eV(self,T,*P):
        """Internal energy of one formula unit of solid, expressed in eV.
        U = solid.U_eV(T)
        Returns a matrix with the same dimensions as T
        """
        U_func = get_potential_aims(self.phonons,'U')
        return (self.pbesol_energy_eV + U_func(T))/self.fu_cell

    def U_J(self,T):
        """Internal energy of one gram-mole of solid, expressed in J/mol
        U = solid.U_J(T)
        Returns a matrix with the same dimensions as T
        """
        return self.U_eV(T) * constants.physical_constants['electron volt-joule relationship'][0] * constants.N_A

    def U_kJ(self,T):
        """Internal energy of one gram-mole of solid, expressed in kJ/mol
        U = solid.U_kJ(T)
        Returns a matrix with the same dimensions as T
        """
        return self.U_J(T)/1000.

    def H_eV(self,T,P):
        """
        Enthalpy of one formula unit of solid, expressed in eV
        H = solid.H_eV(T,P)
 
        T, P may be orthogonal 2D arrays of length m and n, populated in one row/column:
        in this case H is an m x n matrix.

        T, P may instead be equal-length non-orthogonal 1D arrays, in which case H is a vector
        of H values corresponding to T,P pairs.

        Other T, P arrays may result in undefined behaviour.
        """
        U_func = get_potential_aims(self.phonons,'U')
        PV = P * self.volume * 1E-30 * constants.physical_constants['joule-electron volt relationship'][0] / constants.N_A
        return ((self.pbesol_energy_eV + U_func(T)) + PV)/self.fu_cell

    def H_J(self,T,P):
        """Enthalpy of one gram-mole of solid, expressed in J/mol
        H = solid.H_J(T,P)

        T, P may be orthogonal 2D arrays of length m and n, populated in one row/column:
        in this case H is an m x n matrix.

        T, P may instead be equal-length non-orthogonal 1D arrays, in which case H is a vector
        of H values corresponding to T,P pairs.

        Other T, P arrays may result in undefined behaviour.
        """
        return self.H_eV(T,P) * constants.physical_constants['electron volt-joule relationship'][0] * constants.N_A
    
    def H_kJ(self,T,P):
        """Enthalpy of one gram-mole of solid, expressed in kJ/mol
        H = solid.H_kJ(T,P)

        T, P may be orthogonal 2D arrays of length m and n, populated in one row/column:
        in this case H is an m x n matrix.

        T, P may instead be equal-length non-orthogonal 1D arrays, in which case H is a vector
        of H values corresponding to T,P pairs.

        Other T, P arrays may result in undefined behaviour.
        """
        return self.H_J(T,P) * 0.001

    def mu_eV(self,T,P):
        """
        Free energy of one formula unit of solid, expressed in eV
        mu = solid.mu_eV(T,P)
        T, P may be orthogonal 2D arrays of length m and n, populated in one row/column:
        in this case H is an m x n matrix.

        T, P may instead be equal-length non-orthogonal 1D arrays, in which case H is a vector
        of H values corresponding to T,P pairs.

        Other T, P arrays may result in undefined behaviour.
        """
        TS_func = get_potential_aims(self.phonons,'TS')
        H = self.H_eV(T,P)
        return H - TS_func(T)/self.fu_cell

    def mu_J(self,T,P):
        """
        Free energy of one mol of solid, expressed in J/mol
        mu = solid.mu_J(T,P)
        T, P may be orthogonal 2D arrays of length m and n, populated in one row/column:
        in this case H is an m x n matrix.

        T, P may instead be equal-length non-orthogonal 1D arrays, in which case H is a vector
        of H values corresponding to T,P pairs.

        Other T, P arrays may result in undefined behaviour.
        """
        return self.mu_eV(T,P) * constants.physical_constants['electron volt-joule relationship'][0] * constants.N_A

    def mu_kJ(self,T,P):
        """
        Free energy of one mol of solid, expressed in kJ/mol
        mu = solid.mu_kJ(T,P)
        T, P may be orthogonal 2D arrays of length m and n, populated in one row/column:
        in this case H is an m x n matrix.

        T, P may instead be equal-length non-orthogonal 1D arrays, in which case H is a vector
        of H values corresponding to T,P pairs.

        Other T, P arrays may result in undefined behaviour.
        """
        return self.mu_J(T,P) * 0.001

    def Cv_kB(self,T):
        """
        Constant-volume heat capacity of one formula unit of solid, expressed in units
        of the Boltzmann constant kB:
        Cv = solid.Cv_kB(T)
        T may be an array, in which case Cv will be an array of the same dimensions.
        """
        Cv_func = get_potential_aims(self.phonons,'Cv')
        return Cv_func(T)/self.fu_cell

    def Cv_eV(self,T):
        """
        Constant-volume heat capacity of one formula unit of solid, expressed in units
        of the Boltzmann constant kB:
        Cv = solid.Cv_eV(T)
        T may be an array, in which case Cv will be an array of the same dimensions.
        """
        return self.Cv_kB(T) / constants.physical_constants['Boltzmann constant in eV/K']
    
    def Cv_J(self,T):
        """
        Constant-volume heat capacity of solid, expressed in J/mol.
        Cv = solid.Cv_J(T)
        T may be an array, in which case Cv will be an array of the same dimensions.
        """
        return (self.Cv_kB(T) * constants.N_A  / constants.physical_constants['Boltzmann constant'][0] )

    def Cv_kJ(self,T):
       """
        Constant-volume heat capacity of solid, expressed in kJ/mol.
        Cv = solid.Cv_kJ(T)
        T may be an array, in which case Cv will be an array of the same dimensions.
        """
       return self.Cv_J(T) * 0.001
    
class ideal_gas(material):
    """
    Class for ideal gas properties. 

    Sets properties:
    -------------------
    ideal_gas.name             (string)
    ideal_gas.pbesol_energy_eV (DFT total energy in eV with PBEsol XC functional)
    ideal_gas.thermo_data      (String containing path to aims.vibrations output data file)
    ideal_gas.N                (Number of atoms per formula unit)

    Sets methods:
    -------------------
    ideal_gas.U_eV(T), ideal_gas.U_J(T), ideal_gas.U_kJ(T) : Internal energy 
    ideal_gas.H_eV(T), ideal_gas.H_J(T), ideal_gas.H_kJ(T) : Enthalpy H = U + PV
    ideal_gas.mu_eV(T,P), ideal_gas.mu_J(T,P), ideal_gas.mu_kJ(T,P) : Chemical potential mu = U + PV - TS

    Ideal gas law PV=nRT is applied: specifically (dH/dP) at const. T = 0 and int(mu)^P2_P1 dP = kTln(P2/P1)
    Enthalpy has no P dependence as volume is not restricted / expansion step is defined as isothermal
    """

    def __init__(self,name,pbesol_energy_eV,thermo_file,zpe_pbesol=0,zpe_lit=0,N=1):
        material.__init__(self, name, pbesol_energy_eV,N)
        self.thermo_file = materials_directory + thermo_file
        # Initialise ZPE to PBEsol value if provided. 
        # This looks redundant at the moment: the intent is to implement
        # some kind of switch or heirarchy of methods further down the line.
        if zpe_pbesol > 0:
            self.zpe = zpe_pbesol
        elif zpe_lit > 0:
            self.zpe = zpe_lit
        else:
            self.zpe = 0

    def U_eV(self,T):
        """Internal energy of one formula unit of ideal gas, expressed in eV.
        U = ideal_gas.U_eV(T)
        Returns a matrix with the same dimensions as T
        """
        U_func = get_potential_nist_table(self.thermo_file,'U')
        return (self.pbesol_energy_eV + self.zpe +
                U_func(T)*constants.physical_constants['joule-electron volt relationship'][0]/constants.N_A
                )
    def U_J(self,T):
        """Internal energy of one gram-mole of ideal gas, expressed in J/mol
        U = ideal_gas.U_J(T)
        Returns a matrix with the same dimensions as T
        """
        return self.U_eV(T) * constants.physical_constants['electron volt-joule relationship'][0] * constants.N_A

    def U_kJ(self,T):
        """Internal energy of one gram-mole of ideal gas, expressed in kJ/mol
        U = ideal_gas.U_kJ(T)
        Returns a matrix with the same dimensions as T
        """
        return self.U_J(T) * 0.001

    def H_eV(self,T,*P):
        """Enthalpy of one formula unit of ideal gas, expressed in eV
        H = ideal_gas.H_eV(T)
        Returns an array with the same dimensions as T

        Accepts ideal_gas.H_eV(T,P): P is unused
        """
        H_func = get_potential_nist_table(self.thermo_file,'H')
        return (self.pbesol_energy_eV + self.zpe +
                H_func(T)*constants.physical_constants['joule-electron volt relationship'][0]/constants.N_A
                )

    def H_J(self,T,*P):
        """Enthalpy of one gram-mole of ideal gas, expressed in J/mol
        H = ideal_gas.H_J(T)
        Returns an array with the same dimensions as T

        Accepts ideal_gas.H_eV(T,P): P is unused
        """
        return self.H_eV(T) * constants.physical_constants['electron volt-joule relationship'][0] * constants.N_A
    
    def H_kJ(self,T,*P):
        """Enthalpy of one gram-mole of ideal gas, expressed in kJ/mol
        H = ideal_gas.H_kJ(T,P)
        Returns an array with the same dimensions as T

        Accepts ideal_gas.H_eV(T,P): P is unused
        """
        return self.H_J(T) * 0.001

    def mu_eV(self,T,P):
        """
        Free energy of one formula unit of ideal gas, expressed in eV
        mu = ideal_gas.mu_eV(T,P)
        T, P may be orthogonal 2D arrays of length m and n, populated in one row/column:
        in this case H is an m x n matrix.

        T, P may instead be equal-length non-orthogonal 1D arrays, in which case H is a vector
        of H values corresponding to T,P pairs.

        Other T, P arrays may result in undefined behaviour.
        """
        S_func = get_potential_nist_table(self.thermo_file,'S')
        S = S_func(T) * constants.physical_constants['joule-electron volt relationship'][0]/constants.N_A
        H = self.H_eV(T)
        return H - T*S + constants.physical_constants['Boltzmann constant in eV/K'][0] * T * np.log(P/1E5)

    def mu_J(self,T,P):
        """
        Free energy of one mol of ideal gas, expressed in J/mol
        mu = ideal_gas.mu_J(T,P)
        T, P may be orthogonal 2D arrays of length m and n, populated in one row/column:
        in this case H is an m x n matrix.

        T, P may instead be equal-length non-orthogonal 1D arrays, in which case H is a vector
        of H values corresponding to T,P pairs.

        Other T, P arrays may result in undefined behaviour.
        """
        return self.mu_eV(T,P) * constants.physical_constants['electron volt-joule relationship'][0] * constants.N_A

    def mu_kJ(self,T,P):
        """
        Free energy of one mol of ideal gas, expressed in kJ/mol
        mu = ideal_gas.mu_kJ(T,P)
        T, P may be orthogonal 2D arrays of length m and n, populated in one row/column:
        in this case H is an m x n matrix.

        T, P may instead be equal-length non-orthogonal 1D arrays, in which case H is a vector
        of H values corresponding to T,P pairs.

        Other T, P arrays may result in undefined behaviour.
        """
        return self.mu_J(T,P) * 0.001

      
CZTS_kesterite = solid(name='Kesterite CZTS',
                       pbesol_energy_eV=-0.706480597450521e06,
                       fu_cell=2,
                       volume=310.86645888987351,
                       phonons='phonopy_output/czts.dat',
                       N=8
                      )

CZTS = CZTS_kesterite

CZTS_stannite = solid(name='Stannite CZTS',
                      pbesol_energy_eV=-0.353240264472923e06 ,
                      fu_cell=1,
                      volume=155.572938002,
                      phonons='phonopy_output/czts_stannite.dat',
                      N=8
                      )

Cu = solid(name='Cu',
           pbesol_energy_eV=-0.180838109862865e06,
           fu_cell=4,
           volume=45.38855878494433,
           phonons='phonopy_output/Cu.dat'
)

Sn = solid(name='Beta Sn',
           pbesol_energy_eV=-0.340581355063278e06, 
           fu_cell=2,
           volume=69.6092979612,
           phonons='phonopy_output/Sn.dat'
)

Zn = solid(name='Zn',
           pbesol_energy_eV=-0.981596036898606e05, 
           fu_cell=2,
           volume=28.256317,
           phonons='phonopy_output/Zn.dat'
)

alpha_S=solid(
    name='Alpha S',
    pbesol_energy_eV=-0.347575504588933e06,
    fu_cell=32,
    volume= 832.91786077871541,
    phonons='phonopy_output/alpha_S.dat'
)

Cu2S_low=solid(
    name='Low Cu2S',
    pbesol_energy_eV=-0.486150076546942e07,
    fu_cell=48,
    volume=2055.8786918601486,
    phonons='phonopy_output/Cu2S_low.dat',
    N=3
)

Cu2S=Cu2S_low

SnS2=solid(
    name='SnS2',
    pbesol_energy_eV=-0.192015393819437e06,
    fu_cell=1,
    volume=69.55551898436056,
    phonons='phonopy_output/SnS2.dat',
    N=3
)

SnS_pcma=solid(
    name='SnS',
    pbesol_energy_eV=-0.724613674134358E+06,
    fu_cell=4,
    volume=186.605514927,
    phonons='phonopy_output/SnS_pcma.dat',
    N=2
)

SnS=SnS_pcma

ZnS_wurtzite=solid(
    name='ZnS (wurtzite)',
    pbesol_energy_eV=-119886.323698657,
    fu_cell=2,
    volume=76.9580344589,
    phonons='phonopy_output/ZnS_wurtzite.dat',
    N=2
)

ZnS_zincblende=solid(
    name='ZnS (zinc blende)',
    pbesol_energy_eV=-59943.163599041,
    fu_cell=1,
    volume=38.4544005985,
    phonons='phonopy_output/ZnS_zincblende.dat',
    N=2
)

ZnS=ZnS_zincblende

Cu2SnS3_mo1=solid(
    name='Cu2SnS3 (Mo-1)',
    pbesol_energy_eV=-0.117318818763261e07,
    fu_cell=4,
    volume=469.83485422571772,
    phonons='phonopy_output/Cu2SnS3-mo1.dat',
    N=6
)

Cu2SnS3_mo2=solid(
    name='Cu2SnS3 (Mo-2)',
    pbesol_energy_eV=-0.293297062672424e06,
    fu_cell=1,
    volume=117.43775202426687,
    phonons='phonopy_output/Cu2SnS3-mo2.dat',
    N=6
)

S8=ideal_gas(
    name='S8',
    pbesol_energy_eV=-0.868936310037924e05,
    thermo_file='nist_janaf/S8.dat',
    zpe_pbesol=0.32891037,
    N=8
)

S2=ideal_gas(
    name='S2',
    pbesol_energy_eV=-0.217220682510473e05,
    thermo_file='nist_janaf/S2.dat',
    zpe_pbesol=0.04421415,
    N=2
)

def volume_calc(filename):
    """Calculate unit cell volume in cubic angstroms from FHI-aims geometry.in file"""
    import numpy as np
    lattice_vectors = []
    with open(filename, 'r') as f:
        for line in f:
            if line.split()[0] == 'lattice_vector':
                lattice_vectors.append(line.split()[1:4])

    lattice_vectors = np.array(lattice_vectors).astype(float)
    volume = np.dot(lattice_vectors[0],np.cross(lattice_vectors[1],lattice_vectors[2]))

    return abs(volume)

O2=ideal_gas(
    name='O2',
    pbesol_energy_eV=-0.408004839112704e04,
    thermo_file='nist_janaf/O2.dat',
    zpe_lit=0.0976, # Irikura, K. K. (2007). Journal of Physical and 
    #                 Chemical Reference Data, 36(2), 389-397.
    #                 doi:10.1063/1.2436891
    N=2
)
