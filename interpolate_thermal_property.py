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

from scipy.interpolate import interp1d
from numpy import genfromtxt

def get_potential_aims(file,property):
    """Thermodynamic property interpolation function. Requires phonopy-FHI-aims output file.
    Cv in kB/cell. All other properties in eV/cell
    """
    data = genfromtxt(file)
    T = data[:,0]
    if property in ('Cv','Cp','heat_capacity','C'):
        potential = data[:,3]
    elif property in ('U','internal_energy'):
        potential = data[:,2]
    elif property in ('F','A','Helmholtz','free_energy'):
        potential = data[:,1]
    elif property in ('TS'):
        potential = -data[:,4]
    elif property in ('S','Entropy','entropy'):
        potential = -data[:,4]/T
    else:
        raise RuntimeError('Property not found')        

    thefunction = interp1d(T,potential,kind='linear')

    return thefunction

def get_potential_nist_table(file, property):
    """Thermodynamic property interpolation function. Requires NIST-JANAF table. All properties in J, mol and K"""
    data = genfromtxt(file,skip_header=2)
    T = data[:,0]
    if property in ('Cp','C','heat_capacity'):
        potential = data[:,1]
    elif property in ('S','entropy'):
        potential = data[:,2]
    elif property in ('H','enthalpy'): 
        potential = (data[:,4] - data[0,4])*1E3
    elif property in ('U','internal_energy'):  
        # U = H - PV; for ideal gas molar PV = RT so U = H - RT
        from scipy.constants import R as R
        potential = (data[:,4] - data[0,4])*1E3 - R*data[:,0]
    elif property in ('DH','Delta_H','standard_enthalpy_change'):
        potential = data[:,4]*1E3
    else:
        raise RuntimeError('Property not found')        

    thefunction = interp1d(T,potential,kind='cubic')
    return thefunction

