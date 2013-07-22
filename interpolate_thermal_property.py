from scipy.interpolate import interp1d
from numpy import genfromtxt

def get_potential_aims(file,property):
    """Thermodynamic property interpolation function. Requires phonopy-FHI-aims output file."""
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
    """Thermodynamic property interpolation function. Requires NIST-JANAF table."""
    data = genfromtxt(file,skip_header=2)
    T = data[:,0]
    if property in ('Cp','C','heat_capacity'):
        potential = data[:,1]
    elif property in ('S','entropy'):
        potential = data[:,2]
    #
    # Yes, this one looks strange. Both are normalised to U(0) = H(0) = 0 for easy addition to electronic energy.
    # For enthalpy change from standard conditions, use 'DH'
    elif property in ('H','enthalpy','U','internal_energy'): 
        potential = data[:,4] - data[0,4]
    elif property in ('DH','Delta_H','standard_enthalpy_change'):
        potential = data[:,4]
    else:
        raise RuntimeError('Property not found')        

    thefunction = interp1d(T,potential,kind='cubic')
    return thefunction

