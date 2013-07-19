from scipy.interpolate import interp1d
from numpy import genfromtxt
def get_potential_aims(file,property):
    """Thermodynamic property interpolation function. Requires phonopy-FHI-aims output file."""
    data = genfromtxt(file)
    T = data[:,0]
    if property in ('Cv','Cp','heat_capacity','C'):
        thefunction = interp1d(T,data[:,3])
    elif property in ('U','internal_energy'):
        thefunction = interp1d(T,data[:,2])
    elif property in ('F','A','Helmholtz','free_energy'):
        thefunction = interp1d(T,data[:,1])
    elif property in ('TS'):
        thefunction = interp1d(T,-data[:,4])
    elif property in ('S','Entropy','entropy'):
        thefunction = interp1d(T, -data[:,4]/T)
    else:
        raise RuntimeError('Property not found')        

    return thefunction
