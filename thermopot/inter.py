import numpy as np
from scipy.interpolate import interp1d, interp2d
from numpy import genfromtxt
import re
from thermopot import *
import interpolate
import matplotlib.pyplot as plt


def get_potential_aims(file, property):
    """Thermodynamic property interpolation function. Requires phonopy-FHI-aims output file.
    Reads data for S and Cv expressed in J/K/mol, F and U in kJ/mol,
    TS in J/mol.
    Outputs data for S and Cv in kB/cell, U, F and TS in eV/cell.
    """
    data = genfromtxt(file)
    T = data[:, 0]
    if property in ("Cv", "Cp", "heat_capacity", "C"):
        potential = data[:, 3] / kB2JKmol
    elif property in ("U", "internal_energy"):
        potential = data[:, 4] / eV2kJmol
    elif property in ("F", "A", "Helmholtz", "free_energy"):
        potential = data[:, 1] / eV2kJmol
    elif property in ("S", "Entropy", "entropy"):
        potential = data[:, 2] / kB2JKmol
    elif property in ("TS"):
        potential = (T * data[:, 2]) / eV2Jmol
    else:
        raise RuntimeError("Property not found")
    thefunction = interp1d(T, potential, kind="linear")

    return thefunction


num = interpolate.get_potential_aims(
    "/Users/w21013885/QHA_code/ThermoPot/BaZrS3/phonopy_output/BaZrS3_Pnma.dat", "Cv"
)
print(num)
