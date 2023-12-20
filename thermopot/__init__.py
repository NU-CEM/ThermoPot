__version__ = "1.1.0-beta.1"

from scipy import constants

eV2Jmol = (
    constants.physical_constants["electron volt-joule relationship"][0] * constants.N_A
)

eV2kJmol = (
    constants.physical_constants["electron volt-joule relationship"][0]
    * constants.N_A
    / 1000
)

kB2JKmol = constants.physical_constants["Boltzmann constant"][0] * constants.N_A
