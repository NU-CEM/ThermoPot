from thermopot import materials, calculations
import numpy as np
import yaml
from yaml import CLoader as Loader

from phonopy import PhonopyQHA

BaZrS3_calc = calculations.QHACalculation(
    filepath="/Users/w21013885/QHA_code/ThermoPot/BaZrS3/raw_aims_files/ternary/BaZrS3_Pnma/hse06/aims.out",
    ev_filepath="/Users/w21013885/QHA_code/ThermoPot/BaZrS3/qha_output/BaZrS3_Pnma/e-v.dat",
)

BaZrS3 = materials.Solid(
    "BaZrS3",
    {"Ba": 1, "Zr": 1, "S": 3},
    phonon_filepath="../BaZrS3/phonopy_output/BaZrS3_Pnma.dat",
    qha_calculation=BaZrS3_calc,
)

print(
    BaZrS3.mu(
        T=0,
        P=1,
        V_T_filepath="/Users/w21013885/QHA_code/ThermoPot/BaZrS3/qha_output/volume-temperature.dat",
        xc="hse06",
        units="eV",
    )
)