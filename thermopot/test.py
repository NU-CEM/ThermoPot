import materials, calculations, interpolate

BaZrS3_calc = calculations.QHACalculation(
    filepath="/Users/w21013885/QHA_code/ThermoPot/BaZrS3/raw_aims_files/ternary/BaZrS3_Pnma/hse06/aims.out",
    ev_filepath="/Users/w21013885/QHA_code/ThermoPot/BaZrS3/qha_output/e-v.dat",
)

test = interpolate.get_potential_F_V(BaZrS3_calc.volumes, BaZrS3_calc.energies)
print(test)

test2 = interpolate.get_potential_V_T(
    "/Users/w21013885/BaZrS3_analysis/QHA/Ba2ZrS4_I4_mmm/volume-temperature.dat"
)
print(test2)
