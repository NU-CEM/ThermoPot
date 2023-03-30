import materials, calculations, interpolate

BaZrS3_calc = calculations.QHACalculation(
    filepath="/Users/w21013885/QHA_code/ThermoPot/BaZrS3/raw_aims_files/ternary/BaZrS3_Pnma/hse06/aims.out",
    ev_filepath="/Users/w21013885/QHA_code/ThermoPot/BaZrS3/qha_output/e-v.dat")
#BaZrS3_calc = calculations.AimsCalculation(
#    "/Users/w21013885/QHA_code/ThermoPot/BaZrS3/raw_aims_files/ternary/BaZrS3_Pnma/hse06/aims.out")

test = interpolate.get_potential_F_V(BaZrS3_calc.volumes, BaZrS3_calc.energies)
print(test)

test2 = interpolate.get_potential_V_T("/Users/w21013885/BaZrS3_analysis/QHA/Ba2ZrS4_I4_mmm/volume-temperature.dat")
print(test2)

BaZrS3 = materials.Solid("BaZrS3", {"Ba": 1, "Zr": 1, "S": 3}, phonon_filepath="../BaZrS3/phonopy_output/BaZrS3_Pnma.dat", qha_calculation=BaZrS3_calc)
#BaZrS3 = materials.Solid("BaZrS3", {"Ba": 1, "Zr": 1, "S": 3}, phonon_filepath="../BaZrS3/phonopy_output/BaZrS3_Pnma.dat", calculation=BaZrS3_calc)
print(BaZrS3.U(T=298, xc="hse06", units="eV"))
print(BaZrS3.mu(T=0,P=1,V_T_filepath='/Users/w21013885/QHA_code/ThermoPot/BaZrS3/qha_output/volume-temperature.dat',xc="hse06", units="eV"))