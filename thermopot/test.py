import calculations, materials
BaZrS3_calc = calculations.AimsCalculation("/Users/w21013885/QHA_code/ThermoPot/BaZrS3/raw_aims_files/ternary/BaZrS3_Pnma/hse06/aims.out")
#print(vars(BaZrS3_calc))
BaZrS3_calc = calculations.QHACalculation(filepath='/Users/w21013885/QHA_code/ThermoPot/BaZrS3/raw_aims_files/ternary/BaZrS3_Pnma/hse06/aims.out',ev_filepath='/Users/w21013885/QHA_code/ThermoPot/BaZrS3/qha_output/e-v.dat')
#print(vars(BaZrS3_calc))

BaZrS3 = materials.Solid("BaZrS3",{"Ba": 1,"Zr": 1,"S":3},phonon_filepath="/Users/w21013885/QHA_code/ThermoPot/BaZrS3/phonopy_output/BaZrS3_Pnma.dat",QHACalculation=BaZrS3_calc)
print(vars(BaZrS3))
print(BaZrS3.U)
